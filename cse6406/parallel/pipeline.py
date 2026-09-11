from __future__ import annotations

import csv
import os
from dataclasses import asdict, fields
from pathlib import Path

from cse6406.core.errors import PipelineError
from cse6406.core.toolchain import Toolchain
from cse6406.domain.design import ExperimentalDesign
from cse6406.parallel import tasks
from cse6406.parallel.executor import WorkerPool
from cse6406.pipeline import ExperimentPipeline
from cse6406.stage1.branch_scorer import BranchObservation
from cse6406.stage1.calibration import CalibrationAnalyzer
from cse6406.stage1.ils import ILSAnalyzer


def default_workers() -> int:
    """Concurrent tasks to run when the caller does not choose a number."""
    available = getattr(os, "process_cpu_count", None) or os.cpu_count
    return max(1, (available() or 1))


class ParallelExperimentPipeline(ExperimentPipeline):
    """Run the sequential design across a pool of worker processes.

    The phase order, the seeds, the output layout, the resume rules, and the row
    order of every result table are identical to :class:`ExperimentPipeline`.
    The only difference is that the independent units inside each phase - one
    SimPhy replicate, one simulated locus, one gene-tree inference, one
    condition's species trees - run concurrently instead of one at a time.
    """

    def __init__(
        self, design: ExperimentalDesign, tools: Toolchain, *,
        workers: int | None = None, tool_threads: str = "1",
    ):
        super().__init__(design, tools)
        self.workers = default_workers() if workers is None else max(1, int(workers))
        self.tool_threads = str(tool_threads)

    def run(self, *, stage2: bool = True) -> None:
        self.preflight(stage2=stage2)
        design = self.design
        results = design.output / "results"
        conditions = design.conditions()
        print(
            f"[parallel] {self.workers} worker process(es); "
            f"{self.tool_threads} thread(s) per external tool; "
            f"{len(conditions)} condition(s) x {design.loci} loci",
            flush=True,
        )

        with WorkerPool(
            self.workers,
            initializer=tasks.initialize,
            initargs=(self.tools, design, self.tool_threads),
        ) as pool:
            histories = self._simulate_histories(pool)
            self._verify_ils(pool, histories, results)
            alignments = self._simulate_alignments(pool, histories)
            estimates, observations = self._infer_gene_trees(
                pool, conditions, histories, alignments, results
            )
            del alignments
            stage2_rows = (
                self._build_species_trees(pool, conditions, histories, estimates)
                if stage2 else []
            )

        del estimates
        _, calibration_rows, summary, replicate_summary = CalibrationAnalyzer().analyze(
            observations, materialize_branches=False
        )
        del observations
        self.writer.write(results / "stage1_calibration.csv", calibration_rows)
        self.writer.write(results / "stage1_summary.csv", summary)
        self.writer.write(results / "stage1_replicate_summary.csv", replicate_summary)
        self.writer.write(results / "stage2_species_tree_error.csv", stage2_rows)

    # ---------------------------------------------------------------- phases

    def _simulate_histories(self, pool: WorkerPool) -> dict[tuple[str, int], tuple]:
        """Phase 1: true species and gene trees, one task per replicate."""
        design = self.design
        keys = [
            (level.name, replicate)
            for level in design.ils_levels
            for replicate in range(1, design.replicates + 1)
        ]
        histories = {}
        outputs = pool.map_ordered(tasks.simphy, keys, label="simphy")
        for key, (species_tree, true_gene_trees) in zip(keys, outputs, strict=True):
            if len(true_gene_trees) != design.loci:
                raise PipelineError(
                    f"Expected {design.loci} true gene trees; found {len(true_gene_trees)}"
                )
            histories[key] = (species_tree, true_gene_trees)
        return histories

    def _verify_ils(self, pool: WorkerPool, histories: dict, results: Path) -> None:
        """Phase 2: the preregistered ILS gate, before any expensive inference."""
        design = self.design
        analyzer = ILSAnalyzer()
        work = [
            (level.name, replicate, *histories[(level.name, replicate)], level.population_size)
            for level in design.ils_levels
            for replicate in range(1, design.replicates + 1)
        ]
        ils_rows = list(pool.map_ordered(tasks.verify_ils, work, label="ils-check"))
        ils_summary = analyzer.summarize(
            ils_rows,
            low_max=design.low_ils_max_nrf,
            high_min=design.high_ils_min_nrf,
            minimum_gap=design.minimum_ils_nrf_gap,
        )
        self.writer.write(results / "ils_verification.csv", ils_rows)
        self.writer.write(results / "ils_summary.csv", ils_summary)
        if design.enforce_ils_gate:
            analyzer.require_acceptable(ils_summary)

    def _simulate_alignments(
        self, pool: WorkerPool, histories: dict,
    ) -> dict[tuple[str, int, int], list[Path]]:
        """Phase 3: AliSim alignments, one task per locus."""
        design = self.design
        work, keys = [], []
        for level in design.ils_levels:
            for replicate in range(1, design.replicates + 1):
                _species_tree, true_gene_trees = histories[(level.name, replicate)]
                for length in design.sequence_lengths:
                    keys.append((level.name, replicate, length))
                    for locus, true_tree in enumerate(true_gene_trees, 1):
                        work.append((level.name, replicate, length, locus, true_tree))
        produced = pool.map_ordered(tasks.alisim_locus, work, label="alisim")
        alignments: dict[tuple[str, int, int], list[Path]] = {key: [] for key in keys}
        for (ils, replicate, length, _locus, _tree), path in zip(work, produced, strict=True):
            alignments[(ils, replicate, length)].append(path)
        return alignments

    def _infer_gene_trees(
        self, pool: WorkerPool, conditions: list, histories: dict,
        alignments: dict, results: Path,
    ) -> tuple[dict, list[BranchObservation]]:
        """Phase 4: IQ-TREE inference plus branch scoring, one task per locus."""
        design = self.design
        work = []
        for condition in conditions:
            _species_tree, true_gene_trees = histories[(condition.ils, condition.replicate)]
            locus_alignments = alignments[
                (condition.ils, condition.replicate, condition.sequence_length)
            ]
            for locus, (alignment, true_tree) in enumerate(
                zip(locus_alignments, true_gene_trees, strict=True), 1
            ):
                work.append((condition, locus, alignment, true_tree))

        estimates: dict[str, list[Path]] = {condition.key: [] for condition in conditions}
        observations: list[BranchObservation] = []
        produced = pool.map_ordered(tasks.gene_tree_locus, work, label="iqtree")
        # The branch table is by far the largest output, so it is streamed to
        # disk in submission order as workers finish instead of being held as a
        # second full copy in memory.
        branches_path = results / "stage1_branches.csv"
        branches_path.parent.mkdir(parents=True, exist_ok=True)
        header = [field.name for field in fields(BranchObservation)]
        with branches_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=header)
            writer.writeheader()
            for (condition, *_rest), (tree, scored) in zip(work, produced, strict=True):
                estimates[condition.key].append(tree)
                observations.extend(scored)
                writer.writerows(asdict(row) for row in scored)
        return estimates, observations

    def _build_species_trees(
        self, pool: WorkerPool, conditions: list, histories: dict, estimates: dict,
    ) -> list[dict[str, object]]:
        """Phase 5: ASTRAL and wASTRAL species trees, one task per condition."""
        work = [
            (
                condition,
                histories[(condition.ils, condition.replicate)][0],
                estimates[condition.key],
            )
            for condition in conditions
        ]
        return list(pool.map_ordered(tasks.species_trees, work, label="stage2"))
