from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

from cse6406.core.command_runner import CommandRunner
from cse6406.core.errors import PipelineError
from cse6406.core.toolchain import Toolchain
from cse6406.domain.design import ExperimentalDesign
from cse6406.gene_trees.iqtree import IQTreeInferer
from cse6406.reporting.csv_writer import CSVResultWriter
from cse6406.simulation.alisim import AliSimSimulator
from cse6406.simulation.simphy import SimPhySimulator
from cse6406.stage1.branch_scorer import BranchScorer
from cse6406.stage1.calibration import CalibrationAnalyzer
from cse6406.stage1.ils import ILSAnalyzer
from cse6406.stage2.astral import AstralRunner
from cse6406.stage2.contractor import SupportTreeContractor
from cse6406.stage2.hybrid_wastral import HybridWeightedAstralRunner
from cse6406.stage2.species_tree import SpeciesTreeEvaluator
from cse6406.stage2.support_preparer import SupportTreePreparer
from cse6406.stage2.wastral import WeightedAstralRunner
from cse6406.trees.newick import combine_newick, read_tree


class ExperimentPipeline:
    """Orchestrate the full design while preserving every reusable artifact."""

    def __init__(self, design: ExperimentalDesign, tools: Toolchain):
        design.validate()
        self.design, self.tools = design, tools
        self.runner = CommandRunner()
        self.writer = CSVResultWriter()

    def preflight(self, *, stage2: bool = True) -> dict[str, dict[str, str]]:
        required = ["simphy", "iqtree"] + (["astral", "wastral"] if stage2 else [])
        self.tools.require(*required)
        metadata = {name: self.tools.metadata(name) for name in required}
        self.design.output.mkdir(parents=True, exist_ok=True)
        self._validate_or_write_design()
        (self.design.output / "preflight.json").write_text(json.dumps(metadata, indent=2) + "\n")
        return metadata

    def run(self, *, stage2: bool = True) -> None:
        self.preflight(stage2=stage2)
        assert self.tools.simphy and self.tools.iqtree
        simulator = SimPhySimulator(self.tools.simphy, self.runner)
        alisim = AliSimSimulator(self.tools.iqtree, self.runner)
        inferer = IQTreeInferer(self.tools.iqtree, self.runner)
        scorer, calibration, ils_analyzer = BranchScorer(), CalibrationAnalyzer(), ILSAnalyzer()
        all_observations, ils_rows, stage2_rows = [], [], []
        histories = {}

        # Validate the intended biological treatments before starting the much
        # more expensive alignment and gene-tree inference stages.
        for level in self.design.ils_levels:
            for replicate in range(1, self.design.replicates + 1):
                species_tree, true_gene_trees = simulator.run(self.design, level.name, replicate)
                if len(true_gene_trees) != self.design.loci:
                    raise PipelineError(
                        f"Expected {self.design.loci} true gene trees; found {len(true_gene_trees)}"
                    )
                histories[(level.name, replicate)] = (species_tree, true_gene_trees)
                ils_rows.append(ils_analyzer.analyze(
                    ils=level.name, replicate=replicate, species_tree=species_tree,
                    gene_trees=true_gene_trees, population_size=level.population_size,
                ))

        results = self.design.output / "results"
        ils_summary = ils_analyzer.summarize(
            ils_rows,
            low_max=self.design.low_ils_max_nrf,
            high_min=self.design.high_ils_min_nrf,
            minimum_gap=self.design.minimum_ils_nrf_gap,
        )
        self.writer.write(results / "ils_verification.csv", ils_rows)
        self.writer.write(results / "ils_summary.csv", ils_summary)
        if self.design.enforce_ils_gate:
            ils_analyzer.require_acceptable(ils_summary)

        for level in self.design.ils_levels:
            for replicate in range(1, self.design.replicates + 1):
                species_tree, true_gene_trees = histories[(level.name, replicate)]
                for length in self.design.sequence_lengths:
                    alignments = alisim.run(
                        self.design, level.name, replicate, length, true_gene_trees
                    )
                    for model in self.design.models:
                        condition = next(c for c in self.design.conditions()
                                         if c.ils == level.name and c.replicate == replicate
                                         and c.sequence_length == length and c.model == model)
                        estimates = inferer.run(self.design, condition, alignments)
                        for locus, (estimated, truth) in enumerate(
                            zip(estimates, true_gene_trees, strict=True), 1
                        ):
                            all_observations.extend(scorer.score(
                                condition, locus, read_tree(estimated), read_tree(truth)
                            ))
                        if stage2:
                            stage2_rows.append(self._run_stage2(
                                condition, species_tree, estimates
                            ))

        branches, calibration_rows, summary, replicate_summary = calibration.analyze(
            all_observations
        )
        self.writer.write(results / "stage1_branches.csv", branches)
        self.writer.write(results / "stage1_calibration.csv", calibration_rows)
        self.writer.write(results / "stage1_summary.csv", summary)
        self.writer.write(results / "stage1_replicate_summary.csv", replicate_summary)
        self.writer.write(results / "stage2_species_tree_error.csv", stage2_rows)

    def _validate_or_write_design(self) -> None:
        path = self.design.output / "design.json"
        requested = asdict(self.design)
        requested["output"] = str(requested["output"])
        requested = json.loads(json.dumps(requested))
        if path.is_file():
            existing = json.loads(path.read_text())
            if existing != requested:
                raise PipelineError(
                    f"Output directory belongs to a different design: {path}. "
                    "Use the original arguments or select a new --output directory."
                )
        else:
            path.write_text(json.dumps(requested, indent=2) + "\n")

    def _run_stage2(self, condition, species_tree: Path, estimates: list[Path]):
        assert self.tools.astral and self.tools.wastral
        output = self.design.model_dir(condition) / "species_trees"
        input_trees = output / "estimated_gene_trees.tre"
        combine_newick(estimates, input_trees)
        astral_tree = AstralRunner(self.tools.astral, self.runner).run(
            input_trees, output / "astral_unweighted.tre", threads=self.design.threads
        )
        threshold_label = f"{self.design.contracted_abayes_threshold:.2f}"
        contracted_input = output / f"estimated_gene_trees_contracted_abayes_{threshold_label}.tre"
        SupportTreeContractor().contract(
            input_trees, contracted_input,
            threshold=self.design.contracted_abayes_threshold,
        )
        contracted_astral_tree = AstralRunner(self.tools.astral, self.runner).run(
            contracted_input, output / f"astral_contracted_abayes_{threshold_label}.tre",
            threads=self.design.threads,
        )
        weighted_input = output / "estimated_gene_trees_wastral_ready.tre"
        SupportTreePreparer().prepare(input_trees, weighted_input)
        wastral_tree = WeightedAstralRunner(self.tools.wastral, self.runner).run(
            weighted_input, output / "wastral_support.tre", threads=self.design.threads
        )
        hybrid_wastral_tree = HybridWeightedAstralRunner(
            self.tools.wastral, self.runner
        ).run(
            weighted_input, output / "wastral_hybrid.tre", threads=self.design.threads
        )
        return SpeciesTreeEvaluator().evaluate(
            condition=condition, true_species_tree=species_tree,
            astral_tree=astral_tree, contracted_astral_tree=contracted_astral_tree,
            wastral_tree=wastral_tree, hybrid_wastral_tree=hybrid_wastral_tree,
        )
