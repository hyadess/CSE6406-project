from pathlib import Path

from cse6406.core.command_runner import CommandRunner
from cse6406.core.errors import PipelineError
from cse6406.domain.condition import ExperimentCondition
from cse6406.domain.design import ExperimentalDesign


class IQTreeInferer:
    """Infer ML gene trees with aBayes support under one fixed model."""

    def __init__(self, binary: Path, runner: CommandRunner):
        self.binary, self.runner = binary, runner

    def run_locus(
        self, design: ExperimentalDesign, condition: ExperimentCondition,
        locus: int, alignment: Path, *, threads: str | None = None,
    ) -> Path:
        """Infer exactly one gene tree; the unit of work shared with the parallel runner."""
        output = design.model_dir(condition)
        output.mkdir(parents=True, exist_ok=True)
        prefix = output / f"locus_{locus:04d}"
        tree = prefix.with_suffix(".treefile")
        if not tree.is_file():
            self.runner.run(
                [self.binary, "-s", alignment.resolve(), "-m", condition.model.iqtree_model,
                 "--abayes", "-T", design.threads if threads is None else threads,
                 "--seed", design.seed + condition.replicate * 100_000 + locus,
                 "--prefix", prefix, "--redo", "--quiet"],
                log=design.replicate_dir(condition.ils, condition.replicate) / "logs"
                / f"iqtree_L{condition.sequence_length}_{condition.model.name}"
                f"_locus_{locus:04d}.log",
            )
        if not tree.is_file() or not tree.read_text().strip().endswith(";"):
            raise PipelineError(f"Invalid IQ-TREE output: {tree}")
        return tree

    def run(
        self, design: ExperimentalDesign, condition: ExperimentCondition,
        alignments: list[Path], *, threads: str | None = None,
    ) -> list[Path]:
        output = design.model_dir(condition)
        output.mkdir(parents=True, exist_ok=True)
        return [
            self.run_locus(design, condition, locus, alignment, threads=threads)
            for locus, alignment in enumerate(alignments, 1)
        ]
