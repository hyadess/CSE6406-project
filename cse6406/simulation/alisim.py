from pathlib import Path

from cse6406.core.command_runner import CommandRunner
from cse6406.core.errors import PipelineError
from cse6406.domain.design import ExperimentalDesign


class AliSimSimulator:
    """Simulate the shared alignment used by all five analysis models."""

    def __init__(self, iqtree: Path, runner: CommandRunner):
        self.iqtree, self.runner = iqtree, runner

    def run_locus(
        self, design: ExperimentalDesign, ils: str, replicate: int,
        length: int, locus: int, true_tree: Path,
    ) -> Path:
        """Simulate exactly one locus; the unit of work shared with the parallel runner."""
        output = design.alignment_dir(ils, replicate, length)
        output.mkdir(parents=True, exist_ok=True)
        prefix = f"locus_{locus:04d}"
        alignment = output / f"{prefix}.phy"
        if not alignment.is_file():
            self.runner.run(
                [self.iqtree, "--alisim", prefix, "-t", true_tree.resolve(),
                 "-m", design.generating_model, "--length", length,
                 "--seed", design.seed + replicate * 100_000 + length + locus,
                 "--quiet"],
                cwd=output,
                log=design.replicate_dir(ils, replicate) / "logs"
                / f"alisim_L{length}_locus_{locus:04d}.log",
            )
        if not alignment.is_file() or alignment.stat().st_size == 0:
            raise PipelineError(f"AliSim did not create {alignment}")
        return alignment

    def run(
        self, design: ExperimentalDesign, ils: str, replicate: int,
        length: int, true_gene_trees: list[Path],
    ) -> list[Path]:
        output = design.alignment_dir(ils, replicate, length)
        output.mkdir(parents=True, exist_ok=True)
        return [
            self.run_locus(design, ils, replicate, length, locus, true_tree)
            for locus, true_tree in enumerate(true_gene_trees, 1)
        ]
