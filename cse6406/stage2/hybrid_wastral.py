from pathlib import Path

from cse6406.core.command_runner import CommandRunner
from cse6406.stage2.threads import aster_threads
from cse6406.stage2.wastral import WeightedAstralRunner


class HybridWeightedAstralRunner:
    """Run default hybrid wASTRAL with explicit local-Bayesian scaling."""

    def __init__(self, binary: Path, runner: CommandRunner):
        self.binary, self.runner = binary, runner

    def run(self, input_trees: Path, output_tree: Path, *, threads: str = "1") -> Path:
        WeightedAstralRunner._validate_support(input_trees)
        if not output_tree.is_file():
            self.runner.run(
                [self.binary, "-B", "-i", input_trees.resolve(),
                 "-o", output_tree.resolve(), "-t", aster_threads(threads), "-u", 0],
                log=output_tree.with_suffix(".log"),
            )
        return output_tree
