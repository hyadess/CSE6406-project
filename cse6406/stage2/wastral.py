from pathlib import Path

from cse6406.core.command_runner import CommandRunner
from cse6406.core.errors import PipelineError
from cse6406.stage2.threads import aster_threads
from cse6406.trees.newick import informative_splits, read_tree_text


class WeightedAstralRunner:
    """Run wASTRAL support weighting (mode 2) with local Bayesian scaling."""

    def __init__(self, binary: Path, runner: CommandRunner):
        self.binary, self.runner = binary, runner

    def run(self, input_trees: Path, output_tree: Path, *, threads: str = "1") -> Path:
        self._validate_support(input_trees)
        if not output_tree.is_file():
            self.runner.run(
                [self.binary, "--mode", 2, "-B", "-i", input_trees.resolve(),
                 "-o", output_tree.resolve(), "-t", aster_threads(threads), "-u", 0],
                log=output_tree.with_suffix(".log"),
            )
        return output_tree

    @staticmethod
    def _validate_support(input_trees: Path) -> None:
        trees = [line.strip() for line in input_trees.read_text().splitlines() if line.strip()]
        if not trees:
            raise PipelineError(f"No input gene trees: {input_trees}")
        missing = 0
        outside = []
        for text in trees:
            for _split, support in informative_splits(read_tree_text(text), with_support=True):
                if support is None:
                    missing += 1
                elif not 0.0 <= support <= 1.0:
                    outside.append(support)
        if missing or outside:
            raise PipelineError(
                f"wASTRAL input failed support validation: {missing} missing labels, "
                f"{len(outside)} outside [0,1]"
            )
