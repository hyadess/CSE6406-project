from pathlib import Path

from cse6406.core.command_runner import CommandRunner


class AstralRunner:
    """Run unweighted ASTRAL-IV on estimated gene-tree topologies."""

    def __init__(self, binary: Path, runner: CommandRunner):
        self.binary, self.runner = binary, runner

    def run(self, input_trees: Path, output_tree: Path, *, threads: str = "1") -> Path:
        if not output_tree.is_file():
            self.runner.run(
                [self.binary, "-i", input_trees.resolve(), "-o", output_tree.resolve(),
                 "-t", threads, "-u", 0],
                log=output_tree.with_suffix(".log"),
            )
        return output_tree
