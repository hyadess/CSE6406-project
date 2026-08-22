from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

from cse6406.core.command_runner import CommandRunner
from cse6406.core.errors import PipelineError
from cse6406.domain.design import ExperimentalDesign


def _decimal(value: float) -> str:
    return f"{value:.12f}".rstrip("0").rstrip(".")


class SimPhySimulator:
    """Generate true species and per-locus gene trees under the MSC."""

    def __init__(self, binary: Path, runner: CommandRunner):
        self.binary, self.runner = binary, runner

    def run(self, design: ExperimentalDesign, ils: str, replicate: int) -> tuple[Path, list[Path]]:
        destination = design.replicate_dir(ils, replicate) / "simphy"
        existing = self._outputs(destination)
        if existing:
            return existing
        if destination.exists() and any(destination.iterdir()):
            raise PipelineError(f"Partial SimPhy output exists: {destination}")

        level = design.ils_level(ils)
        destination.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="cse6406_simphy_") as temporary:
            generated = Path(temporary) / "output"
            command = [
                self.binary, "-rs", 1, "-rl", f"f:{design.loci}", "-rg", 1,
                "-sp", f"f:{level.population_size}",
                "-su", f"e:{_decimal(design.substitution_rate)}", "-si", "f:1",
                "-cs", design.seed + replicate, "-o", generated, "-ot", 0, "-v", 0,
                "-sb", f"f:{_decimal(design.birth_rate)}", "-sl", f"f:{design.taxa - 1}",
                "-st", f"f:{design.tree_height}", "-so", "f:1",
            ]
            self.runner.run(
                command,
                log=destination.parent / "logs" / "simphy.log",
                timeout=7200,
            )
            produced = generated / "1"
            if not produced.is_dir():
                raise PipelineError(f"SimPhy produced no replicate directory in {generated}")
            shutil.move(str(generated), str(destination))

        outputs = self._outputs(destination)
        if not outputs:
            raise PipelineError(f"SimPhy output is incomplete: {destination}")
        return outputs

    @staticmethod
    def _outputs(destination: Path) -> tuple[Path, list[Path]] | None:
        root = destination / "1"
        species = root / "s_tree.trees"
        genes = sorted(root.glob("g_trees*.trees"))
        return (species, genes) if species.is_file() and genes else None
