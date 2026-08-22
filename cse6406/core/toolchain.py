from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .errors import PipelineError


def _resolve(configured: str | None, names: tuple[str, ...]) -> Path | None:
    if configured:
        path = Path(configured).expanduser()
        return path.resolve() if path.is_file() and os.access(path, os.X_OK) else None
    for name in names:
        found = shutil.which(name)
        if found:
            return Path(found).resolve()
    return None


@dataclass(frozen=True)
class Toolchain:
    """Resolved external programs used by the production workflow."""

    simphy: Path | None
    iqtree: Path | None
    astral: Path | None
    wastral: Path | None

    @classmethod
    def discover(cls) -> "Toolchain":
        return cls(
            simphy=_resolve(os.environ.get("SIMPHY_BIN"), ("simphy",)),
            iqtree=_resolve(os.environ.get("IQTREE_BIN"), ("iqtree3", "iqtree2", "iqtree")),
            astral=_resolve(os.environ.get("ASTRAL_BIN"), ("astral4",)),
            wastral=_resolve(os.environ.get("WASTRAL_BIN"), ("wastral",)),
        )

    def require(self, *names: str) -> None:
        missing = [name for name in names if getattr(self, name) is None]
        if missing:
            variables = ", ".join(f"{name.upper()}_BIN" for name in missing)
            raise PipelineError(
                f"Missing executable(s): {', '.join(missing)}. Install them or set {variables}."
            )

    def version(self, name: str) -> str:
        binary = getattr(self, name)
        if binary is None:
            return "MISSING"
        attempts = ([str(binary), "--version"], [str(binary), "-h"])
        for command in attempts:
            try:
                result = subprocess.run(
                    command, text=True, stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT, timeout=20, check=False,
                )
                line = next((x.strip() for x in result.stdout.splitlines() if x.strip()), "")
                if line:
                    return line
            except (OSError, subprocess.TimeoutExpired):
                pass
        return str(binary)
