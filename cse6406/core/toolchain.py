from __future__ import annotations

import hashlib
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
        fallback = ""
        for command in attempts:
            try:
                result = subprocess.run(
                    command, text=True, stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT, timeout=20, check=False,
                )
                lines = [x.strip() for x in result.stdout.splitlines() if x.strip()]
                explicit = next((x for x in lines if x.lower().startswith("version:")), "")
                if explicit:
                    return explicit.removeprefix("Version:").strip()
                if lines and not fallback:
                    fallback = lines[0]
            except (OSError, subprocess.TimeoutExpired):
                pass
        if name == "simphy":
            return "unreported by executable"
        return fallback or "unreported by executable"

    def metadata(self, name: str) -> dict[str, str]:
        binary = getattr(self, name)
        if binary is None:
            return {"path": "MISSING", "version": "MISSING", "sha256": "MISSING"}
        digest = hashlib.sha256()
        with binary.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return {
            "path": str(binary),
            "version": self.version(name),
            "sha256": digest.hexdigest(),
        }
