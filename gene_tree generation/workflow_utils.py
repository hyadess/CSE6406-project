"""Helpers for finding tools, running commands, and validating files."""
from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
from pathlib import Path


class PipelineError(RuntimeError):
    pass


def find_executable(requested: str | None, names: tuple[str, ...], fallback: Path) -> Path:
    candidates = [Path(requested).expanduser()] if requested else []
    candidates += [Path(found) for name in names if (found := shutil.which(name))]
    candidates.append(fallback)
    for path in candidates:
        if path.is_file() and os.access(path, os.X_OK):
            return path.resolve()
    raise PipelineError(f"Could not find {names[0]}. Checked: {candidates}")


def run_command(command: list[object], *, working_dir: Path | None, log: Path) -> None:
    """Run one program and save its exact command and output in a log."""
    display = " ".join(shlex.quote(str(part)) for part in command)
    result = subprocess.run(
        [str(part) for part in command], cwd=working_dir, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=3600,
    )
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(f"$ {display}\n\n{result.stdout}")
    if result.returncode:
        raise PipelineError(f"Command failed ({result.returncode}); see {log}")


def program_version(binary: Path) -> str:
    result = subprocess.run(
        [str(binary), "--version"], text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=20,
    )
    return next((line.strip() for line in result.stdout.splitlines() if line.strip()), "unknown")


def locus_number(path: Path) -> int:
    match = re.search(r"(\d+)$", path.stem)
    if not match:
        raise PipelineError(f"Cannot read locus number from {path.name}")
    return int(match.group(1))


def true_gene_trees(output: Path) -> list[Path]:
    return sorted((output / "simphy" / "1").glob("g_trees*.trees"), key=locus_number)


def alignment_files(output: Path) -> list[Path]:
    return sorted((output / "alignments").glob("locus_*.phy"), key=locus_number)


def inferred_tree_files(output: Path, condition: str) -> list[Path]:
    return sorted((output / "inferred" / condition).glob("locus_*.treefile"),
                  key=locus_number)


def require_count(files: list[Path], expected: int, description: str) -> None:
    if len(files) != expected:
        raise PipelineError(f"Expected {expected} {description}; found {len(files)}")


def combine_newick(tree_files: list[Path], destination: Path) -> None:
    trees = []
    for path in tree_files:
        tree = path.read_text().strip()
        if not tree.endswith(";"):
            raise PipelineError(f"Incomplete Newick tree: {path}")
        trees.append(tree)
    destination.write_text("\n".join(trees) + "\n")
