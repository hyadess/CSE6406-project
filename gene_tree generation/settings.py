"""Experiment settings shared by every stage."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

HERE = Path(__file__).resolve().parent

# AliSim generates DNA with this heterogeneous GTR model.
GENERATING_MODEL = "GTR{1,2,1,1,2}+F{0.30,0.20,0.20,0.30}+G4{0.20}"

# IQ-TREE analyses the same DNA under these two alternatives.
CORRECT_MODEL = "GTR+F+G4"
MISSPECIFIED_MODEL = "JC"


@dataclass(frozen=True)
class Experiment:
    """Everything that defines one reproducible simulation run."""

    out: Path
    taxa: int
    loci: int
    length: int
    seed: int
    threads: str


def add_experiment_arguments(parser: argparse.ArgumentParser) -> None:
    """Give every stage the same command-line options."""
    parser.add_argument("--out", type=Path, default=HERE / "data" / "practical_50_loci")
    parser.add_argument("--taxa", type=int, default=51)
    parser.add_argument("--loci", type=int, default=50)
    parser.add_argument("--length", type=int, default=500, help="sites per locus")
    parser.add_argument("--seed", type=int, default=6406001)
    parser.add_argument("--threads", default="1", help="IQ-TREE threads per locus")


def experiment_from_args(args: argparse.Namespace) -> Experiment:
    if args.taxa < 4 or args.loci < 1 or args.length < 50:
        raise ValueError("Require taxa >= 4, loci >= 1, and length >= 50")
    return Experiment(
        out=args.out.expanduser().resolve(), taxa=args.taxa, loci=args.loci,
        length=args.length, seed=args.seed, threads=args.threads,
    )


def assert_matching_metadata(experiment: Experiment) -> None:
    """Stop users from mixing two configurations in one output folder."""
    path = experiment.out / "metadata.tsv"
    if not path.exists():
        return
    existing = dict(
        line.split("\t", 1) for line in path.read_text().splitlines() if "\t" in line
    )
    expected = {
        "seed": str(experiment.seed), "taxa": str(experiment.taxa),
        "loci": str(experiment.loci), "sites_per_locus": str(experiment.length),
        "generating_model": GENERATING_MODEL,
        "correct_inference_model": CORRECT_MODEL,
        "misspecified_inference_model": MISSPECIFIED_MODEL,
    }
    conflicts = [
        f"{key}: existing={existing.get(key)!r}, requested={value!r}"
        for key, value in expected.items() if existing.get(key) != value
    ]
    if conflicts:
        raise ValueError(
            "Output directory belongs to another configuration:\n  "
            + "\n  ".join(conflicts) + "\nUse a new --out directory."
        )
