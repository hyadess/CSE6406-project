#!/usr/bin/env python3
"""Stage 2: simulate one DNA alignment per true tree with IQ-TREE AliSim."""
from __future__ import annotations

import argparse
import os
from pathlib import Path

from settings import (GENERATING_MODEL, Experiment, add_experiment_arguments,
                      assert_matching_metadata, experiment_from_args)
from workflow_utils import (PipelineError, alignment_files, find_executable,
                            require_count, run_command, true_gene_trees)


def simulate_alignments(experiment: Experiment, iqtree: Path) -> list[Path]:
    trees = true_gene_trees(experiment.out)
    require_count(trees, experiment.loci, "true gene trees; run simphy_trees.py first")
    alignment_dir = experiment.out / "alignments"
    alignment_dir.mkdir(parents=True, exist_ok=True)
    generated = 0

    for locus, true_tree in enumerate(trees, 1):
        name = f"locus_{locus:04d}"
        alignment = alignment_dir / f"{name}.phy"
        if not alignment.exists():
            command = [
                iqtree, "--alisim", name,
                "-t", true_tree.resolve(),
                "-m", GENERATING_MODEL,
                "--length", experiment.length,
                "--seed", experiment.seed + locus,
                "--quiet",
            ]
            run_command(command, working_dir=alignment_dir,
                        log=experiment.out / "logs" / f"{name}_alisim.log")
            generated += 1
            print(f"AliSim [{locus:04d}/{experiment.loci}]: generated {alignment.name}")
        if alignment.stat().st_size == 0:
            raise PipelineError(f"Empty alignment: {alignment}")

    alignments = alignment_files(experiment.out)
    require_count(alignments, experiment.loci, "alignments")
    reused = len(alignments) - generated
    print(f"AliSim: {len(alignments)} alignments ready ({generated} generated, {reused} reused)")
    return alignments


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_experiment_arguments(parser)
    parser.add_argument("--iqtree", default=os.environ.get("IQTREE_BIN"))
    args = parser.parse_args()
    experiment = experiment_from_args(args)
    assert_matching_metadata(experiment)
    iqtree = find_executable(
        args.iqtree, ("iqtree2", "iqtree"),
        Path.home() / "anaconda3/envs/cse6406/bin/iqtree2",
    )
    simulate_alignments(experiment, iqtree)


if __name__ == "__main__":
    main()
