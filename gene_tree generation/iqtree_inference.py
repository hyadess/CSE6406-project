"""Shared IQ-TREE logic for the correct and misspecified model stages."""
from __future__ import annotations

import argparse
import os
from pathlib import Path

from settings import (Experiment, add_experiment_arguments,
                      assert_matching_metadata, experiment_from_args)
from workflow_utils import (PipelineError, alignment_files, find_executable,
                            inferred_tree_files, require_count, run_command)


def infer_gene_trees(
    experiment: Experiment, iqtree: Path, *, condition: str, model: str,
) -> list[Path]:
    """Infer one maximum-likelihood tree with aBayes support per alignment."""
    alignments = alignment_files(experiment.out)
    require_count(alignments, experiment.loci,
                  "alignments; run simulate_sequences.py first")
    output_dir = experiment.out / "inferred" / condition
    output_dir.mkdir(parents=True, exist_ok=True)
    generated = 0

    for locus, alignment in enumerate(alignments, 1):
        name = f"locus_{locus:04d}"
        prefix = output_dir / name
        tree_file = prefix.with_suffix(".treefile")
        if not tree_file.exists():
            command = [
                iqtree,
                "-s", alignment.resolve(),
                "-m", model,
                "--abayes",
                "-T", experiment.threads,
                "--seed", experiment.seed + locus,
                "--prefix", prefix,
                "--redo", "--quiet",
            ]
            run_command(command, working_dir=experiment.out,
                        log=experiment.out / "logs" / f"{name}_{condition}.log")
            generated += 1
            print(
                f"IQ-TREE {condition} [{locus:04d}/{experiment.loci}]: "
                f"generated {tree_file.name}"
            )
        if not tree_file.read_text().strip().endswith(";"):
            raise PipelineError(f"Invalid IQ-TREE output: {tree_file}")

    trees = inferred_tree_files(experiment.out, condition)
    require_count(trees, experiment.loci, f"{condition} inferred trees")
    reused = len(trees) - generated
    print(
        f"IQ-TREE {condition}: {len(trees)} trees ready "
        f"({generated} generated, {reused} reused)"
    )
    return trees


def model_stage_main(*, condition: str, model: str) -> None:
    """Command-line entry point shared by the two short model scripts."""
    parser = argparse.ArgumentParser(
        description=f"Infer gene trees under the {condition} model ({model})."
    )
    add_experiment_arguments(parser)
    parser.add_argument("--iqtree", default=os.environ.get("IQTREE_BIN"))
    args = parser.parse_args()
    experiment = experiment_from_args(args)
    assert_matching_metadata(experiment)
    iqtree = find_executable(
        args.iqtree, ("iqtree2", "iqtree"),
        Path.home() / "anaconda3/envs/cse6406/bin/iqtree2",
    )
    infer_gene_trees(experiment, iqtree, condition=condition, model=model)
