#!/usr/bin/env python3
"""Run every gene-tree generation stage in its scientific order."""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from finalize_outputs import finalize
from iqtree_inference import infer_gene_trees
from settings import (CORRECT_MODEL, HERE, MISSPECIFIED_MODEL,
                      add_experiment_arguments, assert_matching_metadata,
                      experiment_from_args)
from simphy_trees import generate_true_trees
from simulate_sequences import simulate_alignments
from workflow_utils import PipelineError, find_executable, program_version


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_experiment_arguments(parser)
    parser.add_argument("--simphy", default=os.environ.get("SIMPHY_BIN"))
    parser.add_argument("--iqtree", default=os.environ.get("IQTREE_BIN"))
    parser.add_argument("--check", action="store_true", help="check tools, then exit")
    args = parser.parse_args()

    experiment = experiment_from_args(args)
    simphy = find_executable(
        args.simphy, ("simphy",), HERE.parent.parent / "SimPhy/bin/simphy"
    )
    iqtree = find_executable(
        args.iqtree, ("iqtree2", "iqtree"),
        Path.home() / "anaconda3/envs/cse6406/bin/iqtree2",
    )
    print(f"SimPhy: {simphy}")
    print(f"IQ-TREE: {iqtree} ({program_version(iqtree)})")
    if args.check:
        return 0

    assert_matching_metadata(experiment)

    # Each function below is also available as its own readable script.
    generate_true_trees(experiment, simphy)
    simulate_alignments(experiment, iqtree)
    infer_gene_trees(experiment, iqtree, condition="correct", model=CORRECT_MODEL)
    infer_gene_trees(
        experiment, iqtree, condition="misspecified", model=MISSPECIFIED_MODEL,
    )
    finalize(experiment, simphy, iqtree)

    print(f"\nCompleted: {experiment.out}")
    print(f"Next: python evaluate.py --run {experiment.out}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (PipelineError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
