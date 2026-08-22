#!/usr/bin/env python3
"""Stage 1: generate the species tree and true gene trees with SimPhy."""
from __future__ import annotations

import argparse
import os
import re
import shutil
import tempfile
from pathlib import Path

from settings import (HERE, Experiment, add_experiment_arguments,
                      assert_matching_metadata, experiment_from_args)
from workflow_utils import (PipelineError, find_executable, require_count,
                            run_command, true_gene_trees)


def decimal(value: float) -> str:
    """Write a decimal because SimPhy may reject scientific notation."""
    return f"{value:.12f}".rstrip("0").rstrip(".")


def count_taxa(tree_file: Path) -> int:
    text = tree_file.read_text().strip()
    return len(re.findall(r"(?<=[(,])\s*[^():;,\s]+\s*:", text))


def generate_true_trees(experiment: Experiment, simphy: Path) -> list[Path]:
    trees = true_gene_trees(experiment.out)
    if not trees:
        destination = experiment.out / "simphy"
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists() and any(destination.iterdir()):
            raise PipelineError(f"Partial SimPhy directory exists: {destination}")
        if destination.exists():
            destination.rmdir()

        # SimPhy silently fails when its output path contains a space.
        with tempfile.TemporaryDirectory(prefix="simphy_run_", dir="/tmp") as temp:
            temporary_output = Path(temp) / "output"
            command = [
                simphy,
                "-rs", 1,                         # species-tree replicates
                "-rl", f"f:{experiment.loci}",   # number of loci
                "-rg", 1,                         # gene trees per locus
                "-sp", "f:400000",               # population size
                "-su", "e:10000000", "-si", "f:1",
                "-cs", experiment.seed,
                "-o", temporary_output, "-ot", 0, "-v", 0,
                "-sb", f"f:{decimal(1e-7)}",      # birth rate
                "-sl", f"f:{experiment.taxa - 1}",
                "-st", "f:2500000",              # tree height
                "-so", "f:1",
            ]
            run_command(command, working_dir=None,
                        log=experiment.out / "logs" / "simphy.log")
            generated = sorted((temporary_output / "1").glob("g_trees*.trees"))
            require_count(generated, experiment.loci, "SimPhy gene trees")
            shutil.move(str(temporary_output), str(destination))
        trees = true_gene_trees(experiment.out)

    require_count(trees, experiment.loci, "true gene trees")
    observed_taxa = {count_taxa(path) for path in trees}
    if observed_taxa != {experiment.taxa}:
        raise PipelineError(
            f"Expected {experiment.taxa} taxa in every gene tree; found {observed_taxa}"
        )
    print(f"SimPhy: {len(trees)} true gene trees with {experiment.taxa} taxa")
    return trees


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_experiment_arguments(parser)
    parser.add_argument("--simphy", default=os.environ.get("SIMPHY_BIN"))
    args = parser.parse_args()
    experiment = experiment_from_args(args)
    assert_matching_metadata(experiment)
    simphy = find_executable(
        args.simphy, ("simphy",), HERE.parent.parent / "SimPhy/bin/simphy"
    )
    generate_true_trees(experiment, simphy)


if __name__ == "__main__":
    main()
