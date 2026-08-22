"""Create combined tree files, a locus manifest, and run metadata."""
from __future__ import annotations

import argparse
import csv
import os
from datetime import datetime, timezone
from pathlib import Path

from settings import (CORRECT_MODEL, GENERATING_MODEL, HERE, MISSPECIFIED_MODEL,
                      Experiment, add_experiment_arguments,
                      assert_matching_metadata, experiment_from_args)
from workflow_utils import (alignment_files, combine_newick, inferred_tree_files,
                            find_executable, program_version, require_count,
                            true_gene_trees)


def finalize(experiment: Experiment, simphy: Path, iqtree: Path) -> None:
    truth = true_gene_trees(experiment.out)
    alignments = alignment_files(experiment.out)
    correct = inferred_tree_files(experiment.out, "correct")
    misspecified = inferred_tree_files(experiment.out, "misspecified")

    for files, description in (
        (truth, "true trees"), (alignments, "alignments"),
        (correct, "correct-model trees"),
        (misspecified, "misspecified-model trees"),
    ):
        require_count(files, experiment.loci, description)

    combine_newick(truth, experiment.out / "true_gene_trees.tre")
    combine_newick(correct, experiment.out / "correct_gene_trees.tre")
    combine_newick(misspecified, experiment.out / "misspecified_gene_trees.tre")

    with (experiment.out / "manifest.tsv").open("w", newline="") as handle:
        fields = ["locus", "true_tree", "alignment", "correct_tree",
                  "misspecified_tree"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for index in range(experiment.loci):
            writer.writerow({
                "locus": index + 1,
                "true_tree": truth[index].relative_to(experiment.out),
                "alignment": alignments[index].relative_to(experiment.out),
                "correct_tree": correct[index].relative_to(experiment.out),
                "misspecified_tree": misspecified[index].relative_to(experiment.out),
            })

    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "seed": experiment.seed,
        "taxa": experiment.taxa,
        "loci": experiment.loci,
        "sites_per_locus": experiment.length,
        "generating_model": GENERATING_MODEL,
        "correct_inference_model": CORRECT_MODEL,
        "misspecified_inference_model": MISSPECIFIED_MODEL,
        "simphy_binary": simphy,
        "iqtree_binary": iqtree,
        "iqtree_version": program_version(iqtree),
    }
    (experiment.out / "metadata.tsv").write_text(
        "".join(f"{key}\t{value}\n" for key, value in metadata.items())
    )
    print(f"Finalized combined trees, manifest, and metadata in {experiment.out}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_experiment_arguments(parser)
    parser.add_argument("--simphy", default=os.environ.get("SIMPHY_BIN"))
    parser.add_argument("--iqtree", default=os.environ.get("IQTREE_BIN"))
    args = parser.parse_args()
    experiment = experiment_from_args(args)
    assert_matching_metadata(experiment)
    simphy = find_executable(
        args.simphy, ("simphy",), HERE.parent.parent / "SimPhy/bin/simphy"
    )
    iqtree = find_executable(
        args.iqtree, ("iqtree2", "iqtree"),
        Path.home() / "anaconda3/envs/cse6406/bin/iqtree2",
    )
    finalize(experiment, simphy, iqtree)


if __name__ == "__main__":
    main()
