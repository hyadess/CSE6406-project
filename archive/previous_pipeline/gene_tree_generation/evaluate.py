#!/usr/bin/env python3
"""Score inferred gene-tree bipartitions against each locus's true gene tree."""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import dendropy


def read_trees(path: Path, namespace: dendropy.TaxonNamespace):
    trees = dendropy.TreeList.get(
        path=str(path), schema="newick", taxon_namespace=namespace,
        rooting="force-unrooted", preserve_underscores=True,
        suppress_internal_node_taxa=True,
    )
    for tree in trees:
        tree.is_rooted = False
        tree.deroot()
    return trees


def support(label) -> float | None:
    if label is None:
        return None
    match = re.search(r"(?:^|/)([0-9]*\.?[0-9]+)$", str(label).strip())
    return float(match.group(1)) if match else None


def leafset(node) -> frozenset[str]:
    return frozenset(leaf.taxon.label for leaf in node.leaf_iter())


def canonical(side: frozenset[str], all_taxa: frozenset[str]) -> frozenset[str]:
    other = all_taxa - side
    if len(side) != len(other):
        return side if len(side) < len(other) else other
    return min(side, other, key=lambda x: tuple(sorted(x)))


def branches(tree):
    taxa = frozenset(leaf.taxon.label for leaf in tree.leaf_node_iter())
    for node in tree.preorder_node_iter():
        if node.is_leaf() or node.parent_node is None:
            continue
        side = leafset(node)
        if 1 < len(side) < len(taxa) - 1:
            yield canonical(side, taxa), support(node.label)


def mean(values):
    return sum(values) / len(values) if values else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, type=Path)
    args = ap.parse_args()
    run = args.run.resolve()
    namespace = dendropy.TaxonNamespace()
    truth = read_trees(run / "true_gene_trees.tre", namespace)
    conditions = {
        "correct": read_trees(run / "correct_gene_trees.tre", namespace),
        "misspecified": read_trees(run / "misspecified_gene_trees.tre", namespace),
    }
    if not truth or any(len(x) != len(truth) for x in conditions.values()):
        raise SystemExit("true/correct/misspecified tree counts do not match")

    rows = []
    for condition, estimates in conditions.items():
        for locus, (true_tree, estimated) in enumerate(zip(truth, estimates), 1):
            true_taxa = {x.taxon.label for x in true_tree.leaf_node_iter()}
            est_taxa = {x.taxon.label for x in estimated.leaf_node_iter()}
            if true_taxa != est_taxa:
                raise SystemExit(f"taxon mismatch: {condition}, locus {locus}")
            true_splits = {split for split, _ in branches(true_tree)}
            est = list(branches(estimated))
            matched = sum(split in true_splits for split, _ in est)
            supports = [value for _, value in est if value is not None]
            wrong_supports = [value for split, value in est
                              if split not in true_splits and value is not None]
            rows.append({
                "condition": condition, "locus": locus,
                "taxa": len(true_taxa), "true_branches": len(true_splits),
                "estimated_branches": len(est), "matched_branches": matched,
                "branch_precision": matched / len(est) if est else float("nan"),
                "branch_recall": matched / len(true_splits) if true_splits else float("nan"),
                "normalized_rf": (len(est) + len(true_splits) - 2 * matched) /
                                 (len(est) + len(true_splits)),
                "mean_abayes": mean(supports),
                "mean_wrong_abayes": mean(wrong_supports),
                "wrong_branches_abayes_ge_0.95": sum(x >= 0.95 for x in wrong_supports),
            })

    per_locus = run / "evaluation_per_locus.tsv"
    with per_locus.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    summary_rows = []
    for condition in conditions:
        selected = [r for r in rows if r["condition"] == condition]
        total_branches = sum(r["estimated_branches"] for r in selected)
        wrong_branches = sum(r["estimated_branches"] - r["matched_branches"]
                             for r in selected)
        confident_wrong = sum(r["wrong_branches_abayes_ge_0.95"] for r in selected)
        summary_rows.append({
            "condition": condition,
            "loci": len(selected),
            "estimated_branches": total_branches,
            "wrong_branches": wrong_branches,
            "wrong_branch_fraction": wrong_branches / total_branches,
            "mean_branch_precision": mean([r["branch_precision"] for r in selected]),
            "mean_branch_recall": mean([r["branch_recall"] for r in selected]),
            "mean_normalized_rf": mean([r["normalized_rf"] for r in selected]),
            "mean_abayes": mean([r["mean_abayes"] for r in selected]),
            "mean_wrong_abayes": mean([r["mean_wrong_abayes"] for r in selected]),
            "high_support_wrong_branches": confident_wrong,
            "fraction_wrong_with_abayes_ge_0.95": confident_wrong / wrong_branches,
        })
    summary = run / "evaluation_summary.tsv"
    with summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
    print("\nEvaluation summary")
    for row in summary_rows:
        print(
            f"  {row['condition']}: "
            f"branch accuracy={row['mean_branch_recall']:.3f}, "
            f"normalized RF={row['mean_normalized_rf']:.3f}, "
            f"mean aBayes={row['mean_abayes']:.3f}, "
            f"high-support wrong={row['high_support_wrong_branches']}/"
            f"{row['wrong_branches']} "
            f"({row['fraction_wrong_with_abayes_ge_0.95']:.1%})"
        )
    print(f"\nWrote {per_locus}\nWrote {summary}")


if __name__ == "__main__":
    main()
