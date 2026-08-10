"""PART C -- support-distribution survey on a real gene-tree set.

Runs on the bundled avian UCE data (3,679 loci, 48 taxa, RAxML rapid bootstrap;
Jarvis et al. 2014, redistributed in the TREE-QMC tutorial) or on any Newick
file of gene trees with numeric internal labels.

WHAT THIS CAN AND CANNOT TELL YOU
---------------------------------
It CAN tell you whether support has enough spread for weighting to act on. If
most branches sit at the ceiling, a weighting scheme has almost nothing to
discriminate with, and neither does a capping scheme. That is a real, necessary
preflight before Part A or any weighted analysis.

It CANNOT measure calibration. Biological data has no true gene trees, so
"correct" is undefined, and gene-tree conflict cannot be split into incomplete
lineage sorting versus estimation error. Any claim about accuracy needs Part A
or Part B. Do not let this table stand in for that.

    python3 src/survey_support.py
    python3 src/survey_support.py --trees path/to/genetrees.tre --support-max 1.0
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
from trees import load_trees, informative_branches

QUANTILES = [0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99]


def quantile(sorted_vals, q):
    if not sorted_vals:
        return float("nan")
    i = min(len(sorted_vals) - 1, max(0, int(round(q * (len(sorted_vals) - 1)))))
    return sorted_vals[i]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trees", default=str(config.AVIAN_TREES))
    ap.add_argument("--support-max", type=float,
                    default=config.AVIAN_SUPPORT_MAX)
    ap.add_argument("--limit", type=int, default=None,
                    help="use only the first N gene trees")
    ap.add_argument("--label", default=None)
    a = ap.parse_args()

    path = Path(a.trees)
    if not path.exists():
        raise SystemExit(f"gene tree file not found: {path}")
    label = a.label or path.stem

    print(f"loading {path} ...", flush=True)
    trees = load_trees(path, limit=a.limit)
    print(f"{len(trees)} gene trees", flush=True)

    sups, occupancy, no_label = [], Counter(), 0
    n_branch = 0
    for t in trees:
        occupancy[len(list(t.leaf_node_iter()))] += 1
        for _split, sup in informative_branches(t):
            n_branch += 1
            if sup is None:
                no_label += 1
            else:
                sups.append(sup)

    if not sups:
        raise SystemExit("no numeric support labels found; check the file "
                         "and --support-max")

    sups.sort()
    smax = a.support_max
    ceiling = sum(1 for s in sups if s >= smax - 1e-9)
    conf = sum(1 for s in sups if s / smax >= config.CONFIDENT_FRACTION)

    taxa_counts = sorted(occupancy.items())
    print(f"\ntaxa per gene tree: min={taxa_counts[0][0]} "
          f"max={taxa_counts[-1][0]} "
          f"(complete loci: {occupancy[taxa_counts[-1][0]]})")
    print(f"informative branches: {n_branch}  "
          f"(without a numeric label: {no_label})")
    print(f"\nsupport quantiles (declared maximum {smax:g})")
    rows = []
    for q in QUANTILES:
        v = quantile(sups, q)
        print(f"  p{int(q*100):<3d} {v:>10.4f}")
        rows.append({"dataset": label, "statistic": f"p{int(q*100)}",
                     "value": round(v, 6)})

    print(f"\nat the ceiling (support == {smax:g}): "
          f"{ceiling}/{len(sups)} = {ceiling/len(sups):.3f}")
    print(f"confident (>= {config.CONFIDENT_FRACTION:g} x max): "
          f"{conf}/{len(sups)} = {conf/len(sups):.3f}")
    rows += [
        {"dataset": label, "statistic": "n_branches", "value": len(sups)},
        {"dataset": label, "statistic": "frac_at_ceiling",
         "value": round(ceiling / len(sups), 6)},
        {"dataset": label, "statistic": "frac_confident",
         "value": round(conf / len(sups), 6)},
        {"dataset": label, "statistic": "mean",
         "value": round(sum(sups) / len(sups), 6)},
    ]

    print("\nREADING THIS")
    print("  A large fraction at the ceiling means support carries little")
    print("  discriminating information on this dataset, which limits what any")
    print("  support-weighting or support-capping scheme can do. It says")
    print("  nothing about whether those branches are correct -- that needs")
    print("  true gene trees, i.e. Part A or Part B.")

    outdir = Path(config.RESULT_DIR)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"partC_support_{label}.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["dataset", "statistic", "value"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
