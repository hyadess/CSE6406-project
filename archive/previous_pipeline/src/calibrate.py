"""The measurement.

For every internal branch of an ESTIMATED gene tree, ask whether that branch is
present in the TRUE gene tree for the same locus, and record it against the
branch's support value. Aggregating gives

    P( branch is in the true gene tree | support )

which is the quantity wASTRAL's robustness argument is about. Zhang & Mirarab
(2022) argue that weighting rescues consistency under biased gene-tree
estimation because low-signal gene trees carry low support regardless of bias.
This module measures whether that holds.

Correctness is defined against the TRUE GENE TREE, not the species tree. A
gene-tree branch that conflicts with the species tree because of incomplete
lineage sorting is not an error, and scoring it as one would make the whole
measurement meaningless.

When a gene tree is missing taxa, splits are projected onto the shared taxon
set before comparison (see trees.restrict_split).
"""
from __future__ import annotations

import csv
import math
from collections import defaultdict

from trees import informative_branches, bipartitions, restrict_split

BINS = [(0.0, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 0.95),
        (0.95, 0.99), (0.99, 1.0000001)]


def bin_of(s, smax):
    """Bin a support value after rescaling to [0,1] by its declared maximum."""
    if s is None:
        return None
    x = s / smax
    for lo, hi in BINS:
        if lo <= x < hi:
            return f"[{lo:g},{min(hi,1):g})"
    return None


def score_locus(est_tree, true_tree):
    """Yield (split, support, correct) for each informative branch of est_tree."""
    est_taxa = frozenset(lf.taxon.label for lf in est_tree.leaf_node_iter())
    true_taxa = frozenset(lf.taxon.label for lf in true_tree.leaf_node_iter())
    shared = est_taxa & true_taxa
    if len(shared) < 4:
        return

    true_splits = set()
    for s in bipartitions(true_tree):
        r = restrict_split(s, shared)
        if r is not None:
            true_splits.add(r)

    for split, sup in informative_branches(est_tree):
        r = restrict_split(split, shared)
        if r is None:
            continue
        yield r, sup, (r in true_splits)


def wilson(k, n, z=1.96):
    """Wilson score interval. Correct at k=0 and k=n, unlike the normal
    approximation, which matters because the top support bin is often extreme."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


class Accumulator:
    """Collects branch-level observations and emits the calibration table."""

    def __init__(self, support_max):
        self.support_max = support_max
        self.by_bin = defaultdict(lambda: [0, 0])       # bin -> [correct, total]
        self.rows = []
        self.n_no_support = 0

    def add(self, condition, locus, split, support, correct):
        if support is None:
            self.n_no_support += 1
            return
        b = bin_of(support, self.support_max)
        key = (condition, b)
        self.by_bin[key][0] += int(correct)
        self.by_bin[key][1] += 1
        self.rows.append({
            "condition": condition, "locus": locus,
            "support": support, "correct": correct,
            "split_size": len(split),
        })

    def calibration(self):
        out = []
        for (cond, b), (k, n) in sorted(self.by_bin.items()):
            lo, hi = wilson(k, n)
            out.append({
                "condition": cond, "support_bin": b, "n_branches": n,
                "n_correct": k, "p_correct": round(k / n, 4),
                "ci_lo": round(lo, 4), "ci_hi": round(hi, 4),
            })
        return out

    def summary(self, hi_frac=0.95):
        agg = defaultdict(lambda: [0, 0, 0, 0.0])   # cond -> [n, wrong, wrong_hi, sup_sum_wrong]
        for r in self.rows:
            a = agg[r["condition"]]
            a[0] += 1
            if not r["correct"]:
                a[1] += 1
                a[3] += r["support"]
                if r["support"] / self.support_max >= hi_frac:
                    a[2] += 1
        out = []
        for cond, (n, wrong, wrong_hi, sup_sum) in sorted(agg.items()):
            lo, hip = wilson(wrong_hi, n)
            out.append({
                "condition": cond, "n_branches": n,
                "branch_error_rate": round(wrong / n, 4) if n else None,
                "mean_support_on_wrong": (round(sup_sum / wrong, 4) if wrong else None),
                "wrong_and_confident": round(wrong_hi / n, 4) if n else None,
                "wrong_and_confident_lo": round(lo, 4),
                "wrong_and_confident_hi": round(hip, 4),
            })
        return out

    def write(self, outdir, prefix="calibration"):
        cal, summ = self.calibration(), self.summary()
        for name, rows in (("calibration", cal), ("summary", summ)):
            if not rows:
                continue
            p = outdir / f"{prefix}_{name}.csv"
            with open(p, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
                w.writeheader()
                w.writerows(rows)
        if self.rows:
            p = outdir / f"{prefix}_branches.csv"
            with open(p, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(self.rows[0].keys()))
                w.writeheader()
                w.writerows(self.rows)
        return cal, summ


def print_tables(cal, summ, support_max, hi_frac=0.95):
    print(f"\n{'condition':<28}{'n':>9}{'err':>8}{'sup|wrong':>11}{'wrong&conf':>12}")
    print("-" * 68)
    for s in summ:
        mw = "n/a" if s["mean_support_on_wrong"] is None else f"{s['mean_support_on_wrong']:.3f}"
        print(f"{s['condition']:<28}{s['n_branches']:>9}"
              f"{s['branch_error_rate']:>8.3f}{mw:>11}{s['wrong_and_confident']:>12.3f}")

    print(f"\n{'condition':<28}{'support bin':>14}{'n':>8}{'P(correct)':>12}{'95% CI':>16}")
    print("-" * 78)
    for c in cal:
        print(f"{c['condition']:<28}{c['support_bin']:>14}{c['n_branches']:>8}"
              f"{c['p_correct']:>12.3f}   [{c['ci_lo']:.2f},{c['ci_hi']:.2f}]")

    print(f"\nsupport scale: values divided by {support_max} before binning; "
          f"'confident' means >= {hi_frac:g} of that maximum")
    print("P(correct) = probability the branch is in the TRUE GENE TREE for "
          "that locus.")
    print("wASTRAL's stated assumption implies wrong branches concentrate at "
          "LOW support;")
    print("a large 'wrong&conf' is the direct counter-evidence.")
