"""PART A -- the measurement on the real S100 benchmark. THE MAIN EXPERIMENT.

S100 is the simulated benchmark of Zhang, Rabiee, Sayyari & Mirarab (2018) and
one of the evaluation sets for wASTRAL (Zhang & Mirarab 2022). It is the right
dataset for this question because it ships BOTH true gene trees and estimated
gene trees with branch support, at four sequence lengths spanning a wide range
of gene-tree estimation error.

The dataset is NOT bundled -- it is large and not ours to redistribute. See
INSTALL.md. Run `--check` first: it reports exactly which files it found and
which it did not, and refuses to guess.

    python3 src/run_s100.py --check
    python3 src/run_s100.py --reps 01 02 --lengths 200 1600
    python3 src/run_s100.py --support ""            # bootstrap instead of aBayes
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
from trees import load_trees, load_tree, strip_simphy_suffix
from calibrate import Accumulator, score_locus, print_tables


# --------------------------------------------------------------- discovery
S100_ROOT = None          # set by --s100-dir; defaults to config.S100_DIR


def s100_root():
    return Path(S100_ROOT) if S100_ROOT else Path(config.S100_DIR)


def replicate_dir(rep):
    return s100_root() / rep


def est_path(rep, length, suffix):
    return replicate_dir(rep) / config.S100_EST_GENE_TEMPLATE.format(
        length=length, suffix=suffix)


def find_true_gene_trees(rep):
    """Locate the true gene trees for a replicate.

    Two layouts are accepted, because the distribution has been packaged more
    than one way:
      * many files matching g_trees*.trees (SimPhy's native output)
      * one multi-tree file named truegenetrees / true.gt.trees / g_trees.trees
    Returns (kind, paths). kind is "per_locus" or "concatenated" or None.
    """
    d = replicate_dir(rep)
    if not d.is_dir():
        return None, []
    per = sorted(d.glob(config.S100_TRUE_GENE_GLOB),
                 key=lambda p: int("".join(c for c in p.stem if c.isdigit()) or 0))
    if len(per) > 1:
        return "per_locus", per
    for name in ("truegenetrees", "true.gt.trees", "g_trees.trees",
                 "truegenetrees.trees"):
        p = d / name
        if p.exists():
            return "concatenated", [p]
    if len(per) == 1:
        return "concatenated", per
    return None, []


def align_labels(true_tree, est_tree):
    """Make the two trees share a taxon vocabulary.

    S100 gene trees may carry SimPhy's <species>_<locus>_<individual> labels or
    bare species labels depending on how they were packaged. Try the identity
    first, then stripping the suffix from either side. Returns a short string
    naming what was done, or None if the label sets cannot be reconciled.
    """
    t = {lf.taxon.label for lf in true_tree.leaf_node_iter()}
    e = {lf.taxon.label for lf in est_tree.leaf_node_iter()}
    if t & e:
        return "identity"
    if strip_simphy_suffix(true_tree):
        t = {lf.taxon.label for lf in true_tree.leaf_node_iter()}
        if t & e:
            return "stripped_true"
    if strip_simphy_suffix(est_tree):
        e = {lf.taxon.label for lf in est_tree.leaf_node_iter()}
        if t & e:
            return "stripped_est"
    return None


# ------------------------------------------------------------------ check
def check(reps, lengths, suffix):
    root = s100_root()
    print(f"S100_DIR = {root}")
    if not root.is_dir():
        print(f"[MISSING] directory does not exist. See INSTALL.md.")
        return False
    ok = True
    for rep in reps:
        d = replicate_dir(rep)
        if not d.is_dir():
            print(f"[MISSING] replicate {rep}: {d}")
            ok = False
            continue
        sp = d / config.S100_SPECIES_TREE
        print(f"[{'ok' if sp.exists() else 'MISSING'}] {rep} species tree: {sp.name}")
        kind, paths = find_true_gene_trees(rep)
        if kind is None:
            print(f"[MISSING] {rep} TRUE gene trees. Looked for "
                  f"{config.S100_TRUE_GENE_GLOB} and truegenetrees in {d}.")
            print(f"          Without true gene trees this measurement cannot "
                  f"be made: 'correct' would be undefined and species-tree "
                  f"discordance from ILS would be miscounted as error.")
            print(f"          The Dryad record (doi:10.6076/D1WK5R) and the "
                  f"GitLab tarball differ in what they include -- check both.")
            ok = False
        else:
            print(f"[ok] {rep} true gene trees: {kind}, {len(paths)} file(s)")
        for L in lengths:
            p = est_path(rep, L, suffix)
            print(f"[{'ok' if p.exists() else 'MISSING'}] {rep} estimated "
                  f"{L}bp: {p.relative_to(root) if p.exists() else p}")
            ok &= p.exists()

    if ok:
        rep, L = reps[0], lengths[0]
        est = load_trees(est_path(rep, L, suffix), limit=1)[0]
        kind, paths = find_true_gene_trees(rep)
        true = (load_tree(paths[0]) if kind == "per_locus"
                else load_trees(paths[0], limit=1)[0])
        mode = align_labels(true, est)
        n_est = len(list(est.leaf_node_iter()))
        n_true = len(list(true.leaf_node_iter()))
        print(f"[{'ok' if mode else 'FAIL'}] label reconciliation: {mode} "
              f"(estimated {n_est} taxa, true {n_true} taxa)")
        sups = [s for _sp, s in __import__("trees").informative_branches(est)]
        good = [s for s in sups if s is not None]
        print(f"[{'ok' if good else 'FAIL'}] support parsed on "
              f"{len(good)}/{len(sups)} informative branches; "
              f"range {min(good) if good else '-'}..{max(good) if good else '-'}")
        print(f"     declared maximum for suffix {suffix!r}: "
              f"{config.S100_SUPPORT[suffix][0]} -- CONFIRM this matches the "
              f"range above before trusting any binning.")
        ok &= bool(mode) and bool(good)
    return ok


# -------------------------------------------------------------------- run
def run(reps, lengths, suffix, max_genes):
    smax, sname = config.S100_SUPPORT[suffix]
    acc = Accumulator(support_max=smax)
    modes = set()

    for rep in reps:
        kind, tpaths = find_true_gene_trees(rep)
        if kind is None:
            print(f"skipping replicate {rep}: no true gene trees")
            continue
        true_trees = ([load_tree(p) for p in tpaths[:max_genes]]
                      if kind == "per_locus"
                      else list(load_trees(tpaths[0], limit=max_genes)))
        for L in lengths:
            ests = load_trees(est_path(rep, L, suffix), limit=max_genes)
            n = min(len(ests), len(true_trees))
            cond = f"{sname}_L{L}"
            for i in range(n):
                est, true = ests[i], true_trees[i]
                mode = align_labels(true, est)
                if mode is None:
                    continue
                modes.add(mode)
                for split, sup, correct in score_locus(est, true):
                    acc.add(cond, f"{rep}_g{i}", split, sup, correct)
            print(f"  {rep} {cond}: {n} loci", flush=True)

    outdir = Path(config.RESULT_DIR)
    outdir.mkdir(parents=True, exist_ok=True)
    cal, summ = acc.write(outdir, prefix=f"partA_S100_{sname}")
    if not cal:
        raise SystemExit("no branches scored -- run --check")
    print_tables(cal, summ, acc.support_max, config.CONFIDENT_FRACTION)
    print(f"\nlabel reconciliation used: {sorted(modes)}")
    if acc.n_no_support:
        print(f"warning: {acc.n_no_support} branches had no parseable support")
    print(f"wrote {outdir}/partA_S100_{sname}_*.csv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", nargs="*", default=config.S100_REPLICATES)
    ap.add_argument("--lengths", nargs="*", type=int, default=config.S100_LENGTHS)
    ap.add_argument("--support", default=".abayes",
                    choices=list(config.S100_SUPPORT.keys()))
    ap.add_argument("--max-genes", type=int, default=config.S100_MAX_GENES)
    ap.add_argument("--s100-dir", default=None,
                    help="override config.S100_DIR (used by the self-test)")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()

    global S100_ROOT
    S100_ROOT = a.s100_dir

    ok = check(a.reps, a.lengths, a.support)
    if a.check:
        sys.exit(0 if ok else 1)
    if not ok:
        sys.exit("preflight failed; fix the items above (see INSTALL.md)")
    run(a.reps, a.lengths, a.support, a.max_genes)


if __name__ == "__main__":
    main()
