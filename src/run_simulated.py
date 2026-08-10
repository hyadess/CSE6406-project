"""PART B -- regenerate S100-protocol data with the published toolchain.

  SimPhy   -> species tree + TRUE gene trees under the MSC
  AliSim   -> sequences along each true gene tree
  IQ-TREE  -> estimated gene tree with aBayes support
  calibrate-> P(branch in true gene tree | support)

Use this when you need a condition the published S100 dataset does not contain.
The one it is built for is MODEL MISSPECIFICATION: S100 was simulated and
analysed under broadly matching models, so its gene-tree error is largely
stochastic. Under a correctly specified model ML is consistent and error
shrinks with sequence length; only misspecification makes error persist. Both
arms use the SAME generated sequences and differ only in the analysis model.

    python3 src/run_simulated.py --check
    python3 src/run_simulated.py --pilot
    python3 src/run_simulated.py --loci 100 --reps 2
"""
from __future__ import annotations

import argparse
import shutil
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
from tools import simphy_run, alisim, iqtree_abayes, check_tools
from trees import load_tree
from calibrate import Accumulator, score_locus, print_tables


def run(args):
    acc = Accumulator(support_max=1.0)          # aBayes
    work = Path(config.WORK_DIR) / "simulated"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    t0 = time.time()
    for rep in range(1, args.reps + 1):
        seed = config.MASTER_SEED + rep
        sp_dir = work / f"rep{rep}" / "simphy"
        print(f"[rep {rep}] SimPhy: {args.taxa} taxa, {args.loci} loci ...",
              flush=True)
        sp_tree, true_genes = simphy_run(
            config.SIMPHY_BIN, sp_dir,
            n_taxa=args.taxa, n_loci=args.loci,
            birth_rate=config.SIM_BIRTH_RATE,
            tree_height=config.SIM_TREE_HEIGHT,
            pop_size=config.SIM_POP_SIZE,
            subst_rate=config.SIM_SUBST_RATE,
            seed=seed)
        print(f"          {len(true_genes)} true gene trees", flush=True)

        for length in args.lengths:
            for sim_m, ana_m, label in config.SIM_MODEL_PAIRS:
                cond = f"{label}_L{length}"
                wd = work / f"rep{rep}" / f"{label}_L{length}"
                wd.mkdir(parents=True, exist_ok=True)
                t1 = time.time()
                n_ok = 0
                for i, gpath in enumerate(true_genes):
                    s = seed * 1000 + i
                    try:
                        aln = alisim(config.IQTREE_BIN, wd, gpath, sim_m,
                                     length, s, f"a{i}")
                        tre = iqtree_abayes(config.IQTREE_BIN, wd, aln, ana_m,
                                            s, f"e{i}")
                        est = load_tree(tre)
                        true = load_tree(gpath)
                    except Exception as exc:                    # noqa: BLE001
                        print(f"    locus {i} skipped: "
                              f"{type(exc).__name__}: {exc}", flush=True)
                        continue
                    for split, sup, correct in score_locus(est, true):
                        acc.add(cond, f"r{rep}g{i}", split, sup, correct)
                    n_ok += 1
                dt = time.time() - t1
                print(f"  {cond:<22} {n_ok} loci  {dt:6.1f}s  "
                      f"({dt/max(1,n_ok):.2f}s/locus)", flush=True)
                if not args.keep:
                    shutil.rmtree(wd, ignore_errors=True)
        if not args.keep:
            shutil.rmtree(sp_dir.parent, ignore_errors=True)

    outdir = Path(config.RESULT_DIR)
    outdir.mkdir(parents=True, exist_ok=True)
    cal, summ = acc.write(outdir, prefix="partB_simulated")
    print_tables(cal, summ, acc.support_max, config.CONFIDENT_FRACTION)
    if acc.n_no_support:
        print(f"\nwarning: {acc.n_no_support} branches had no parseable "
              f"support label and were skipped")
    print(f"\ntotal {time.time()-t0:.1f}s -> {outdir}/partB_simulated_*.csv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--taxa", type=int, default=config.SIM_N_TAXA)
    ap.add_argument("--loci", type=int, default=config.SIM_N_LOCI)
    ap.add_argument("--reps", type=int, default=config.SIM_REPLICATES)
    ap.add_argument("--lengths", nargs="*", type=int,
                    default=config.SIM_SEQ_LENGTHS)
    ap.add_argument("--pilot", action="store_true")
    ap.add_argument("--keep", action="store_true")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()

    ok = check_tools(config.IQTREE_BIN, config.SIMPHY_BIN)
    if a.check:
        sys.exit(0 if ok else 1)
    if not ok:
        sys.exit(1)
    if a.pilot:
        a.taxa, a.loci, a.reps, a.lengths = 21, 10, 1, [200]
    run(a)


if __name__ == "__main__":
    main()
