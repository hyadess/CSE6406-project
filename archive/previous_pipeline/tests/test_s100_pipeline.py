"""Integration test for the PART A code path.

Part A reads the real S100 distribution, which is too large to bundle. This
test builds a directory with the SAME layout and the SAME label conventions
using the real tools -- SimPhy for the species tree and true gene trees,
AliSim + IQ-TREE for the estimated gene trees -- and then runs run_s100.py
against it.

That exercises everything Part A does on real S100 files: the g_trees*.trees
discovery, the SimPhy <species>_<locus>_<individual> label reconciliation, the
IQ-TREE leading-slash support labels, and the calibration accumulator.

It does NOT prove the real S100 archive contains true gene trees. Nothing
short of downloading it can. `run_s100.py --check` is what tells you that.

    python3 tests/test_s100_pipeline.py
"""
import shutil
import subprocess
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parent.parent / "src"
sys.path.insert(0, str(SRC))
import config
from tools import simphy_run, alisim, iqtree_abayes
from trees import load_tree, load_trees, informative_branches

N_TAXA, N_LOCI, LENGTH = 15, 6, 300
ROOT = Path(config.ROOT) / "work" / "s100_selftest"
REP = "01"


def build():
    if ROOT.exists():
        shutil.rmtree(ROOT)
    rep_dir = ROOT / REP
    rep_dir.mkdir(parents=True)

    sp_out = ROOT / "_simphy"
    sp_tree, true_genes = simphy_run(
        config.SIMPHY_BIN, sp_out,
        n_taxa=N_TAXA, n_loci=N_LOCI,
        birth_rate=config.SIM_BIRTH_RATE, tree_height=config.SIM_TREE_HEIGHT,
        pop_size=config.SIM_POP_SIZE, subst_rate=config.SIM_SUBST_RATE,
        seed=99)
    assert len(true_genes) == N_LOCI, (len(true_genes), N_LOCI)

    # lay the SimPhy output out exactly like an S100 replicate
    shutil.copy(sp_tree, rep_dir / config.S100_SPECIES_TREE)
    for p in true_genes:
        shutil.copy(p, rep_dir / p.name)

    est_dir = rep_dir / "bestMLestimatedgenetree"
    est_dir.mkdir(parents=True)
    scratch = ROOT / "_scratch"
    scratch.mkdir()

    lines = []
    for i, g in enumerate(true_genes):
        aln = alisim(config.IQTREE_BIN, scratch, g, "GTR+G4{0.3}",
                     LENGTH, 500 + i, f"a{i}")
        tre = iqtree_abayes(config.IQTREE_BIN, scratch, aln, "GTR+G4",
                            500 + i, f"e{i}")
        lines.append(tre.read_text().strip())

    target = est_dir / Path(config.S100_EST_GENE_TEMPLATE.format(
        length=LENGTH, suffix=".abayes")).name
    target.write_text("\n".join(lines) + "\n")
    shutil.rmtree(scratch, ignore_errors=True)
    shutil.rmtree(sp_out, ignore_errors=True)
    return target


def main():
    print(f"building an S100-layout replicate with SimPhy + IQ-TREE ...")
    target = build()
    print(f"  estimated gene trees -> {target.relative_to(ROOT)}")

    ests = load_trees(target)
    assert len(ests) == N_LOCI, f"expected {N_LOCI} estimated trees, got {len(ests)}"
    sups = [s for t in ests for _sp, s in informative_branches(t)]
    assert sups and all(s is not None for s in sups), \
        "IQ-TREE support labels did not parse (leading-slash handling?)"
    assert all(0.30 <= s <= 1.0 for s in sups), \
        f"aBayes outside [1/3,1]: min={min(sups)} max={max(sups)}"
    print(f"  ok  {len(sups)} support values parsed, "
          f"range {min(sups):.3f}..{max(sups):.3f}")

    truth = load_tree(sorted((ROOT / REP).glob("g_trees*.trees"))[0])
    est_lbl = {lf.taxon.label for lf in ests[0].leaf_node_iter()}
    true_lbl = {lf.taxon.label for lf in truth.leaf_node_iter()}
    assert est_lbl == true_lbl, "labels should match before reconciliation here"
    print(f"  ok  {len(est_lbl)} taxa, SimPhy-style labels "
          f"e.g. {sorted(est_lbl)[0]}")

    for stage in (["--check"], []):
        cmd = [sys.executable, str(SRC / "run_s100.py"), "--s100-dir", str(ROOT),
               "--reps", REP, "--lengths", str(LENGTH), "--support", ".abayes",
               *stage]
        print(f"\n$ {' '.join(cmd[1:])}")
        r = subprocess.run(cmd, capture_output=True, text=True)
        print(r.stdout[-2200:])
        if r.returncode != 0:
            print(r.stderr[-1500:])
            raise SystemExit(f"run_s100.py failed at stage {stage or 'run'}")

    print("\nPART A code path verified against a real S100-layout directory")


if __name__ == "__main__":
    main()
