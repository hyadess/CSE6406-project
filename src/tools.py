"""Wrappers for the three published tools this project drives.

  SimPhy   Mallo, De Oliveira Martins & Posada (2016), Syst Biol 65:334-344.
           The simulator used to generate the S100 benchmark.
  AliSim   Ly-Trong, Naser-Khdour, Lanfear & Minh (2023), MBE 40:msac092.
           Ships inside IQ-TREE 2.
  IQ-TREE  Minh et al. (2020), MBE 37:1530-1534. Used here for ML gene-tree
           estimation with approximate Bayes support (Anisimova et al. 2011).

Every flag below was checked against the tool's own usage output and then
executed. Verified versions: SimPhy built from the adamallo/SimPhy master
branch; IQ-TREE 2.4.0.

DEVIATION FROM THE PUBLISHED S100 PROTOCOL, stated plainly:
S100 used INDELible for sequence simulation and FastTree2 for gene-tree
estimation. This project substitutes AliSim and IQ-TREE, because IQ-TREE also
supplies the aBayes support that the measurement needs. If you want the S100
protocol exactly, use the real downloaded S100 files (run_s100.py) instead of
regenerating them.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


class ToolError(RuntimeError):
    pass


def _run(cmd, cwd=None, timeout=1800):
    proc = subprocess.run([str(c) for c in cmd], cwd=(str(cwd) if cwd else None),
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          timeout=timeout, text=True)
    if proc.returncode != 0:
        raise ToolError(f"failed (rc={proc.returncode}): "
                        f"{' '.join(str(c) for c in cmd)}\n{proc.stdout[-2500:]}")
    return proc.stdout


# ------------------------------------------------------------------ SimPhy
def simphy_run(simphy_bin, outdir, *, n_taxa, n_loci, birth_rate, tree_height,
               pop_size, subst_rate, seed, species_tree=None, timeout=3600):
    """Simulate a species tree and TRUE gene trees under the MSC.

    Sampling notation (from SimPhy's own usage text): f:x fixed, u:a,b uniform,
    ln:mu,sigma lognormal, e:rate exponential.

    Writes <outdir>/1/s_tree.trees and <outdir>/1/g_trees<N>.trees, which is the
    same layout the real S100 distribution uses.
    """
    outdir = Path(outdir)
    if outdir.exists():
        shutil.rmtree(outdir)
    # SimPhy will not create nested output directories itself; it errors with
    # "The output folder, is not accesible. Errno 2". Make the parent first.
    outdir.parent.mkdir(parents=True, exist_ok=True)

    # SimPhy's parser did not accept scientific notation in testing, so every
    # rate is formatted as a plain decimal.
    def dec(x):
        return f"{float(x):.12f}".rstrip("0").rstrip(".") if float(x) < 1 else f"{float(x):.6f}".rstrip("0").rstrip(".")

    cmd = [simphy_bin, "-rs", 1, "-rl", f"f:{int(n_loci)}", "-rg", 1,
           "-sp", f"f:{int(pop_size)}", "-su", f"e:{dec(subst_rate)}",
           "-si", "f:1", "-cs", int(seed), "-o", str(outdir), "-ot", 0, "-v", 0]
    if species_tree:
        cmd += ["-s", species_tree]
    else:
        cmd += ["-sb", f"f:{dec(birth_rate)}", "-sl", f"f:{int(n_taxa)}",
                "-st", f"f:{int(tree_height)}", "-so", "f:1"]
    _run(cmd, timeout=timeout)

    rep = outdir / "1"
    sp = rep / "s_tree.trees"
    genes = sorted(rep.glob("g_trees*.trees"),
                   key=lambda p: int("".join(c for c in p.stem if c.isdigit()) or 0))
    if not sp.exists() or not genes:
        raise ToolError(f"SimPhy produced no usable output in {rep}")
    return sp, genes


# ------------------------------------------------------------------- AliSim
def alisim(iqtree_bin, workdir, gene_tree_path, model, length, seed, prefix,
           timeout=600):
    """Simulate one alignment along one gene tree. Writes <prefix>.phy."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    _run([iqtree_bin, "--alisim", prefix, "-t", str(Path(gene_tree_path).resolve()),
          "-m", model, "--length", length, "--seed", seed, "--quiet"],
         cwd=workdir, timeout=timeout)
    aln = workdir / f"{prefix}.phy"
    if not aln.exists():
        raise ToolError(f"AliSim wrote no {aln.name}; check model string {model!r}")
    return aln


# ------------------------------------------------------------------ IQ-TREE
def iqtree_abayes(iqtree_bin, workdir, aln_path, model, seed, prefix,
                  timeout=1800):
    """ML gene tree with approximate Bayes branch support.

    `model` is the ANALYSIS model. Deliberately setting it simpler than the
    generating model is how a systematic-bias condition is created; that is a
    choice to report, not to hide.
    """
    workdir = Path(workdir)
    _run([iqtree_bin, "-s", Path(aln_path).name, "-m", model, "--abayes",
          "--prefix", prefix, "--seed", seed, "--quiet", "--redo"],
         cwd=workdir, timeout=timeout)
    tre = workdir / f"{prefix}.treefile"
    if not tre.exists():
        raise ToolError(f"IQ-TREE wrote no {tre.name}")
    return tre


def check_tools(iqtree_bin, simphy_bin=None):
    ok = True
    p = Path(iqtree_bin)
    print(f"[{'ok' if p.exists() else 'MISSING'}] IQ-TREE: {p}")
    if p.exists():
        v = _run([p, "--version"]).splitlines()[0]
        print(f"[ok] {v}")
    else:
        ok = False
    if simphy_bin is not None:
        q = Path(simphy_bin)
        print(f"[{'ok' if q.exists() else 'MISSING'}] SimPhy: {q}")
        ok &= q.exists()
    try:
        import dendropy
        print(f"[ok] dendropy {dendropy.__version__}")
    except ImportError:
        print("[MISSING] dendropy -- pip install -r requirements.txt")
        ok = False
    return ok
