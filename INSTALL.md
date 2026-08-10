# INSTALL — exact steps

Every command below was executed on Ubuntu 24.04 (x86-64) while building this
project. Where something failed the first time, the fix is included.

---

## 0. Python

```bash
cd stage1
python3 -m pip install -r requirements.txt      # dendropy >= 5.0
python3 -c "import dendropy; print(dendropy.__version__)"
```
Verified with dendropy 5.0.10, Python 3.11.

---

## 1. IQ-TREE 2 — required for Parts A, B, C

Provides both AliSim (sequence simulation) and ML gene-tree estimation with
aBayes support. **Verified: 2.4.0.**

```bash
cd ..
curl -L -O https://github.com/iqtree/iqtree2/releases/download/v2.4.0/iqtree-2.4.0-Linux-intel.tar.gz
tar xzf iqtree-2.4.0-Linux-intel.tar.gz
./iqtree-2.4.0-Linux-intel/bin/iqtree2 --version
```
On Apple Silicon or ARM use the `-Linux-arm` asset, or install via
`conda install -c bioconda iqtree`.

Then set `IQTREE_BIN` in `src/config.py`.

> **If you use a version other than 2.4.0, run `tests/test_s100_pipeline.py`
> before anything else.** The `--abayes` treefile writes the internal label
> with a leading slash (`(C:0.0,D:0.67)/0.994:0.04`). That format is what the
> parser handles; if it changes, the parser must change with it.

---

## 2. SimPhy — required for Part B only

The simulator used to generate S100 (Mallo, De Oliveira Martins & Posada 2016).

```bash
cd ..
git clone https://github.com/adamallo/SimPhy.git
cd SimPhy
```

The bundled `bin/simphy` is a **macOS Mach-O binary** and will not run on
Linux (`Exec format error`). Delete it and build from source. The build needs
three libraries that are not installed by default:

```bash
sudo apt-get update
sudo apt-get install -y libgsl-dev libmpfr-dev libsqlite3-dev
rm -f bin/simphy
make
./bin/simphy | head -3        # should print the SimPhy banner
```

If `apt-get install` returns `404 Not Found`, run `apt-get update` first —
that was the failure encountered here.

Then set `SIMPHY_BIN` in `src/config.py`.

---

## 3. Datasets

### Part C — avian UCE (bundled, nothing to do)

`data/avian/avian_uce_trees_3679.tre` — 3,679 UCE loci, 48 taxa, RAxML rapid
bootstrap. Avian Phylogenomics Project (Jarvis et al. 2014, *Science*
346:1320-1331), redistributed in the TREE-QMC tutorial (Han & Molloy 2023,
*Genome Res* 33:1042-1052):
<https://github.com/molloy-lab/TREE-QMC/tree/main/tutorial/gene-trees>

### Part A — S100 (you must download it)

The main experiment. Not bundled: it is large and not ours to redistribute.
Zhang, Rabiee, Sayyari & Mirarab (2018), *BMC Bioinformatics* 19:153; also an
evaluation set for wASTRAL (Zhang & Mirarab 2022, *MBE* 39:msac215).

Two sources, and **they do not contain the same things** — check both:

| Source | URL |
|---|---|
| GitLab tarball `S101.tar.gz` | <https://gitlab.com/esayyari/ASTRALIII> |
| Dryad, aBayes + bootstrap gene trees | doi:10.6076/D1WK5R |
| True species trees + published species-tree estimates | <https://github.com/chaoszhang/Weighted-ASTRAL_data> (verified: contains `S100/{01..50}/s_tree.trees`, replicate directories are **zero-padded**) |

Extract so that `data/S100/01/`, `data/S100/02/`, … exist, then:

```bash
python3 src/run_s100.py --check --reps 01 --lengths 200
```

`--check` prints every file it found and every one it did not. It will not
guess. In particular it will tell you whether **true gene trees** are present —
without them the measurement cannot be made at all, because "correct" would be
undefined and gene-tree discordance caused by incomplete lineage sorting would
be wrongly scored as estimation error.

**Two things to confirm by eye before trusting any output**, both printed by
`--check`:

1. **Support scale.** `--check` prints the observed support range and the
   declared maximum. aBayes should fall in [1/3, 1]; bootstrap in [0, 100]. If
   they disagree, fix `S100_SUPPORT` in `src/config.py`.
2. **Label reconciliation.** SimPhy writes gene-tree leaves as
   `<species>_<locus>_<individual>` (`6_0_0`) and species-tree leaves as bare
   integers (`6`). `--check` reports which reconciliation it applied.

---

## 4. Verify the install

```bash
python3 tests/test_trees.py            # no external tools
python3 tests/test_s100_pipeline.py    # needs IQ-TREE and SimPhy
python3 src/run_simulated.py --check
```
