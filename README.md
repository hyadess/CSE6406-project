# Stage 1 — Is gene-tree branch support calibrated when the error is systematic?

## The question

Weighted summary methods (wASTRAL, weighted TREE-QMC, wASTRID) down-weight
low-support gene-tree branches. Zhang & Mirarab (2022) address biased gene-tree
estimation explicitly, including long-branch attraction, and argue that
weighting can rescue consistency where unweighted ASTRAL is positively
misleading. That argument rests on an assumption they state in words:

> gene trees with lower signal have lower support **regardless of bias**, so
> down-weighting low-support branches removes the biased quartets

Stage 1 measures the quantity that assumption is about:

> **P( a gene-tree branch is in the TRUE gene tree | its support )**

**Correctness is defined against the true gene tree, not the species tree.** A
branch that conflicts with the species tree because of incomplete lineage
sorting is not an error. Scoring it as one would make the measurement
meaningless.

This is a **measurement, not a hypothesis test**. If wrong branches carry low
support, the assumption holds and you report that. If they carry high support,
you have found where it fails. Both are results, which is why the project is
staged this way rather than betting on one outcome.

## What is deliberately not claimed

- **No novelty claim.** Someone may have measured something equivalent. The
  framing — *we test a stated assumption* — is defensible either way.
- **No claim about how common any regime is in real data.** Part B's
  misspecified condition is chosen to be severe. Showing a regime exists is not
  showing it is typical.
- **No claim that Part B reproduces S100 exactly.** It substitutes AliSim for
  INDELible and IQ-TREE for FastTree2. Stated in `src/tools.py` and below.
- **No claim that the S100 archive contains true gene trees.** It could not be
  downloaded from the machine this was built on. `run_s100.py --check` is what
  tells you, and it refuses to guess.

---

## The three parts

| | What | Data | Tools | Runs here? |
|---|---|---|---|---|
| **A** | **The main experiment.** Calibration on the real S100 benchmark | S100 (you download) | — | code path verified, data not bundled |
| **B** | Same measurement on freshly generated data, including a systematic-bias condition S100 does not contain | generated | SimPhy, AliSim, IQ-TREE | yes |
| **C** | Support-distribution survey on real biological gene trees | avian UCE, **bundled** | — | yes |

### Datasets

- **S100** — Zhang, Rabiee, Sayyari & Mirarab (2018), *BMC Bioinformatics*
  19:153. Simulated, 101 taxa, 1000 loci, four sequence lengths, ships true
  gene trees and estimated gene trees with aBayes and bootstrap support. Also
  an evaluation set for wASTRAL. **This is the right dataset for the question**
  because it has both truth and support.
- **Avian UCE** — 3,679 loci, 48 taxa, RAxML rapid bootstrap. Avian
  Phylogenomics Project (Jarvis et al. 2014, *Science* 346:1320-1331), taken
  from the TREE-QMC tutorial (Han & Molloy 2023, *Genome Res* 33:1042-1052).
  Bundled in `data/avian/`.

### Tools — all published, none hand-rolled

- **SimPhy** — Mallo, De Oliveira Martins & Posada (2016), *Syst Biol*
  65:334-344. The simulator that generated S100.
- **AliSim** — Ly-Trong et al. (2023), *MBE* 40:msac092. Ships inside IQ-TREE 2.
- **IQ-TREE 2** — Minh et al. (2020), *MBE* 37:1530-1534, with approximate
  Bayes support (Anisimova et al. 2011).

**Deviation from the published S100 protocol, stated plainly:** S100 used
INDELible and FastTree2. Part B substitutes AliSim and IQ-TREE, because
IQ-TREE also supplies the aBayes support the measurement needs. If you want the
S100 protocol exactly, use the real files via Part A.

---

## How to run

Install first — see **INSTALL.md**. It has exact commands and the three
failures encountered while building (macOS binary shipped in the SimPhy repo,
missing GSL/MPFR/SQLite headers, stale apt index).

```bash
cd stage1
python3 -m pip install -r requirements.txt
```

### Verify the install

```bash
python3 tests/test_trees.py            # parsing unit tests, no external tools
python3 tests/test_s100_pipeline.py    # builds an S100-layout replicate with
                                       # SimPhy + IQ-TREE and runs Part A on it
```

### Part C — support survey (runs immediately, bundled data)

```bash
python3 src/survey_support.py
python3 src/survey_support.py --trees mydata.tre --support-max 1.0
```

### Part B — generate data and measure

```bash
python3 src/run_simulated.py --check
python3 src/run_simulated.py --pilot            # ~35 s
python3 src/run_simulated.py --taxa 51 --loci 100 --reps 2 --lengths 200 1600
```

### Part A — the main experiment, on real S100

```bash
# after downloading S100 into data/S100/ per INSTALL.md
python3 src/run_s100.py --check --reps 01 --lengths 200
python3 src/run_s100.py --reps 01 02 03 --lengths 200 400 800 1600
python3 src/run_s100.py --support ""            # bootstrap instead of aBayes
```

`--check` prints every file found and missing and will not proceed on guesses.
Run it before anything else.

### Outputs

Written to `results/`:

- `*_calibration.csv` — P(correct | support bin) with **Wilson 95% intervals**
  (correct at 0 and 1, unlike the normal approximation, which matters because
  the top bin is usually extreme). **This is the headline table.**
- `*_summary.csv` — per condition: branch error rate, mean support on wrong
  branches, and `wrong_and_confident` = fraction of branches that are wrong
  **and** carry support ≥ 95% of maximum. The assumption says bias sits at
  *low* support, so a large value here is the direct counter-evidence.
- `*_branches.csv` — one row per branch, for your own analysis.

---

## Results produced on this machine

### Part C — avian UCE, 3,679 real loci

```
3679 gene trees, 48 taxa each (all loci complete)
165,555 informative branches, 0 without a numeric label

support quantiles (max 100)
  p10    1     p25    5     p50   22     p75   88     p90  100

at the ceiling (== 100):     33,204 / 165,555 = 0.201
confident (>= 95):           38,447 / 165,555 = 0.232
```

Median bootstrap support is **22**. Only about a fifth of branches are
confident. Whatever else is true of this dataset, gene-tree support here has
plenty of spread for a weighting scheme to act on — and a *capping* scheme has
correspondingly little to bite on, since only 20% of branches sit at the
ceiling.

This says nothing about whether those branches are correct. Biological data has
no true gene trees. That is exactly why Part A needs S100.

### Part B — pilot, 21 taxa, 10 loci, 200 bp

Both arms use the **same generated sequences** and differ only in the analysis
model.

| condition | support bin | n | P(correct) | 95% CI |
|---|---|---|---|---|
| correct | [0.99, 1) | 33 | **1.000** | [0.90, 1.00] |
| misspecified | [0.99, 1) | 64 | **0.625** | [0.50, 0.73] |

Under a correctly specified model, maximum likelihood is consistent: branches
at the top of the support scale are right. Under misspecification the support
still saturates at the ceiling — twice as many branches land there — but more
than a third of them are wrong.

That is the shape of the failure the assumption would need to exclude. **It is
a 10-locus pilot in a deliberately severe condition and proves only that the
regime is reachable.** Scale it up and vary the misspecification before drawing
any conclusion.

### Part A code path

Verified against a directory built with SimPhy + IQ-TREE in the real S100
layout — `g_trees*.trees` discovery, SimPhy `6_0_0` label reconciliation,
IQ-TREE leading-slash support labels, and the accumulator. Sanity check passed:
with a correctly specified model, P(correct) in the top support bin was 0.966
[0.83, 0.99].

---

## Things that will silently break this if you skip them

1. **IQ-TREE writes aBayes support with a leading slash** — `/0.994`. A naive
   `float(label)` raises, and the file is invalid input for wASTRAL, whose
   documented requirement is a bare non-negative number. Handled in
   `trees.clean_support`; covered by `tests/test_trees.py`.
2. **SimPhy ships a macOS binary** in `bin/`. On Linux it fails with
   `Exec format error`. Delete and rebuild.
3. **SimPhy will not create nested output directories** and rejects scientific
   notation in its rate flags. Both handled in `tools.simphy_run`.
4. **Label conventions differ between files** — SimPhy gene trees use
   `6_0_0`, species trees use `6`. Reconciled and reported, never assumed.
5. **Missing taxa.** Splits are projected onto the shared taxon set, and any
   split that becomes trivial is dropped rather than scored.
6. **Support scale.** Always confirm the observed range against the declared
   maximum. `--check` prints both.

## Layout

```
INSTALL.md                    exact build steps, with the failures encountered
src/config.py                 paths, grids, documented parameter caveats
src/trees.py                  Newick I/O, splits, label reconciliation
src/calibrate.py              the measurement and Wilson intervals
src/tools.py                  SimPhy / AliSim / IQ-TREE wrappers
src/run_s100.py               PART A — real S100
src/run_simulated.py          PART B — generate and measure
src/survey_support.py         PART C — support survey
tests/test_trees.py           parsing unit tests
tests/test_s100_pipeline.py   Part A integration test on a real-layout dir
data/avian/                   3,679 real avian UCE gene trees
results/                      outputs from the runs above
```

## If the measurement comes out the other way

If wrong branches carry low support across S100, the assumption holds, the
mechanism for expecting weighted methods to fail disappears, and Stage 2 is a
formality. Write that up — it is a clean negative result, and getting one is
the point of measuring rather than assuming.
