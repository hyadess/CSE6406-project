# Stage 1 — Gene-tree branch-support calibration

## 1. Purpose

Stage 1 tests a key assumption behind support-weighted summary methods such as
wASTRAL:

> Gene-tree branches with less phylogenetic signal should receive lower support,
> including when gene-tree estimation is biased. Consequently, down-weighting
> low-support branches should reduce the influence of erroneous quartets.

The quantity measured is

\[
P(\text{estimated branch is in the true gene tree}\mid\text{branch support}).
\]

This is a calibration and discrimination measurement at the **gene-tree branch
level**. It is not yet a comparison of estimated species trees. Stage 2 asks
whether the support behavior found here propagates into species-tree accuracy.

### Why correctness uses the true gene tree

Each locus has its own evolutionary history. Because of incomplete lineage
sorting (ILS), a true gene tree can legitimately disagree with the species
tree. Such disagreement is biological and is not a gene-tree estimation error.

Therefore, Stage 1 compares

```text
estimated gene tree for locus i  <->  true gene tree for locus i
```

and deliberately does not compare an estimated gene tree directly with the
species tree. The latter would incorrectly count ILS as estimation error.

## 2. Experimental structure

Stage 1 has three complementary parts.

| Part | Dataset | Ground truth available? | Purpose |
|---|---|---:|---|
| A | Published S100 simulation benchmark | Yes, true gene trees | Main calibration experiment |
| B | Newly simulated data | Yes, true gene trees | Compare correctly specified and deliberately misspecified inference |
| C | Biological avian UCE gene trees | No | Measure the support distribution, not correctness |

The parts answer different questions:

- **Part A:** Is branch support informative and calibrated on the real S100
  benchmark?
- **Part B:** Can systematic model misspecification create highly supported
  wrong branches?
- **Part C:** Does a real biological gene-tree collection have enough support
  variation for support weighting to change the analysis?

## 3. Core measurement

### 3.1 Trees as bipartitions

Every informative internal branch in an unrooted tree divides the taxa into two
sets. For example,

```text
{A, B} | {C, D, E}
```

describes the same split as

```text
{C, D, E} | {A, B}.
```

`src/trees.py` canonicalizes the two representations so the same biological
branch compares equal regardless of which side is encountered first. Leaf
edges, the artificial Newick root, and trivial splits are excluded.

A fully resolved unrooted binary tree with `n` taxa has `n - 3` informative
internal branches. The S100 trees have 101 taxa, so a fully resolved tree has
98 informative branches.

### 3.2 Per-locus scoring

For every estimated tree and its corresponding true tree, `score_locus()`:

1. Finds their shared taxon set.
2. Extracts the true tree's informative splits.
3. Projects splits onto the shared taxa when taxa are missing.
4. Drops a projected split if it becomes trivial.
5. Extracts each informative estimated split and its support.
6. Marks the estimated split `correct=True` if it occurs in the true split set.

This produces one observation per supported estimated branch:

```text
condition, locus, support, correct, split_size
```

### 3.3 Support parsing

The parser accepts support formats encountered in the datasets and tools:

| Raw label | Interpretation |
|---|---:|
| `/0.994` | IQ-TREE aBayes support, 0.994 |
| `0.994` | Plain aBayes support, 0.994 |
| `95` | Bootstrap support, 95 |
| `88/0.97` | Last slash-separated numeric field, 0.97 |

An absent or nonnumeric label is recorded as missing support and excluded from
calibration. The number excluded is always reported.

### 3.4 Support bins

Support is divided by its declared maximum and placed in these bins:

```text
[0.00, 0.50)
[0.50, 0.70)
[0.70, 0.90)
[0.90, 0.95)
[0.95, 0.99)
[0.99, 1.00]
```

The last printed label appears as `[0.99,1)`, but the implementation adds a
small numerical margin and includes support exactly equal to 1.0.

The declared maxima are:

- aBayes: `1.0`
- Bootstrap: `100.0`

For each bin, the program calculates

\[
\widehat P(\text{correct}) =
\frac{\text{number of correct branches}}{\text{number of scored branches}}.
\]

### 3.5 Confidence intervals

The calibration CSV reports Wilson 95% binomial intervals. Wilson intervals
behave sensibly when all branches in a bin are correct or all are wrong, unlike
the simple normal approximation.

These intervals treat branches as independent and are useful as descriptive
intervals. They are not the final inferential uncertainty because branches are
clustered within loci and replicates. Final uncertainty should resample whole
loci or use replicate-level estimates.

## 4. Repository components

| File | Responsibility |
|---|---|
| `src/config.py` | Dataset paths, support scales, simulation parameters, seeds and grids |
| `src/trees.py` | Newick loading, support parsing, split extraction and taxon-label normalization |
| `src/calibrate.py` | Branch scoring, support bins, Wilson intervals and CSV writing |
| `src/run_s100.py` | Part A: published S100 experiment |
| `src/run_simulated.py` | Part B: SimPhy/AliSim/IQ-TREE simulation experiment |
| `src/survey_support.py` | Part C: biological support-distribution survey |
| `src/tools.py` | Wrappers for SimPhy, AliSim and IQ-TREE |
| `tests/test_trees.py` | Unit tests for parsing and split behavior |
| `tests/test_s100_pipeline.py` | Integration test using an S100-like directory |

## 5. Installation

Create and activate a project-local Python environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Confirm DendroPy and run the parsing tests:

```bash
python -c "import dendropy; print(dendropy.__version__)"
python tests/test_trees.py
```

Part A and Part C do not run external phylogenetic programs after their input
trees have been downloaded. Part B additionally requires SimPhy and IQ-TREE 2.
See `INSTALL.md` for their build and configuration instructions.

## 6. Part A — S100 calibration

### 6.1 Dataset

S100 is a published simulation benchmark containing 101 taxa, 1,000 loci per
replicate, and estimated gene trees for sequence lengths of 200, 400, 800 and
1,600 base pairs. The downloaded material used here combines:

- Dryad `S100.zip`: true species trees and estimated gene trees with bootstrap
  and aBayes support.
- GitLab `S101.tar.gz`: the missing true per-locus gene trees in nested
  `genetrees.tar.gz`.

The working dataset root is currently:

```text
data/downloads/S100_dryad
```

Each replicate has this effective structure:

```text
01/
├── truegenetrees
├── s_tree.trees
└── bestMLestimatedgenetree/
    ├── estimatedgenetre_200.gtr.rerooted.final.contracted.non
    ├── estimatedgenetre_200.gtr.rerooted.final.contracted.non.abayes
    ├── ...400...
    ├── ...800...
    └── ...1600...
```

`truegenetrees` contains 1,000 Newick trees, one per locus. The estimated files
also contain one tree per locus. Trees are paired by their order in these files.

### 6.2 Download and extraction commands

Download the Dryad archive:

```bash
mkdir -p data/downloads

curl -L --fail --retry 3 -C - \
  -o data/downloads/S100.zip \
  https://datadryad.org/downloads/file_stream/2374124

unzip -t data/downloads/S100.zip
unzip -q data/downloads/S100.zip -d data/downloads/S100_dryad
```

Download the original GitLab bundle:

```bash
curl -L --fail --retry 3 -C - \
  -o data/downloads/S101_gitlab.tar.gz \
  'https://gitlab.com/esayyari/ASTRALIII/-/raw/master/S101.tar.gz?inline=false'
```

The GitLab archive contains three nested archives:

```text
genetrees.tar.gz
speciestrees.tar.gz
alignments.tar.gz
```

Only the true gene-tree files are needed. Extract them directly into the Dryad
layout without unpacking alignments or duplicate results:

```bash
tar -xOzf data/downloads/S101_gitlab.tar.gz genetrees.tar.gz \
  | tar -xzf - -C data/downloads/S100_dryad {01..50}/truegenetrees
```

Verify that all true-gene files were recovered:

```bash
find data/downloads/S100_dryad -type f -name truegenetrees | wc -l
wc -l data/downloads/S100_dryad/01/truegenetrees
```

Expected values are 50 replicate files and 1,000 trees in replicate 01.

The alignment archive is not needed for either the current Stage 1 calibration
or Stage 2's species-tree analysis.

### 6.3 Preflight command

```bash
python src/run_s100.py \
  --check \
  --s100-dir data/downloads/S100_dryad \
  --reps 01 \
  --lengths 200 \
  --support .abayes
```

Argument meanings:

| Argument | Meaning |
|---|---|
| `--check` | Inspect inputs and exit without performing calibration |
| `--s100-dir` | Override the S100 directory configured in `src/config.py` |
| `--reps 01` | Check replicate 01 only |
| `--lengths 200` | Check the 200-bp estimated gene-tree condition |
| `--support .abayes` | Select the `.abayes` files and support maximum 1.0 |

The completed preflight confirmed:

```text
true gene trees: concatenated, 1 file
estimated taxa: 101
true taxa: 101
label reconciliation: identity
support range: 0.333 to 1.0
informative branches in the first tree: 98/98 labelled
```

The 98 branches match the expected `101 - 3` internal branches.

### 6.4 Commands completed

One replicate, 200 loci, 200 bp:

```bash
python src/run_s100.py \
  --s100-dir data/downloads/S100_dryad \
  --reps 01 \
  --lengths 200 \
  --support .abayes \
  --max-genes 200
```

One-replicate sequence-length sweep:

```bash
python src/run_s100.py \
  --s100-dir data/downloads/S100_dryad \
  --reps 01 \
  --lengths 200 400 800 1600 \
  --support .abayes \
  --max-genes 200
```

Ten-replicate sequence-length sweep, first 200 loci per replicate:

```bash
python src/run_s100.py \
  --s100-dir data/downloads/S100_dryad \
  --reps {01..10} \
  --lengths 200 400 800 1600 \
  --support .abayes \
  --max-genes 200
```

Shell brace expansion turns `{01..10}` into the ten arguments `01`, `02`, ...,
`10`.

### 6.5 Completed full-locus command

The final Part A run used all 1,000 loci in each of 10 replicates:

```bash
python src/run_s100.py \
  --s100-dir data/downloads/S100_dryad \
  --reps {01..10} \
  --lengths 200 400 800 1600 \
  --support .abayes \
  --max-genes 1000
```

This command completed successfully. Its theoretical maximum branch count was:

```text
10 replicates × 4 lengths × 1,000 loci × 98 branches = 3,920,000
```

It scored 3,879,835 supported branches. A further 40,165 informative branches
had no numeric support label and were excluded from calibration.

To analyze bootstrap instead of aBayes support, use an empty suffix:

```bash
python src/run_s100.py \
  --s100-dir data/downloads/S100_dryad \
  --reps {01..10} \
  --lengths 200 400 800 1600 \
  --support "" \
  --max-genes 200
```

Here the declared support maximum becomes 100 rather than 1.

## 7. Final Part A findings

These final Part A results come from 10 replicates, four sequence lengths and
all 1,000 loci per replicate. They are real S100 results, not the earlier
synthetic integration test.

### 7.1 Pooled calibration

| Length | Support bin | Branches | P(correct) |
|---:|---:|---:|---:|
| 200 | <0.50 | 390,607 | 0.147 |
| 200 | 0.50–0.70 | 92,739 | 0.328 |
| 200 | 0.70–0.90 | 109,814 | 0.478 |
| 200 | 0.90–0.95 | 46,007 | 0.605 |
| 200 | 0.95–0.99 | 76,302 | 0.709 |
| 200 | ≥0.99 | 236,821 | 0.941 |
| 400 | <0.50 | 319,305 | 0.185 |
| 400 | 0.50–0.70 | 83,416 | 0.407 |
| 400 | 0.70–0.90 | 96,755 | 0.571 |
| 400 | 0.90–0.95 | 40,820 | 0.702 |
| 400 | 0.95–0.99 | 69,271 | 0.795 |
| 400 | ≥0.99 | 361,493 | 0.961 |
| 800 | <0.50 | 245,237 | 0.224 |
| 800 | 0.50–0.70 | 71,094 | 0.470 |
| 800 | 0.70–0.90 | 80,648 | 0.646 |
| 800 | 0.90–0.95 | 34,167 | 0.773 |
| 800 | 0.95–0.99 | 58,863 | 0.851 |
| 800 | ≥0.99 | 487,266 | 0.972 |
| 1600 | <0.50 | 180,825 | 0.262 |
| 1600 | 0.50–0.70 | 57,559 | 0.514 |
| 1600 | 0.70–0.90 | 63,656 | 0.698 |
| 1600 | 0.90–0.95 | 27,081 | 0.822 |
| 1600 | 0.95–0.99 | 47,956 | 0.884 |
| 1600 | ≥0.99 | 602,133 | 0.978 |

At every sequence length, correctness increases monotonically with support.
Longer sequences also improve correctness within every corresponding support
bin.

### 7.2 Replicate-level results

The current program prints pooled branch-level results, but the branch CSV's
`locus` field retains the replicate identifier. Replicate-level calculations
give:

| Length | Median overall error | Range | Median P(wrong \| support ≥0.95) | Range | Median P(wrong \| support ≥0.99) | Range |
|---:|---:|---:|---:|---:|---:|---:|
| 200 | 0.536 | 0.412–0.639 | 0.109 | 0.081–0.183 | 0.058 | 0.045–0.104 |
| 400 | 0.389 | 0.296–0.546 | 0.059 | 0.047–0.142 | 0.033 | 0.021–0.098 |
| 800 | 0.268 | 0.211–0.460 | 0.033 | 0.021–0.122 | 0.021 | 0.010–0.096 |
| 1600 | 0.184 | 0.146–0.392 | 0.021 | 0.010–0.117 | 0.014 | 0.006–0.099 |

Replicate 01 is a consistent high-error outlier across all four lengths. This
heterogeneity is why pooled branch-level intervals should not be treated as the
final inferential uncertainty.

### 7.3 Pooled headline summaries

| Length | Scored branches | Error rate | Mean support on wrong branches | P(wrong and support ≥0.95) | P(wrong \| support ≥0.95) |
|---:|---:|---:|---:|---:|---:|
| 200 | 952,290 | 0.533 | 0.496 | 0.038 | 0.116 |
| 400 | 971,060 | 0.403 | 0.492 | 0.029 | 0.065 |
| 800 | 977,275 | 0.293 | 0.494 | 0.023 | 0.041 |
| 1600 | 979,210 | 0.208 | 0.500 | 0.019 | 0.029 |

The program's `wrong_and_confident` column is a **joint probability**:

\[
P(\text{wrong and support}\ge 0.95).
\]

It is not the conditional error among confident branches. The conditional
quantity shown in the last column above is often easier to interpret:

\[
P(\text{wrong}\mid\text{support}\ge 0.95).
\]

### 7.4 Missing support

The theoretical maximum is 980,000 branches per length for the completed
10-replicate, 1,000-locus run. Missing-support exclusions were:

| Length | Potential branches | Scored | Missing support |
|---:|---:|---:|---:|
| 200 | 980,000 | 952,290 | 27,710 |
| 400 | 980,000 | 971,060 | 8,940 |
| 800 | 980,000 | 977,275 | 2,725 |
| 1600 | 980,000 | 979,210 | 790 |

The total is 40,165 excluded branches. Direct inspection of the raw Newick
files confirmed that these are genuinely unlabelled internal nodes, not labels
the parser failed to understand. Their frequency decreases sharply with longer
sequences, suggesting that missing support is associated with difficult or
poorly resolved gene trees rather than being random.

### 7.5 Interpretation

The current S100 evidence supports the **directional discrimination assumption**
used by support weighting:

- Wrong branches strongly concentrate at lower support.
- Correctness rises monotonically with support at every sequence length.
- More than 90% of wrong branches are below the 0.95 confidence threshold.
- Confident-branch error decreases as sequence length increases.

However, support is not numerically perfect:

- At 200 bp, the pooled ≥0.99 bin is only 94.1% correct.
- The 0.95–0.99 bin is only 70.9% correct at 200 bp.
- Even at 1,600 bp, the ≥0.99 bin is 97.8% rather than approximately 99%.
- Some replicates, especially replicate 01, retain considerably more
  high-support error than the pooled result suggests.

Thus two statements must be kept separate:

1. **Support is useful for ranking branches by reliability:** strongly
   supported on S100.
2. **The numeric support value is a perfectly calibrated correctness
   probability:** not supported.

Under the original Stage 2 gate, the completed S100 experiment selects the
**reduced D1 Stage 2 experiment**: wrong branches predominantly concentrate at
low support, so the directional support-weighting assumption holds on S100.
This decision was recorded before inspecting any Stage 2 species-tree results.

## 8. Part B — generated correct and misspecified conditions

### 8.1 Pipeline

Part B runs:

```text
SimPhy
  -> species tree and true gene trees under the multispecies coalescent
AliSim
  -> sequence alignment along each true gene tree
IQ-TREE
  -> estimated maximum-likelihood gene tree with aBayes support
calibrate.py
  -> correctness as a function of branch support
```

The generating model is GTR with strong gamma-distributed rate heterogeneity.
Two analysis arms are defined:

| Condition | Analysis model | Purpose |
|---|---|---|
| `correct` | GTR+G4 | Stochastic estimation-error control |
| `misspecified` | JC | Deliberately removes model complexity and creates systematic bias |

Both arms use the same true gene tree, generating model, sequence length and
random seed. The retained pilot alignments were compared byte-for-byte, and all
10 correct/misspecified pairs were identical. The arms therefore differ only
in the IQ-TREE analysis model.

### 8.2 Commands

Check external tools:

```bash
python src/run_simulated.py --check
```

Run the small pilot:

```bash
python src/run_simulated.py --pilot
```

Run a larger grid:

```bash
python src/run_simulated.py \
  --taxa 51 \
  --loci 100 \
  --reps 2 \
  --lengths 200 1600
```

Use `--keep` to retain alignments and inferred gene trees. Without `--keep`,
the script removes intermediate working directories after scoring them.

### 8.3 Pilot findings

The completed pilot used 21 taxa, 10 loci, one replicate and 200-bp sequences.

| Condition | Branches | Error rate | Mean support on wrong branches | P(wrong and confident) |
|---|---:|---:|---:|---:|
| Correct | 190 | 0.516 | 0.573 | 0.021 |
| Misspecified | 190 | 0.674 | 0.841 | 0.274 |

In the ≥0.99 support bin:

- Correct model: 32/33 branches correct.
- Misspecified model: 40/61 branches correct.
- Therefore, 21/61 near-maximum-support branches were wrong under
  misspecification.

This demonstrates that severe model misspecification can produce confidently
wrong branches. It is only a 10-locus pilot in an intentionally severe
condition, so it establishes feasibility of the regime, not its prevalence.

Part B must be scaled to more independent replicates before it can serve as the
formal gate for Stage 2's systematic-error mechanism arm.

## 9. Part C — avian UCE support survey

### 9.1 Purpose and limitation

The bundled avian dataset contains 3,679 biological gene trees with 48 taxa and
RAxML bootstrap support. It has no true gene trees, so it cannot determine
whether a branch is correct. Part C measures only whether branch support has
enough spread for weighting to act on.

### 9.2 Commands

Full bundled dataset:

```bash
python src/survey_support.py
```

First `N` trees only:

```bash
python src/survey_support.py --limit 200
```

Custom Newick file whose support maximum is 1:

```bash
python src/survey_support.py \
  --trees path/to/gene_trees.tre \
  --support-max 1.0 \
  --label my_dataset
```

Argument meanings:

| Argument | Meaning |
|---|---|
| `--trees` | Multi-tree Newick input |
| `--support-max` | Maximum possible support used for normalization |
| `--limit` | Analyze only the first N trees |
| `--label` | Dataset name used in the output filename and CSV rows |

### 9.3 Completed full-data result

The full 3,679-tree survey completed successfully and wrote:

```text
results/partC_support_avian_uce_trees_3679_full.csv
```

Its results are:

- 165,555 informative branches (`3,679 × 45`)
- 0 branches without numeric support
- All 3,679 loci contain all 48 taxa
- Median bootstrap support 22
- 33,204 branches at support 100 (`20.1%`)
- 38,447 branches at or above support 95 (`23.2%`)

The support quantiles are:

| Quantile | Support |
|---:|---:|
| 1% | 0 |
| 5% | 0 |
| 10% | 1 |
| 25% | 5 |
| 50% | 22 |
| 75% | 88 |
| 90% | 100 |
| 95% | 100 |
| 99% | 100 |

This shows substantial support variation, so weighting can materially change
the contribution of branches.

### 9.4 Current CSV caveat

The older `partC_support_avian_uce_trees_3679.csv` contains only 9,000 branches,
exactly `200 × 45`, and came from a `--limit 200` run. It must not be presented
as the full result. The authoritative full result is the newer file whose name
ends in `_full.csv`.

## 10. Output files and column definitions

All outputs are written under `results/`. Re-running the same part overwrites
files with the same prefix.

### 10.1 `*_branches.csv`

One row per scored estimated branch.

| Column | Meaning |
|---|---|
| `condition` | Analysis condition, such as `aBayes_L200` or `misspecified_L200` |
| `locus` | Locus identifier; Part A uses values such as `01_g17`, containing replicate 01 and gene index 17 |
| `support` | Raw branch-support value before normalization |
| `correct` | `True` if the estimated split occurs in the true gene tree; otherwise `False` |
| `split_size` | Number of taxa on the canonical, usually smaller, side of the split |

The exact taxon membership of the split is not saved, only its size. This
limits branch-level auditing and should be improved before a final research
release.

Branches without numeric support do not appear in this CSV. They contribute to
the printed `n_no_support` warning instead.

### 10.2 `*_calibration.csv`

One row per condition and support bin.

| Column | Meaning |
|---|---|
| `condition` | Dataset/model/length condition |
| `support_bin` | Normalized support interval |
| `n_branches` | Number of scored branches in the bin |
| `n_correct` | Number whose split occurs in the true gene tree |
| `p_correct` | `n_correct / n_branches` |
| `ci_lo` | Lower Wilson 95% interval bound |
| `ci_hi` | Upper Wilson 95% interval bound |

### 10.3 `*_summary.csv`

One row per condition.

| Column | Meaning |
|---|---|
| `condition` | Dataset/model/length condition |
| `n_branches` | Total supported branches scored |
| `branch_error_rate` | Fraction with `correct=False` |
| `mean_support_on_wrong` | Mean raw support among wrong branches |
| `wrong_and_confident` | Fraction of all branches that are wrong and have support at least 95% of the maximum |
| `wrong_and_confident_lo` | Lower Wilson bound for the joint fraction |
| `wrong_and_confident_hi` | Upper Wilson bound for the joint fraction |

The summary does not currently include
`P(wrong | confident)`. It can be calculated from the branch CSV as:

```text
number of wrong branches with support >= 0.95
------------------------------------------------
number of all branches with support >= 0.95
```

### 10.4 Part C support CSV

`partC_support_<label>.csv` has:

| Column | Meaning |
|---|---|
| `dataset` | Input label |
| `statistic` | Quantile or summary-statistic name |
| `value` | Calculated value |

Statistics include `p1`, `p5`, `p10`, `p25`, `p50`, `p75`, `p90`, `p95`,
`p99`, `n_branches`, `frac_at_ceiling`, `frac_confident` and `mean`.

## 11. Meaning of the main command-line options

### `run_s100.py`

| Option | Meaning |
|---|---|
| `--reps` | Replicate directory names to analyze |
| `--lengths` | Sequence-length conditions to analyze |
| `--support .abayes` | Read aBayes files and normalize by 1.0 |
| `--support ""` | Read bootstrap files and normalize by 100.0 |
| `--max-genes N` | Use only the first N true and estimated trees per replicate and length |
| `--s100-dir PATH` | Override the configured S100 root |
| `--check` | Validate required files, labels and support parsing without scoring the dataset |

### `run_simulated.py`

| Option | Meaning |
|---|---|
| `--taxa N` | Number of simulated species/taxa |
| `--loci N` | Number of true gene trees and alignments per replicate |
| `--reps N` | Number of independent simulation replicates |
| `--lengths ...` | Alignment lengths to simulate |
| `--pilot` | Override the grid with 21 taxa, 10 loci, one replicate and 200 bp |
| `--keep` | Preserve intermediate simulation, alignment and inferred-tree files |
| `--check` | Check external binaries and DendroPy, then exit |

### `survey_support.py`

| Option | Meaning |
|---|---|
| `--trees PATH` | Input multi-tree Newick file |
| `--support-max X` | Declared maximum support |
| `--limit N` | Use only the first N trees |
| `--label NAME` | Label used in output rows and filename |

## 12. Current conclusions and Stage 2 gate

The evidence currently supports the following statements:

1. On S100, branch support strongly discriminates correct from incorrect
   branches.
2. S100 aBayes support is numerically overconfident, especially for short
   sequences and in the 0.95–0.99 range.
3. High-support errors exist, but most S100 errors remain at low support.
4. Severe model misspecification can generate a different regime in which many
   wrong branches receive near-maximum support, as demonstrated by the small
   Part B pilot.
5. The biological avian data have substantial support variation, but no
   correctness conclusion is possible without true gene trees.

The Stage 2 decision should be separated by mechanism:

- **D1/S100:** the completed 1,000-locus experiment selects the reduced Stage 2
  replication because support is strongly discriminative. This gate is now
  recorded.
- **D2/systematic misspecification:** scale Part B to more independent
  replicates. If high-support error persists, run the full support-weighted
  versus unweighted mechanism comparison.
- **D3/avian:** use only as a real-data disagreement analysis; do not claim one
  species tree is more accurate.

## 13. Important limitations

1. Branch-level Wilson intervals ignore correlation within loci and replicates.
2. Trees are paired by file order rather than an explicit locus identifier.
3. Missing support is more common at short sequence lengths and is unlikely to
   be missing completely at random.
4. Part B's current result is a 10-locus pilot with intentionally severe model
   misspecification.
5. Part B substitutes AliSim and IQ-TREE for S100's INDELible and FastTree2
   toolchain.
6. The older non-`_full` Part C CSV is a limited 200-tree artifact and should
   not be confused with the completed `_full.csv` survey.
7. Fixed output prefixes allow later runs to overwrite earlier results.
8. The branch CSV does not retain exact split membership.
9. Numerical calibration and useful reliability ranking are related but
   distinct properties and must not be conflated.

## 14. Immediate next steps

1. Calculate locus-clustered uncertainty for the completed Part A result and
   save a formal replicate-level summary CSV.
2. Preserve and checksum the completed Part C `_full.csv` artifact.
3. Scale Part B to enough independent replicates to support its own Stage 2
   gate.
4. Preserve the completed Part A artifacts and their run parameters before any
   later command can overwrite them.
5. Proceed to the reduced D1 Stage 2 experiment; decide D2 separately after its
   expanded Stage 1 experiment.
