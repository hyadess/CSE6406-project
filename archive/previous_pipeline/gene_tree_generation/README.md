# Gene-tree generation: SimPhy to IQ-TREE

This folder contains a complete, reproducible gene-tree experiment—not a toy
tree example. The standard run generates **50 independent loci, 51 taxa per
locus, and 500 nucleotide sites per locus**. Increase `--loci` and `--length`
for the final production grid after the standard run succeeds.

## What is generated

The workflow is:

```text
SimPhy species tree + 50 true gene trees under the multispecies coalescent
                              |
             AliSim: one DNA alignment per true gene tree
                              |
                    the same alignment
                       /              \
       IQ-TREE GTR+F+G4                IQ-TREE JC
       correct model family            misspecified model
```

The generating model is
`GTR{1,2,1,1,2}+F{0.30,0.20,0.20,0.30}+G4{0.20}`. The correct arm estimates
the parameters of the same `GTR+F+G4` family. The misspecified arm uses JC,
which incorrectly assumes equal substitution rates, equal base frequencies,
and no among-site rate variation. Both arms use identical alignments and seeds;
therefore their difference is caused by the inference model, not different data.

SimPhy creates the biological **true gene trees**. "Correct" and
"misspecified" refer to the models used by IQ-TREE to estimate gene trees from
the simulated DNA; they are not two different SimPhy histories.

## Requirements

- Python 3.10+
- [SimPhy](https://github.com/adamallo/SimPhy)
- IQ-TREE 2 with AliSim (`iqtree2`)
- DendroPy for evaluation

From the project root:

```bash
python3 -m pip install -r "gene_tree generation/requirements.txt"
export SIMPHY_BIN=/absolute/path/to/SimPhy/bin/simphy
export IQTREE_BIN=/absolute/path/to/iqtree2
python3 "gene_tree generation/run_pipeline.py" --check
```

The script also searches `PATH`, the project's neighboring `SimPhy` checkout,
and the existing `cse6406` Conda environment.

## How to run it

Run the practical dataset and its evaluation:

```bash
bash "gene_tree generation/run_practical.sh"
```

### Code layout

Each file has one clear responsibility:

| Code file | What it does |
|---|---|
| `settings.py` | experiment defaults and the three evolutionary models |
| `workflow_utils.py` | command logging, file discovery, and validation helpers |
| `simphy_trees.py` | Stage 1: SimPhy species tree and true gene trees |
| `simulate_sequences.py` | Stage 2: AliSim DNA generation on each true tree |
| `iqtree_inference.py` | shared IQ-TREE command used by both inference conditions |
| `correct_model.py` | Stage 3: selects `GTR+F+G4` and creates correct-model trees |
| `misspecified_model.py` | Stage 4: selects `JC` and creates misspecified-model trees |
| `finalize_outputs.py` | combines trees and writes the manifest and metadata |
| `evaluate.py` | compares each inferred topology with its true gene tree |
| `run_pipeline.py` | calls the stages above in order |

The important scientific choices are near the top of `settings.py`. The two
model files are deliberately short, making it immediately visible which model
each condition uses. `iqtree_inference.py` contains only the mechanics common
to both conditions.

### Run one stage at a time

All stages accept the same experiment arguments. From the project root:

```bash
# 1. SimPhy: species tree and true gene trees
python3 "gene_tree generation/simphy_trees.py" \
  --taxa 51 --loci 50 --length 500 --seed 6406001

# 2. AliSim: one shared DNA alignment per true gene tree
python3 "gene_tree generation/simulate_sequences.py" \
  --taxa 51 --loci 50 --length 500 --seed 6406001

# 3. IQ-TREE under the matching GTR+F+G4 model family
python3 "gene_tree generation/correct_model.py" \
  --taxa 51 --loci 50 --length 500 --seed 6406001

# 4. IQ-TREE under the deliberately wrong JC model
python3 "gene_tree generation/misspecified_model.py" \
  --taxa 51 --loci 50 --length 500 --seed 6406001

# 5. Combine per-locus trees and record metadata
python3 "gene_tree generation/finalize_outputs.py" \
  --taxa 51 --loci 50 --length 500 --seed 6406001

# 6. Compare both sets of inferred trees with the SimPhy truth
python3 "gene_tree generation/evaluate.py" \
  --run "gene_tree generation/data/practical_50_loci"
```

Use exactly the same `--out`, `--taxa`, `--loci`, `--length`, and `--seed`
values at every stage. The standard defaults shown above may be omitted.

The pipeline is resumable. Existing valid alignments and tree files are reused,
so an interrupted run can be launched again with the same command. Do not
change `--taxa`, `--loci`, `--length`, models, or seed inside an existing output
directory; use a new `--out` directory for a different experiment.

For a larger final analysis, run replicates in separate output directories:

```bash
python3 "gene_tree generation/run_pipeline.py" \
  --out "gene_tree generation/data/replicate_01_L1000" \
  --taxa 51 --loci 200 --length 1000 --seed 6406101 --threads AUTO

python3 "gene_tree generation/evaluate.py" \
  --run "gene_tree generation/data/replicate_01_L1000"
```

Use multiple independent SimPhy seeds for biological replicates. Sequence
lengths such as 200, 500, 1000, and 1600 bp are useful treatments, but every
condition being compared must reuse its alignment in the two IQ-TREE arms.

## Output files and what they tell you

Inside `data/practical_50_loci/`:

| Output | Meaning |
|---|---|
| `simphy/1/s_tree.trees` | true species tree |
| `simphy/1/g_trees*.trees` | one true gene tree per locus; the scoring truth |
| `alignments/locus_*.phy` | DNA evolved on each corresponding true gene tree |
| `inferred/correct/*.treefile` | IQ-TREE ML trees under `GTR+F+G4`, with aBayes support |
| `inferred/misspecified/*.treefile` | IQ-TREE ML trees under JC, with aBayes support |
| `true_gene_trees.tre` | all true trees as analysis-ready multi-Newick |
| `correct_gene_trees.tre` | all correct-model estimates as multi-Newick |
| `misspecified_gene_trees.tre` | all misspecified estimates as multi-Newick |
| `manifest.tsv` | exact locus-to-alignment-to-tree mapping |
| `metadata.tsv` | seed, dimensions, models, binaries, and IQ-TREE version |
| `evaluation_per_locus.tsv` | topology accuracy and support for each locus |
| `evaluation_summary.tsv` | mean accuracy and high-support errors per condition |
| `logs/` | exact command and full program output for audit/debugging |

`branch_precision` is the fraction of estimated internal splits found in the
corresponding true gene tree. `branch_recall` is the fraction of true splits
recovered. `normalized_rf` is normalized Robinson–Foulds distance: 0 means the
topologies match and larger values mean more disagreement. `mean_abayes`
describes confidence; `high_support_wrong_branches` counts estimated branches
that are absent from the true gene tree despite aBayes support at least 0.95.
`fraction_wrong_with_abayes_ge_0.95` gives that count as a fraction of all
incorrect branches, so runs with different locus counts remain comparable.

Judge model misspecification by comparing conditions across the same loci. A
higher normalized RF distance or more highly supported wrong branches in the
JC arm is direct evidence that the wrong analysis model harms gene-tree
estimation. This comparison concerns estimated versus **true gene trees**, not
estimated gene trees versus the species tree; the latter would incorrectly
count genuine incomplete lineage sorting as estimation error.

## Results of the included practical run

The checked-in `data/practical_50_loci/` run completed all 50 loci (2,400
informative estimated branches per condition) with no skipped loci:

| Metric | Correct `GTR+F+G4` | Misspecified `JC` |
|---|---:|---:|
| True branches recovered | 72.125% | 61.875% |
| Mean normalized RF | 0.27875 | 0.38125 |
| Mean aBayes over all branches | 0.7791 | 0.9507 |
| Mean aBayes on wrong branches | 0.5039 | 0.9008 |
| Wrong branches with aBayes >= 0.95 | 22 / 669 (3.29%) | 534 / 915 (58.36%) |

For this dataset, JC misspecification both reduces topology accuracy and makes
wrong branches look much more certain. The unusually high overall support in
the JC arm is not evidence of better trees: comparison with the known SimPhy
truth shows the opposite. This is a practical 50-locus result and a strong
pipeline check. Treat it as one simulation replicate, not a population-level
claim; final inference should compare multiple independent seeds and sequence
lengths.

## Reproducibility checks performed by the code

The workflow stops rather than silently continuing when the number of SimPhy
gene trees, number of taxa, number of inferred trees, taxon labels, or Newick
termination is wrong. Every locus must have exactly one shared alignment and
one tree in each inference arm. SimPhy 1.0.2 also silently fails when its
output path contains a space; the code works around that upstream limitation
with a space-free temporary directory and moves the validated output here.
