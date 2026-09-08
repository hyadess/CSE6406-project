# Paper-Aligned Subset Experiment

This runbook provides two computationally smaller but scientifically structured
subsets of the experiment in `overall_plan.md`. They are real experiments, not
software demos or smoke tests.

| Option | Taxa | Loci | Replicates | IQ-TREE analyses | Output directory |
|---|---:|---:|---:|---:|---|
| Fast subset | 21 | 10 | 2 | 600 | `work/paper_subset_10x2` |
| Larger subset | 51 | 50 | 5 | 7,500 | `work/paper_subset` |

Use the fast subset first when runtime is the main constraint. The larger
subset instructions remain below for a later confirmation run.

## Subset design

| Factor | Subset setting | Reason retained |
|---|---:|---|
| Taxa | 51 | Keeps the planned project taxon setting |
| Loci per replicate | 50 | Matches one gene-count condition in Zhang and Mirarab (2022) |
| Replicates | 5 | Preserves independent replication while reducing runtime |
| ILS treatments | Low and high | Preserves the biological-discordance contrast |
| Sequence lengths | 200, 800, and 1,600 bp | Preserves the planned error gradient; all occur in the paper's S100 study |
| IQ-TREE models | GTR+G4, GTR, HKY+G4, HKY, and JC | Preserves the complete model-misspecification contrast |
| Species-tree methods | Resolved/contracted ASTRAL-IV and support-only/hybrid wASTRAL | Preserves all four production comparisons |
| Random seed | 6,406,001 | Makes the run reproducible |

The subset performs 7,500 IQ-TREE analyses:

```text
2 ILS levels x 5 replicates x 3 lengths x 5 models x 50 loci = 7,500
```

It also produces 1,500 alignments and 150 paired condition sets, each with all
four species-tree methods. This is 40 times smaller than the 300,000-analysis
production grid.

## Scientific scope

This subset preserves every experimental axis needed to address the project's
question. Only the numbers of loci and replicates are reduced. Analyze the
results as an exploratory subset and report uncertainty across the five
replicates; do not present it as an exact replication of the published S100 or
S200 datasets.

The original S100 study used 100 ingroup taxa plus one outgroup, 1,000
simulated genes, 50 replicates, lengths of 200/400/800/1,600 bp, FastTree-2,
bootstrap support, and multiple gene-count subsets. The present code uses its
own preregistered 51-taxon model-misspecification design, AliSim, IQ-TREE
aBayes support, ASTRAL-IV, and support-weighted wASTRAL. Therefore,
"paper-aligned subset" is the accurate description.

# Fast Subset: 21 Taxa x 10 Loci x 2 Replicates

This is the smallest recommended structured run for this project. It retains
both ILS levels, all three planned sequence lengths, all five gene-tree models,
and all four species-tree methods. It reduces taxa, loci, and independent
replicates to shorten runtime.

| Factor | Fast setting |
|---|---:|
| Taxa | 21 |
| Loci per replicate | 10 |
| Replicates | 2 |
| ILS treatments | Low and high |
| Sequence lengths | 200, 800, and 1,600 bp |
| IQ-TREE models | GTR+G4, GTR, HKY+G4, HKY, and JC |
| Species-tree methods | Resolved/contracted ASTRAL-IV and support-only/hybrid wASTRAL |
| Random seed | 6,406,001 |

The fast subset performs:

```text
2 ILS levels x 2 replicates x 3 lengths x 5 models x 10 loci = 600 IQ-TREE analyses
2 ILS levels x 2 replicates x 3 lengths x 10 loci = 120 alignments
2 ILS levels x 2 replicates x 3 lengths x 5 models = 60 paired species-tree conditions
```

This is paper-faithful in experimental logic, not paper-identical in sample
size. The published S100 study used 101 taxa and did not use a 10-gene,
2-replicate condition. The 21-taxon setting is valid for the implemented
pipeline and matches the repository's documented smaller-production setting,
but it cannot reproduce the taxon-scale behavior of S100 or S200. Treat these
results as a fast exploratory check of the complete experimental structure,
not as publication-strength evidence.

The commands use `--skip-ils-gate` because two replicates and ten loci cannot
reliably certify the production nRF ranges. ILS summaries are still written
and must be inspected; only the automatic stop is disabled.

## A. Prepare and validate the fast subset

Move into the project directory.

```bash
cd /Users/sayemshahad/Downloads/CSE6406-project
```

Activate the project's Conda environment.

```bash
conda activate cse6406
```

Set all required executable paths for this shell.

```bash
export SIMPHY_BIN=/Users/sayemshahad/Downloads/SimPhy/bin/simphy IQTREE_BIN=/Users/sayemshahad/anaconda3/envs/cse6406/bin/iqtree2 ASTRAL_BIN=/Users/sayemshahad/Downloads/ASTER/bin/astral4 WASTRAL_BIN=/Users/sayemshahad/Downloads/ASTER/bin/wastral
```

Run the active unit-test suite before the experiment.

```bash
python -m unittest discover -s tests -v
```

Check the toolchain and freeze the 21-taxon, 10 x 2 design.

```bash
python run_project.py --check --skip-ils-gate --output work/paper_subset_10x2 --taxa 21 --loci 10 --replicates 2 --lengths 200 800 1600 --models gtr_g4 gtr hky_g4 hky jc --seed 6406001 --threads AUTO
```

Review the frozen design before starting the run.

```bash
cat work/paper_subset_10x2/design.json
```

## B. Run the fast subset

Run the experiment in the current terminal.

```bash
python run_project.py --skip-ils-gate --output work/paper_subset_10x2 --taxa 21 --loci 10 --replicates 2 --lengths 200 800 1600 --models gtr_g4 gtr hky_g4 hky jc --seed 6406001 --threads AUTO
```

Alternatively, run it in the background after the preflight command succeeds.

```bash
nohup python run_project.py --skip-ils-gate --output work/paper_subset_10x2 --taxa 21 --loci 10 --replicates 2 --lengths 200 800 1600 --models gtr_g4 gtr hky_g4 hky jc --seed 6406001 --threads AUTO > work/paper_subset_10x2/run.log 2>&1 & echo $! > work/paper_subset_10x2/run.pid
```

Use only one of the two run commands.

## C. Monitor, stop, or resume the fast subset

Check whether the background process is still running.

```bash
ps -p "$(cat work/paper_subset_10x2/run.pid)" -o pid,etime,%cpu,%mem,command
```

Count inferred gene trees; completion requires 600.

```bash
find work/paper_subset_10x2 -type f -name '*.treefile' | wc -l
```

Count simulated alignments; completion requires 120.

```bash
find work/paper_subset_10x2 -type f -name 'locus_*.phy' | wc -l
```

Count each species-tree method; completion requires 60 per command.

```bash
find work/paper_subset_10x2 -type f -name 'astral_unweighted.tre' | wc -l
find work/paper_subset_10x2 -type f -name 'astral_contracted_abayes_0.90.tre' | wc -l
find work/paper_subset_10x2 -type f -name 'wastral_support.tre' | wc -l
find work/paper_subset_10x2 -type f -name 'wastral_hybrid.tre' | wc -l
```

Stop the active child tool and then the background pipeline without deleting outputs.

```bash
pipeline_pid="$(cat work/paper_subset_10x2/run.pid)"; pkill -TERM -P "$pipeline_pid" 2>/dev/null; kill "$pipeline_pid"
```

Resume later with exactly the same parameters.

```bash
python run_project.py --skip-ils-gate --output work/paper_subset_10x2 --taxa 21 --loci 10 --replicates 2 --lengths 200 800 1600 --models gtr_g4 gtr hky_g4 hky jc --seed 6406001 --threads AUTO
```

## D. Verify and inspect the fast subset

List the completed result tables.

```bash
ls -lh work/paper_subset_10x2/results
```

Verify that all seven result tables exist and are nonempty.

```bash
for file in stage1_branches.csv stage1_calibration.csv stage1_summary.csv stage1_replicate_summary.csv ils_verification.csv ils_summary.csv stage2_species_tree_error.csv; do test -s "work/paper_subset_10x2/results/$file" && echo "OK: $file" || echo "MISSING: $file"; done
```

Inspect whether the realized ILS treatments are separated.

```bash
column -s, -t < work/paper_subset_10x2/results/ils_summary.csv
```

Inspect the Stage 1 branch-error summaries.

```bash
column -s, -t < work/paper_subset_10x2/results/stage1_summary.csv | less -S
```

Inspect the paired Stage 2 ASTRAL/wASTRAL errors.

```bash
column -s, -t < work/paper_subset_10x2/results/stage2_species_tree_error.csv | less -S
```

# Larger Subset: 50 Loci x 5 Replicates

## 1. Prepare the shell

Move into the project directory.

```bash
cd /Users/sayemshahad/Downloads/CSE6406-project
```

Activate the project's Conda environment.

```bash
conda activate cse6406
```

Set all four executable paths for this shell.

```bash
export SIMPHY_BIN=/Users/sayemshahad/Downloads/SimPhy/bin/simphy IQTREE_BIN=/Users/sayemshahad/anaconda3/envs/cse6406/bin/iqtree2 ASTRAL_BIN=/Users/sayemshahad/Downloads/ASTER/bin/astral4 WASTRAL_BIN=/Users/sayemshahad/Downloads/ASTER/bin/wastral
```

Confirm that every configured executable is usable.

```bash
for bin in "$SIMPHY_BIN" "$IQTREE_BIN" "$ASTRAL_BIN" "$WASTRAL_BIN"; do test -x "$bin" && echo "OK: $bin" || echo "MISSING: $bin"; done
```

## 2. Validate the subset configuration

Run the current unit-test suite before spending compute time.

```bash
python -m unittest discover -s tests -v
```

Check the complete toolchain and freeze the subset design.

```bash
python run_project.py --check --output work/paper_subset --taxa 51 --loci 50 --replicates 5 --lengths 200 800 1600 --models gtr_g4 gtr hky_g4 hky jc --seed 6406001 --threads AUTO
```

Review the exact frozen experimental design.

```bash
cat work/paper_subset/design.json
```

Review the detected program versions.

```bash
cat work/paper_subset/preflight.json
```

Do not start the experiment if any test or preflight check fails.

## 3. Run the subset experiment

Start the subset in the current terminal.

```bash
python run_project.py --output work/paper_subset --taxa 51 --loci 50 --replicates 5 --lengths 200 800 1600 --models gtr_g4 gtr hky_g4 hky jc --seed 6406001 --threads AUTO
```

Alternatively, start it in the background so it survives closing the terminal.

```bash
nohup python run_project.py --output work/paper_subset --taxa 51 --loci 50 --replicates 5 --lengths 200 800 1600 --models gtr_g4 gtr hky_g4 hky jc --seed 6406001 --threads AUTO > work/paper_subset/run.log 2>&1 & echo $! > work/paper_subset/run.pid
```

Use only one of the two commands above.

## 4. Monitor or stop the run

Check whether the background experiment is still running.

```bash
ps -p "$(cat work/paper_subset/run.pid)" -o pid,etime,%cpu,%mem,command
```

Count completed gene-tree analyses; the final expected count is 7,500.

```bash
find work/paper_subset -type f -name '*.treefile' | wc -l
```

Count completed alignments; the final expected count is 1,500.

```bash
find work/paper_subset -type f -name 'locus_*.phy' | wc -l
```

Count each species-tree method; the final expected count is 150 per command.

```bash
find work/paper_subset -type f -name 'astral_unweighted.tre' | wc -l
find work/paper_subset -type f -name 'astral_contracted_abayes_0.90.tre' | wc -l
find work/paper_subset -type f -name 'wastral_support.tre' | wc -l
find work/paper_subset -type f -name 'wastral_hybrid.tre' | wc -l
```

Stop the active child tool first and then stop the background pipeline, without deleting completed outputs.

```bash
pipeline_pid="$(cat work/paper_subset/run.pid)"; pkill -TERM -P "$pipeline_pid" 2>/dev/null; kill "$pipeline_pid"
```

Check for any child tool that remained after stopping the pipeline.

```bash
pgrep -af "run_project.py|iqtree2|simphy|astral4|wastral"
```

## 5. Resume an interrupted run

Resume with the identical directory and parameters so completed work is reused.

```bash
python run_project.py --output work/paper_subset --taxa 51 --loci 50 --replicates 5 --lengths 200 800 1600 --models gtr_g4 gtr hky_g4 hky jc --seed 6406001 --threads AUTO
```

Do not change the parameters while reusing `work/paper_subset`; its
`design.json` intentionally prevents mixing incompatible experiments.

## 6. Verify the completed experiment

List the final result tables.

```bash
ls -lh work/paper_subset/results
```

Verify that every expected result table is nonempty.

```bash
for file in stage1_branches.csv stage1_calibration.csv stage1_summary.csv stage1_replicate_summary.csv ils_verification.csv ils_summary.csv stage2_species_tree_error.csv; do test -s "work/paper_subset/results/$file" && echo "OK: $file" || echo "MISSING: $file"; done
```

Inspect whether the simulated low- and high-ILS treatments are separated.

```bash
column -s, -t < work/paper_subset/results/ils_summary.csv
```

Inspect gene-tree branch error and high-support incorrect branches by condition.

```bash
column -s, -t < work/paper_subset/results/stage1_summary.csv | less -S
```

Inspect support calibration across models, lengths, and ILS levels.

```bash
column -s, -t < work/paper_subset/results/stage1_calibration.csv | less -S
```

Inspect the four paired ASTRAL and wASTRAL species-tree errors.

```bash
column -s, -t < work/paper_subset/results/stage2_species_tree_error.csv | less -S
```

Interpret each Stage 2 weighted-minus-baseline `delta` column as:

- `delta < 0`: weighting helped;
- `delta = 0`: weighting made no topological difference;
- `delta > 0`: weighting hurt.

Use the replicate-level summaries as the independent units for uncertainty.
Individual branches and loci within the same replicate are not independent
biological replicates. With only five replicates, emphasize effect sizes and
consistency rather than strong claims based only on p-values.
