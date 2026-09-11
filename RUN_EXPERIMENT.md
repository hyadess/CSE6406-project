# Full Production Experiment Runbook

This runbook executes the complete experiment implemented in
`overall_plan.md`. It is not a demo or smoke test.

The production design is:

| Factor | Production setting |
|---|---|
| Taxa | 50 ingroup + 1 outgroup; one individual per species |
| Loci per replicate | 200 |
| Replicates | 50 |
| Species tree | Conditioned on 50 ingroup leaves and 2,500,000 generations; birth rate `1e-7` |
| Substitution rate | Exponential across species-tree replicates; mean `1e-7` substitutions/site/generation |
| ILS treatments | Haploid Ne 100,000 and 1,000,000; empirical nRF gate required |
| Sequence lengths | 200, 800, and 1,600 bp |
| Generating model | `GTR{1,2,1,1,2}+F{0.30,0.20,0.20,0.30}+G4{0.20}` |
| IQ-TREE analysis models | GTR+G4 control; GTR, HKY+G4, HKY, and JC misspecified conditions |
| Species-tree methods | Resolved/contracted ASTRAL-IV; support-only/hybrid wASTRAL |
| Random seed | 6,406,001 |

> **Running this on more than one core:** every command below drives the
> single-process runner `run_project.py`. `run_project_parallel.py` executes
> the identical design across worker processes and produces byte-identical
> results in a fraction of the wall-clock time. See `RUN_PARALLEL.md`; the
> setup, monitoring, and verification steps here apply unchanged.

This produces 300,000 IQ-TREE gene-tree analyses
(`2 ILS x 50 replicates x 3 lengths x 5 models x 200 loci`) and 1,500 paired
condition rows, each containing four species-tree methods. Expect substantial
CPU time and disk usage.

`+F` is useful because it represents unequal DNA composition. Here it fixes
A/C/G/T frequencies at 0.30/0.20/0.20/0.30 instead of forcing 0.25 each.
`+G4{0.20}` places sites into four speed categories with strong rate variation.
The exact generating parameters and their plain-language interpretation are in
`README.md` and `overall_plan.md`.

## Relationship to the published wASTRAL paper

This project tests the paper-motivated question in `overall_plan.md`: whether
gene-tree model misspecification changes support reliability and therefore the
benefit of support-weighted ASTRAL. It is not an exact reconstruction of the
paper's S100 or S200 datasets. Zhang and Mirarab (2022) used, among other
conditions, 50 replicates, S100/S200 taxon sets, FastTree-2, sequence lengths
of 200/400/800/1,600 bp for S100, gene subsets of 50/200/500/1,000, and several
weighting schemes. The production replicate count is aligned with those
reference simulations, but the current pipeline otherwise implements the controlled
design above using SimPhy, AliSim, IQ-TREE, ASTRAL-IV, and support-weighted
wASTRAL. Do not describe its results as an exact replication of S100 or S200.

## 1. Enter the project and activate its environment

Move into the repository so all relative paths point to the correct files.

```bash
cd /Users/sayemshahad/Downloads/CSE6406-project
```

Activate the Conda environment containing Python, DendroPy, and IQ-TREE.

```bash
conda activate cse6406
```

Install the required Python dependency into the active environment.

```bash
python -m pip install -r requirements.txt
```

## 2. Configure the external programs

Tell the pipeline where the locally installed SimPhy executable is located.

```bash
export SIMPHY_BIN=/Users/sayemshahad/Downloads/SimPhy/bin/simphy
```

Tell the pipeline where the IQ-TREE executable in the Conda environment is located.

```bash
export IQTREE_BIN=/Users/sayemshahad/anaconda3/envs/cse6406/bin/iqtree2
```

Tell the pipeline where the unweighted ASTRAL-IV executable is located.

```bash
export ASTRAL_BIN=/Users/sayemshahad/Downloads/ASTER/bin/astral4
```

Tell the pipeline where the weighted ASTRAL executable is located.

```bash
export WASTRAL_BIN=/Users/sayemshahad/Downloads/ASTER/bin/wastral
```

Confirm that all four configured files exist and are executable.

```bash
for bin in "$SIMPHY_BIN" "$IQTREE_BIN" "$ASTRAL_BIN" "$WASTRAL_BIN"; do test -x "$bin" && echo "OK: $bin" || echo "MISSING: $bin"; done
```

## 3. Validate the software before spending compute time

Confirm that Python is coming from the intended Conda environment.

```bash
which python
```

Confirm that DendroPy imports successfully.

```bash
python -c "import dendropy; print('DendroPy', dendropy.__version__)"
```

Run the current unit-test suite; there is no active `tests/test_trees.py`.

```bash
python -m unittest discover -s tests -v
```

Check all four external tools and record their versions using the exact production design.

```bash
python run_project.py --check --output work/full_experiment_r50 --taxa 51 --loci 200 --replicates 50 --lengths 200 800 1600 --models gtr_g4 gtr hky_g4 hky jc --seed 6406001 --threads AUTO
```

Inspect the recorded tool versions before starting the experiment.

```bash
cat work/full_experiment_r50/preflight.json
```

Inspect the frozen design and confirm that it matches the table above.

```bash
cat work/full_experiment_r50/design.json
```

Do not continue if a test fails, a tool is missing, or the recorded design is
wrong. Use a new output directory if you intentionally change any design
parameter; the existing directory is tied to its `design.json`.

## 4. Run the complete experiment

Start the full production grid in the current terminal.

```bash
python run_project.py --output work/full_experiment_r50 --taxa 51 --loci 200 --replicates 50 --lengths 200 800 1600 --models gtr_g4 gtr hky_g4 hky jc --seed 6406001 --threads AUTO
```

For a long unattended run, use this command instead to keep it alive after the terminal closes.

```bash
nohup python run_project.py --output work/full_experiment_r50 --taxa 51 --loci 200 --replicates 50 --lengths 200 800 1600 --models gtr_g4 gtr hky_g4 hky jc --seed 6406001 --threads AUTO > work/full_experiment_r50/run.log 2>&1 &
```

Immediately record the background process identifier after starting with `nohup`.

```bash
echo $! > work/full_experiment_r50/run.pid
```

The foreground and `nohup` commands are alternatives; do not start both.

## 5. Monitor the production run

Check whether the recorded background process is still running.

```bash
ps -p "$(cat work/full_experiment_r50/run.pid)" -o pid,etime,%cpu,%mem,command
```

Follow the top-level output from the unattended run.

```bash
tail -f work/full_experiment_r50/run.log
```

Count completed IQ-TREE gene trees; the final expected count is 300,000.

```bash
find work/full_experiment_r50 -type f -name '*.treefile' | wc -l
```

Count completed AliSim alignments; the final expected count is 60,000.

```bash
find work/full_experiment_r50 -type f -name 'locus_*.phy' | wc -l
```

Count fully resolved unweighted species trees; the final expected count is 1,500.

```bash
find work/full_experiment_r50 -type f -name 'astral_unweighted.tre' | wc -l
```

Count contracted ASTRAL trees; the final expected count is 1,500.

```bash
find work/full_experiment_r50 -type f -name 'astral_contracted_abayes_0.90.tre' | wc -l

Count support-only and hybrid wASTRAL trees; each final expected count is 1,500.

```bash
find work/full_experiment_r50 -type f -name 'wastral_support.tre' | wc -l
find work/full_experiment_r50 -type f -name 'wastral_hybrid.tre' | wc -l
```
```

Show recently written external-tool logs when diagnosing progress or a failure.

```bash
find work/full_experiment_r50 -type f -name '*.log' -print0 | xargs -0 ls -lt | head -20
```

## 6. Resume an interrupted experiment

Re-export the executable paths in every new shell before resuming.

```bash
export SIMPHY_BIN=/Users/sayemshahad/Downloads/SimPhy/bin/simphy IQTREE_BIN=/Users/sayemshahad/anaconda3/envs/cse6406/bin/iqtree2 ASTRAL_BIN=/Users/sayemshahad/Downloads/ASTER/bin/astral4 WASTRAL_BIN=/Users/sayemshahad/Downloads/ASTER/bin/wastral
```

Resume using exactly the same output directory and experimental parameters.

```bash
python run_project.py --output work/full_experiment_r50 --taxa 51 --loci 200 --replicates 50 --lengths 200 800 1600 --models gtr_g4 gtr hky_g4 hky jc --seed 6406001 --threads AUTO
```

The pipeline reuses validated SimPhy trees, alignments, inferred gene trees,
and species trees already present. Do not delete partial scientific outputs
unless a log demonstrates that a particular file is corrupt.

## 7. Verify and inspect the final results

List every final result table with its file size.

```bash
ls -lh work/full_experiment_r50/results
```

Confirm that the expected seven result tables exist.

```bash
for file in stage1_branches.csv stage1_calibration.csv stage1_summary.csv stage1_replicate_summary.csv ils_verification.csv ils_summary.csv stage2_species_tree_error.csv; do test -s "work/full_experiment_r50/results/$file" && echo "OK: $file" || echo "MISSING: $file"; done
```

Inspect the realized low/high ILS separation before interpreting treatment effects.

```bash
column -s, -t < work/full_experiment_r50/results/ils_summary.csv
```

Inspect the Stage 1 condition-level branch-error and high-support-error summaries.

```bash
column -s, -t < work/full_experiment_r50/results/stage1_summary.csv | less -S
```

Inspect the Stage 1 support-calibration bins.

```bash
column -s, -t < work/full_experiment_r50/results/stage1_calibration.csv | less -S
```

Inspect the paired Stage 2 ASTRAL/wASTRAL species-tree errors.

```bash
column -s, -t < work/full_experiment_r50/results/stage2_species_tree_error.csv | less -S
```

In Stage 2, interpret every `delta = weighted error - baseline error` as follows:

- negative: wASTRAL improved the topology;
- zero: weighting made no topological difference;
- positive: wASTRAL performed worse.

Use replicate-level rows for uncertainty estimates and statistical tests. Do
not treat loci or individual branches from the same replicate as independent
biological replicates. Follow `ANALYSIS_PLAN.md`, including its paired-block
analysis and multiplicity rules. The pipeline now stops before sequence
inference unless `ils_summary.csv` meets all preregistered ILS criteria.
