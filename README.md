# Model misspecification, branch support, and weighted ASTRAL

This repository implements the controlled experiment in `overall_plan.md`:

```text
SimPhy low/high ILS true histories
        -> AliSim GTR+G4 DNA
        -> IQ-TREE under GTR+G4, GTR, HKY+G4, HKY, and JC
        -> Stage 1 branch correctness and support calibration
        -> resolved and 0.90-aBayes-contracted ASTRAL-IV
        -> support-only and default-hybrid wASTRAL
        -> normalized RF error against the true species tree
```

Correctness in Stage 1 is always measured against the corresponding **true
gene tree**, not the species tree. This keeps biological discordance from ILS
separate from gene-tree estimation error.

## Run

See `INSTALL.md` for complete installation and platform setup.

```bash
conda activate cse6406
python run_project.py --check
python -m unittest discover -s tests -v
python run_project.py --threads AUTO
```

The default grid is a production experiment and requires substantial compute.
It uses 50 replicate blocks, matching the simulation replicate count in the
reference wASTRAL S100/S200 analyses. A smaller pilot command is provided in
`INSTALL.md`; confirmatory inference follows `ANALYSIS_PLAN.md`.

## What GTR+F+G4 means

AliSim generates every alignment with the explicit model

```text
GTR{1,2,1,1,2}+F{0.30,0.20,0.20,0.30}+G4{0.20}
```

In plain language:

- **GTR** allows the six possible DNA changes (A-C, A-G, A-T, C-G, C-T,
  and G-T) to occur at different relative speeds. IQ-TREE fixes the final G-T
  rate to 1 as the reference, so the five supplied values define the other
  rates: `1, 2, 1, 1, 2`.
- **`+F{0.30,0.20,0.20,0.30}`** sets the long-run DNA composition to 30% A,
  20% C, 20% G, and 30% T. The `F` is useful because real DNA need not contain
  the four bases equally. Without an unequal-frequency component, a model such
  as JC assumes 25% of each base and deliberately ignores this signal.
- **`+G4{0.20}`** says that sites evolve at different speeds. AliSim divides a
  gamma distribution into four rate categories; alpha 0.20 creates strong
  variation, with many slowly changing sites and a smaller group of fast sites.

In an IQ-TREE *analysis* model, plain `+F` means frequencies counted from the
alignment. In this AliSim *generating* model, braces provide the actual fixed
frequencies, so `+F{...}` is a complete simulation parameter rather than an
undefined instruction. The same generating parameters are held fixed across
loci to isolate inference-model effects; lack of across-locus parameter
heterogeneity is therefore a deliberate limitation.

## Tree-simulation parameters

The 51 taxa are 50 ingroup species plus one outgroup, with one sampled
individual per species. SimPhy conditions the ingroup tree on 50 leaves and a
2,500,000-generation height, using birth rate `1e-7`. The substitution rate is
drawn once per species-tree replicate from an exponential distribution with
mean `1e-7` substitutions per site per generation. In SimPhy syntax this is
`-su e:10000000`, because the exponential argument is the reciprocal mean.

Low and high ILS use haploid effective population sizes of 100,000 and
1,000,000 on paired species trees generated with the same replicate seed. The
pipeline requires pooled true-gene/species nRF <=0.25 for low ILS, >=0.55 for
high ILS, and a gap >=0.30 before starting sequence inference.

## Code layout

Each major class is in its own readable file:

```text
cse6406/core/          command execution, errors, executable discovery
cse6406/domain/        models, conditions, and the factorial design
cse6406/inference_models/ one file for each GTR/HKY/JC inference option
cse6406/trees/         Newick parsing, support handling, topology distance
cse6406/simulation/    SimPhy and AliSim classes
cse6406/gene_trees/    IQ-TREE inference class
cse6406/stage1/        branch scorer, calibration, and ILS verification
cse6406/stage2/        ASTRAL, wASTRAL, input preparation, species-tree scoring
cse6406/reporting/     stable CSV writer
cse6406/pipeline.py    resumable orchestration
run_project.py         command-line entry point
```

Inference choices are intentionally separated and easy to inspect:

```text
cse6406/inference_models/gtr_g4.py  GTR+G4 reference
cse6406/inference_models/gtr.py     GTR without gamma
cse6406/inference_models/hky_g4.py  HKY+G4
cse6406/inference_models/hky.py     HKY without gamma
cse6406/inference_models/jc.py      severe JC stress test
cse6406/inference_models/registry.py fixed run order and name lookup
```

Run all five models by default, or select individual files by name:

```bash
python run_project.py --models gtr_g4 gtr hky_g4 hky jc
python run_project.py --models gtr_g4 jc
```

Superseded code and historical outputs are isolated under
`archive/previous_pipeline/`; nothing there is imported by the active project.

## Scientific outputs

Stage 1 reports overall branch error, calibration in fixed support bins,
mean support per bin, Brier score, expected calibration error,
`P(wrong | support >= 0.95)`, `P(wrong | support >= 0.99)`, missing support,
and replicate-level summaries. ILS is checked independently as true-gene-tree
versus true-species-tree normalized RF discordance.

Stage 2 runs four reference-aligned analyses on the same inferred gene trees:

- ASTRAL-IV on fully resolved gene trees;
- ASTRAL-IV after contracting branches with aBayes below 0.90;
- support-only wASTRAL using `--mode 2 -B`;
- default hybrid wASTRAL using `-B`.

The primary deltas compare weighted methods with both unweighted baselines:

```text
delta = weighted normalized RF error - baseline normalized RF error
```

Negative values mean weighting helped; positive values mean it hurt.

Completed findings and explicit limitations are in `result.md`.

## Reference alignment

Implementation choices were checked against the local Zhang and Mirarab
(2022) paper in `references/wastral.pdf` and the current official documentation:

- [ASTER weighted ASTRAL tutorial](https://github.com/chaoszhang/ASTER/blob/master/tutorial/wastral.md)
- [ASTER ASTRAL-IV tutorial](https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral4.md)
- [IQ-TREE substitution models](https://iqtree.github.io/doc/Substitution-Models)
- [IQ-TREE AliSim](https://iqtree.github.io/doc/AliSim)
- [SimPhy source and manual links](https://github.com/adamallo/SimPhy)

These references support the five IQ-TREE analysis conditions (one correct-
model reference and four misspecified conditions), explicit AliSim parameters,
ASTRAL input polytomies, the 0.90-aBayes contraction comparison used for S200,
wASTRAL support mode 2, default hybrid weighting, and `-B` local-Bayesian
scaling.
