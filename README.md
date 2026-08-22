# Model misspecification, branch support, and weighted ASTRAL

This repository implements the controlled experiment in `overall_plan.md`:

```text
SimPhy low/high ILS true histories
        -> AliSim GTR+G4 DNA
        -> IQ-TREE under GTR+G4, GTR, HKY+G4, HKY, and JC
        -> Stage 1 branch correctness and support calibration
        -> ASTRAL-IV versus support-weighted wASTRAL
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
A smaller pilot command is provided in `INSTALL.md`.

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
`P(wrong | support >= 0.95)`, `P(wrong | support >= 0.99)`, missing support,
and replicate-level summaries. ILS is checked independently as true-gene-tree
versus true-species-tree normalized RF discordance.

Stage 2 runs unweighted ASTRAL-IV and `wastral --mode 2 -B` on the same inferred
gene trees and reports

```text
delta = weighted normalized RF error - unweighted normalized RF error
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

These references support the five IQ-TREE models, discrete-gamma syntax,
AliSim tree-based simulation, ASTRAL topology input, wASTRAL support mode 2,
and `-B` scaling for aBayes/local Bayesian support.
