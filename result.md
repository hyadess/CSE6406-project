# Current project status

## Active implementation

The reorganized project implements the complete experimental design in
`overall_plan.md` through `run_project.py` and the `cse6406/` package.

The five IQ-TREE inference setups are now separate files:

| File | IQ-TREE option | Experimental role |
|---|---|---|
| `cse6406/inference_models/gtr_g4.py` | `GTR+G4` | Correct-model reference |
| `cse6406/inference_models/gtr.py` | `GTR` | Removes gamma rate heterogeneity |
| `cse6406/inference_models/hky_g4.py` | `HKY+G4` | Simpler substitution model |
| `cse6406/inference_models/hky.py` | `HKY` | Simpler model without gamma |
| `cse6406/inference_models/jc.py` | `JC` | Severe misspecification stress test |

`registry.py` imports these five definitions in the planned order. The CLI
accepts the same names through `--models`, making full and reduced grids
explicit and reproducible.

## Verification completed

- The focused active test suite passes, including contracted-ASTRAL, hybrid
  wASTRAL, calibration-score, and ILS-gate tests.
- All active Python modules compile.
- SimPhy 1.0.2 and IQ-TREE 2.4.0 pass Stage 1 preflight.
- A real-tool smoke run exercised low/high ILS and all five inference models.
- The smoke run produced branch observations, calibration, pooled summaries,
  replicate summaries, and empirical ILS separation outputs.
- Resolved and 0.90-aBayes-contracted ASTRAL-IV, support-only and hybrid
  wASTRAL, support preparation, RF scoring, and paired delta direction are
  covered by tests.

The smoke run verifies the code path only; its two loci and 60 sites are not a
scientific sample and are not reported as experimental evidence.

## Production status

The revised 50-replicate production grid has not been run. Earlier partial
artifacts under `work/full_experiment` belong to the superseded 10-replicate
design and must not be mixed with `work/full_experiment_r50`. Therefore, the
project does not
yet claim a Stage 1 model-ordering result or a Stage 2 weighted-versus-
unweighted result.

The current machine has ASTRAL-IV and wASTRAL v1.25.4.8. Preflight now records
each executable's path, reported version, and SHA-256 digest. `INSTALL.md`
contains the installation, preflight, pilot, and production commands.

## Previous pipeline

The superseded implementations, historical tests, and old generated artifacts
are outside the active project under `archive/previous_pipeline/`. They are
preserved only for recovery and provenance and are not imported or referenced
by the current workflow.
