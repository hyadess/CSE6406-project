# Installation and reproducible setup

This file contains the complete setup for the full experiment in
`overall_plan.md`. The production workflow requires Python, SimPhy, IQ-TREE,
ASTRAL-IV, and wASTRAL. Stage 1 can run without the two ASTER executables.

The implementation was smoke-tested on macOS ARM64 with Python 3.13.15,
DendroPy 5.0.10, IQ-TREE 2.4.0, and SimPhy 1.0.2.

## 1. Create the Python environment

From the repository root:

```bash
conda create -n cse6406 -c conda-forge python=3.13 pip
conda activate cse6406
python -m pip install -r requirements.txt
```

Verify Python and DendroPy:

```bash
python --version
python -c "import dendropy; print(dendropy.__version__)"
```

## 2. Install IQ-TREE with AliSim

The validated version is IQ-TREE 2.4.0:

```bash
conda install -c conda-forge -c bioconda iqtree=2.4.0
iqtree2 --version
```

IQ-TREE performs both AliSim sequence generation and ML gene-tree inference
with aBayes support. The fixed analysis grid is `GTR+G4`, `GTR`, `HKY+G4`,
`HKY`, and `JC`.

AliSim generates data with
`GTR{1,2,1,1,2}+F{0.30,0.20,0.20,0.30}+G4{0.20}`. Here `+F{...}` fixes
unequal A/C/G/T frequencies at 0.30/0.20/0.20/0.30; it is useful because it
does not force the equal-frequency assumption used by JC. `+G4{0.20}` models
strong variation in evolutionary speed among sites using four gamma
categories. See `README.md` for the complete plain-language explanation.

## 3. Build SimPhy

The production command conditions the ingroup tree on 50 leaves and a height
of 2,500,000 generations, adds one outgroup, and samples one individual per
species. `-su e:10000000` is an exponential substitution-rate distribution
with mean `1e-7`; SimPhy uses the reciprocal mean as the exponential argument.
The pipeline stores this scientifically meaningful mean as
`substitution_rate_mean` and converts it to SimPhy syntax.

Install native libraries in the same Conda environment:

```bash
conda install -c conda-forge gsl mpfr gmp sqlite
```

Clone SimPhy beside this repository and compile it:

```bash
cd ..
git clone https://github.com/adamallo/SimPhy.git
cd SimPhy
mv bin/simphy bin/simphy.prebuilt
make \
  CFLAGS="-I${CONDA_PREFIX}/include" \
  LDFLAGS="-L${CONDA_PREFIX}/lib -Wl,-rpath,${CONDA_PREFIX}/lib"
cd ../CSE6406-project
```

On macOS, confirm that the result matches the machine architecture:

```bash
file ../SimPhy/bin/simphy
otool -L ../SimPhy/bin/simphy
../SimPhy/bin/simphy | head -3
```

The prebuilt file in the upstream repository may target a different platform;
the native rebuild is therefore intentional.

## 4. Install ASTER (ASTRAL-IV and wASTRAL)

The current ASTER documentation offers Conda installation:

```bash
conda install -c conda-forge -c bioconda aster
astral4 -h
wastral -h
```

If the package does not provide binaries for the platform, build from source:

```bash
cd ..
git clone https://github.com/chaoszhang/ASTER.git
cd ASTER
make
bin/astral4 -h
bin/wastral -h
cd ../CSE6406-project
```

The project uses:

- `astral4 -u 0` for the unweighted topology analysis;
- `astral4 -u 0` after contracting aBayes branches below 0.90;
- `wastral --mode 2 -B -u 0` for support-only weighting with local Bayesian
  support scaled from 0.333 to 1;
- default `wastral -B -u 0` for hybrid support-plus-length weighting.

These choices follow the current official ASTER tutorials. The weighted input
preparer converts IQ-TREE labels such as `/0.994` to numeric values. If
IQ-TREE leaves a zero-length internal branch unlabeled, the preparer records
the event and assigns 1/3, the documented local-Bayesian lower bound and thus
the minimum support weight.

## 5. Configure executable locations

Executables on `PATH` are discovered automatically. Otherwise set absolute
paths:

```bash
export SIMPHY_BIN=/Users/sayemshahad/Downloads/SimPhy/bin/simphy
export IQTREE_BIN=/Users/sayemshahad/anaconda3/envs/cse6406/bin/iqtree2
export ASTRAL_BIN=/Users/sayemshahad/Downloads/ASTER/bin/astral4
export WASTRAL_BIN=/Users/sayemshahad/Downloads/ASTER/bin/wastral
```

Do not edit source files to configure machine-specific paths.

## 6. Verify the installation

Run all Python tests:

```bash
python -m unittest discover -s tests -v
python tests/test_trees.py
```

Preflight the full toolchain:

```bash
python run_project.py --check
```

`preflight.json` records the resolved path, reported version, and SHA-256
digest of every executable. SimPhy 1.0.2 does not report its version through a
standard flag, so its pinned installation instructions plus the executable
digest identify the exact build.

If ASTER is not installed yet, verify the Stage 1 toolchain only:

```bash
python run_project.py --check --stage1-only
```

Run a small real-tool smoke test. This is for integration verification, not a
scientific result:

```bash
python run_project.py \
  --stage1-only \
  --output work/smoke \
  --taxa 6 \
  --loci 2 \
  --replicates 1 \
  --lengths 60 \
  --skip-ils-gate \
  --threads 1
```

## 7. Run the experiment

The preregistered default grid is 51 taxa, 200 loci, 50 replicates, two ILS
levels, three sequence lengths, and five IQ-TREE models:

```bash
python run_project.py --check
python run_project.py --threads AUTO
```

This implies 300,000 IQ-TREE locus analyses
(`2 × 50 × 3 × 5 × 200`) and 1,500 paired species-tree conditions. Run it on an
appropriate workstation or scheduler. The workflow is resumable: validated
trees and alignments are reused on restart.

For a smaller production pilot:

```bash
python run_project.py \
  --output work/pilot \
  --taxa 21 \
  --loci 50 \
  --replicates 3 \
  --lengths 200 800 \
  --threads AUTO
```

For Stage 1 alone:

```bash
python run_project.py --stage1-only --threads AUTO
```

To run only selected inference setups, use the names of their files under
`cse6406/inference_models/`:

```bash
python run_project.py --models gtr_g4 jc --threads AUTO
```

## 8. Output layout

Each ILS/replicate directory contains:

```text
simphy/1/                         true species and per-locus gene trees
alignments/L<length>/             shared AliSim alignments
estimated_gene_trees/L*/<model>/  IQ-TREE outputs for one model
logs/                             exact commands and captured program output
```

`<output>/results/` contains:

```text
ils_verification.csv              per-replicate true-gene/species discordance
ils_summary.csv                   preregistered low/high ILS acceptance checks
stage1_branches.csv               one row per estimated internal branch
stage1_calibration.csv            calibration and mean support by bin
stage1_summary.csv                branch error, Brier/ECE, high-support errors
stage1_replicate_summary.csv      replicate-level Stage 1 summaries
stage2_species_tree_error.csv     four methods, nRF errors, and paired deltas
```

`design.json`, `preflight.json`, per-command logs, and support-imputation logs
make every run auditable.

## 9. Clean-project boundary

Only `run_project.py` and `cse6406/` form the active pipeline. Superseded code,
old generated outputs, and historical tests are isolated in
`archive/previous_pipeline/` and are not part of installation or validation.
