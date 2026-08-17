# Installation

Validated on macOS ARM64 with Python 3, DendroPy 5.0.10, IQ-TREE 2.4.0 and a
native ARM64 build of SimPhy.

## Requirements by experiment

| Experiment | Requirements |
|---|---|
| Part A — S100 | Python packages and S100 trees |
| Part B — simulation | Python packages, IQ-TREE and SimPhy |
| Part C — avian survey | Python packages only; data are bundled |

## 1. Python

```bash
conda create -n cse6406 python=3.11
conda activate cse6406
python -m pip install -r requirements.txt
```

Verify:

```bash
python -c "import dendropy; print(dendropy.__version__)"
python tests/test_trees.py
```

## 2. IQ-TREE 2.4.0

Install the version used to validate this project:

```bash
conda install -c conda-forge -c bioconda iqtree=2.4.0
which iqtree2
iqtree2 --version
```

IQ-TREE supplies both AliSim and maximum-likelihood inference with aBayes
support.

## 3. SimPhy on macOS ARM64

Install build libraries into the active Conda environment:

```bash
conda install -c conda-forge gsl mpfr gmp sqlite
```

From the directory containing this project, clone and build SimPhy:

```bash
git clone https://github.com/adamallo/SimPhy.git
cd SimPhy
mv bin/simphy bin/simphy.prebuilt

make \
  CFLAGS="-I${CONDA_PREFIX}/include" \
  LDFLAGS="-L${CONDA_PREFIX}/lib -Wl,-rpath,${CONDA_PREFIX}/lib"
```

Verify the native binary and its libraries:

```bash
file bin/simphy
otool -L bin/simphy
./bin/simphy | head -3
```

`src/config.py` discovers `iqtree2` from `PATH` and uses
`../SimPhy/bin/simphy` as the default SimPhy location. Either path can be
overridden:

```bash
export IQTREE_BIN=/path/to/iqtree2
export SIMPHY_BIN=/path/to/simphy
```

Return to the project and check the toolchain:

```bash
python src/run_simulated.py --check
```

For Linux or Intel macOS, install the matching IQ-TREE binary and build SimPhy
with the same required libraries. The environment-variable overrides avoid any
need to edit `src/config.py`.

## 4. S100 data for Part A

Download and extract the supported estimated trees and species trees:

```bash
mkdir -p data/downloads

curl -L --fail --retry 3 -C - \
  -o data/downloads/S100.zip \
  https://datadryad.org/downloads/file_stream/2374124

unzip -q data/downloads/S100.zip -d data/downloads/S100_dryad
```

Download the original archive containing the true gene trees:

```bash
curl -L --fail --retry 3 -C - \
  -o data/downloads/S101_gitlab.tar.gz \
  'https://gitlab.com/esayyari/ASTRALIII/-/raw/master/S101.tar.gz?inline=false'

tar -xOzf data/downloads/S101_gitlab.tar.gz genetrees.tar.gz \
  | tar -xzf - -C data/downloads/S100_dryad {01..50}/truegenetrees
```

Verify the merged dataset:

```bash
find data/downloads/S100_dryad -type f -name truegenetrees | wc -l
wc -l data/downloads/S100_dryad/01/truegenetrees

python src/run_s100.py \
  --check \
  --s100-dir data/downloads/S100_dryad \
  --reps 01 \
  --lengths 200 \
  --support .abayes
```

Expected: 50 `truegenetrees` files, 1,000 trees per replicate, 101 taxa and
aBayes support in `[0.333, 1.0]`.

The nested `alignments.tar.gz` archive is not required for the current Stage 1
or Stage 2 experiments.

## 5. Final checks

```bash
python tests/test_trees.py
python src/run_simulated.py --check
python src/run_s100.py \
  --check \
  --s100-dir data/downloads/S100_dryad \
  --reps 01 \
  --lengths 200
```
