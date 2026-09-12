# Parallel Runbook

`run_project_parallel.py` runs the experiment defined in `overall_plan.md`
across many worker processes instead of one at a time. It is a drop-in
replacement for `run_project.py`: the design, the seeds, the directory layout,
the resume rules, and every result table are unchanged.

## Why the sequential runner is slow

`run_project.py` issues every external command one after another in a single
process. The production grid is 100 SimPhy replicates, 60,000 AliSim
alignments, 300,000 IQ-TREE gene-tree analyses, and 1,500 conditions of
ASTRAL/wASTRAL. Each of those commands is small, so a single core spends the
whole run starting one short program at a time while the other cores idle.

The units of work are already independent — a locus never reads another locus's
output — so they can run concurrently without changing a single result.

## What runs in parallel

| Phase | Task granularity | Tasks in the production design |
|---|---|---|
| SimPhy true trees | one replicate | 100 |
| ILS gate verification | one replicate | 100 |
| AliSim alignments | one locus | 60,000 |
| IQ-TREE inference + branch scoring | one locus | 300,000 |
| ASTRAL / wASTRAL species trees | one condition | 1,500 |

Branch scoring runs inside the same worker as the inference that produced the
tree, so the Newick parsing for 300,000 loci is parallel too rather than a
single-threaded pass at the end.

Phases stay in the original order, so the preregistered ILS gate still stops the
run before any alignment or inference work when the treatments miss their nRF
criteria.

## Results are identical, not merely equivalent

Every seed is a pure function of the design and the replicate/length/locus
indices, so concurrency cannot change what any tool computes. Results are
collected back in submission order, so every CSV keeps the same row order the
sequential runner produces.

This is verified end to end: for a complete small design the two runners produce
byte-identical `results/*.csv`, byte-identical alignments, gene trees, and
species trees.

```bash
diff -r work/serial_run/results work/parallel_run/results
```

## Measured speedup

Benchmarked on this project's 10-core machine using the paper-subset shape
(`--taxa 21 --loci 10 --replicates 2 --lengths 200 800`, all five models: 40
conditions, 400 gene-tree analyses). Every row below produced results identical
to the sequential run.

| Runner | Wall clock | Speedup |
|---|---|---|
| `run_project.py` (sequential) | 268.5 s | 1.00x |
| `run_project_parallel.py --workers 4` | 123.2 s | 2.18x |
| `run_project_parallel.py --workers 6` | 96.6 s | 2.78x |
| `run_project_parallel.py --workers 8` | 64.2 s | 4.18x |
| `run_project_parallel.py --workers 10` | 60.1 s | 4.47x |

Speedup is below the core count because not all cores are equal on this machine
and because sustained all-core load reduces clock speed. More workers still won
at every step, so the default of one worker per detected CPU is the right
starting point; lower it only to leave room for other work.

The production grid has far more work per task than this subset - 51 taxa and up
to 1,600 bp instead of 21 taxa and 200-800 bp - so per-task overhead matters
proportionally less there.

## Threads versus workers

Two different knobs control CPU use, and only one of them should be large.

- `--workers N` — how many tasks run at once. Defaults to the detected CPU
  count. This is where the speedup comes from.
- `--tool-threads T` — threads handed to each individual IQ-TREE/ASTER
  invocation. Defaults to `1` whenever more than one worker is used.

Total cores in use is roughly `N x T`, so leave `--tool-threads` at `1`.
Internal threading is a poor fit here anyway: these alignments are 200-1,600 bp,
and IQ-TREE cannot keep several threads busy on a tree that small. Running many
single-threaded analyses at once scales far better than running one analysis
with many threads.

`--threads` is kept only so `design.json` stays identical to a sequential run
and the two runners can resume each other's output directories. It does not
control the pool. With `--workers 1` the runner reproduces sequential behaviour
exactly, including passing `--threads` through to the tools.

## Memory

The Stage 1 branch table is the largest output — roughly 14 million rows for the
production design. The parallel runner streams it to
`results/stage1_branches.csv` as workers finish instead of building a second
full copy in memory, which keeps peak memory far below the sequential runner's
on the same design.

## Running it

Set up the environment and tool paths exactly as in `RUN_EXPERIMENT.md`
sections 1-3, then substitute this runner.

Preflight, including the resolved worker and thread counts:

```bash
python run_project_parallel.py --check --output work/full_experiment_r50 --taxa 51 --loci 200 --replicates 50 --lengths 200 800 1600 --models gtr_g4 gtr hky_g4 hky jc --seed 6406001 --threads AUTO
```

Full production grid:

```bash
python run_project_parallel.py --output work/full_experiment_r50 --taxa 51 --loci 200 --replicates 50 --lengths 200 800 1600 --models gtr_g4 gtr hky_g4 hky jc --seed 6406001 --threads AUTO
```

Unattended:

```bash
nohup python run_project_parallel.py --output work/full_experiment_r50 --taxa 51 --loci 200 --replicates 50 --lengths 200 800 1600 --models gtr_g4 gtr hky_g4 hky jc --seed 6406001 --threads AUTO > work/full_experiment_r50/run.log 2>&1 &
echo $! > work/full_experiment_r50/run.pid
```

Paper subset, matching `RUN_PAPER_SUBSET.md`:

```bash
python run_project_parallel.py --skip-ils-gate --output work/paper_subset_10x2 --taxa 21 --loci 10 --replicates 2 --lengths 200 800 1600 --models gtr_g4 gtr hky_g4 hky jc --seed 6406001 --threads AUTO
```

Leave cores free for other work with `--workers`:

```bash
python run_project_parallel.py --workers 8 ...
```

## Monitoring

Each phase prints its own progress line with completion count, elapsed time,
ETA, and throughput:

```
[iqtree] 41200/300000 (13.7%) elapsed 02:41:07 eta 16d 21:03:44 rate 4.26/s
```

Every `find | wc -l` progress check in `RUN_EXPERIMENT.md` section 5 still
applies unchanged.

## Interrupting and resuming

Resume with the identical command; completed artifacts are detected on disk and
skipped, exactly as in the sequential runner. A run started with
`run_project.py` can be resumed with `run_project_parallel.py` and vice versa,
because both write and validate the same `design.json`.

If a worker fails, the pool cancels the queued tasks and the error propagates
with the same `PipelineError` message and log path as the sequential runner.
