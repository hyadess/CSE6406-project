#!/usr/bin/env python3
"""Run the experiment from overall_plan.md across parallel worker processes.

This is the multi-process counterpart to ``run_project.py``. Every experimental
parameter, seed, output path, and result table is unchanged; only the execution
strategy differs. The two runners share the same output directory format and the
same ``design.json`` contract, so a run started by one can be resumed by the
other.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from cse6406.core.errors import PipelineError
from cse6406.core.toolchain import Toolchain
from cse6406.domain.design import ExperimentalDesign
from cse6406.inference_models.registry import ANALYSIS_MODELS, MODEL_BY_NAME
from cse6406.parallel.pipeline import ParallelExperimentPipeline, default_workers


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--output", type=Path, default=Path("work/full_experiment_r50"))
    parser.add_argument("--taxa", type=int, default=51)
    parser.add_argument("--loci", type=int, default=200)
    parser.add_argument("--replicates", type=int, default=50)
    parser.add_argument("--lengths", type=int, nargs="+", default=[200, 800, 1600])
    parser.add_argument(
        "--models", nargs="+", choices=tuple(MODEL_BY_NAME),
        default=[model.name for model in ANALYSIS_MODELS],
        help="IQ-TREE model files to run; defaults to the complete five-model grid",
    )
    parser.add_argument("--seed", type=int, default=6_406_001)
    parser.add_argument(
        "--threads", default="1",
        help="recorded in design.json for provenance; see --tool-threads for the "
             "threads actually handed to each external tool",
    )
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--stage1-only", action="store_true")
    parser.add_argument(
        "--skip-ils-gate", action="store_true",
        help="allow diagnostic smoke/pilot runs outside the preregistered ILS ranges",
    )
    parser.add_argument(
        "--workers", type=int, default=None,
        help=f"concurrent worker processes (default: detected CPUs, {default_workers()} here)",
    )
    parser.add_argument(
        "--tool-threads", default=None,
        help="threads given to each IQ-TREE/ASTER invocation (default: 1 when running "
             "more than one worker, otherwise the --threads value)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    design = ExperimentalDesign(
        output=args.output.resolve(), taxa=args.taxa, loci=args.loci,
        replicates=args.replicates, sequence_lengths=tuple(args.lengths),
        models=tuple(MODEL_BY_NAME[name] for name in args.models),
        seed=args.seed, threads=args.threads, enforce_ils_gate=not args.skip_ils_gate,
    )
    workers = default_workers() if args.workers is None else max(1, args.workers)
    tool_threads = args.tool_threads
    if tool_threads is None:
        # Parallelism already comes from running many independent tools at once,
        # so each tool takes a single thread unless the caller overrides it.
        tool_threads = design.threads if workers == 1 else "1"

    pipeline = ParallelExperimentPipeline(
        design, Toolchain.discover(), workers=workers, tool_threads=tool_threads,
    )
    if args.check:
        for name, metadata in pipeline.preflight(stage2=not args.stage1_only).items():
            print(f"[ok] {name}: {metadata['version']} ({metadata['path']})")
        print(f"[ok] workers: {workers}; threads per external tool: {tool_threads}")
        return 0
    pipeline.run(stage2=not args.stage1_only)
    print(f"Completed results: {design.output / 'results'}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (PipelineError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
