#!/usr/bin/env python3
"""Run or preflight the complete experiment defined in overall_plan.md."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from cse6406.core.errors import PipelineError
from cse6406.core.toolchain import Toolchain
from cse6406.domain.design import ExperimentalDesign
from cse6406.inference_models.registry import ANALYSIS_MODELS, MODEL_BY_NAME
from cse6406.pipeline import ExperimentPipeline


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("work/full_experiment"))
    parser.add_argument("--taxa", type=int, default=51)
    parser.add_argument("--loci", type=int, default=200)
    parser.add_argument("--replicates", type=int, default=10)
    parser.add_argument("--lengths", type=int, nargs="+", default=[200, 800, 1600])
    parser.add_argument(
        "--models", nargs="+", choices=tuple(MODEL_BY_NAME),
        default=[model.name for model in ANALYSIS_MODELS],
        help="IQ-TREE model files to run; defaults to the complete five-model grid",
    )
    parser.add_argument("--seed", type=int, default=6_406_001)
    parser.add_argument("--threads", default="1")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--stage1-only", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    design = ExperimentalDesign(
        output=args.output.resolve(), taxa=args.taxa, loci=args.loci,
        replicates=args.replicates, sequence_lengths=tuple(args.lengths),
        models=tuple(MODEL_BY_NAME[name] for name in args.models),
        seed=args.seed, threads=args.threads,
    )
    pipeline = ExperimentPipeline(design, Toolchain.discover())
    if args.check:
        for name, version in pipeline.preflight(stage2=not args.stage1_only).items():
            print(f"[ok] {name}: {version}")
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
