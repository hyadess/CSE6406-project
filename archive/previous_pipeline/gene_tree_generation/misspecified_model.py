#!/usr/bin/env python3
"""Stage 4: infer gene trees under the deliberately misspecified JC model."""

from iqtree_inference import model_stage_main
from settings import MISSPECIFIED_MODEL


if __name__ == "__main__":
    model_stage_main(condition="misspecified", model=MISSPECIFIED_MODEL)
