#!/usr/bin/env python3
"""Stage 3: infer gene trees under the correct GTR+F+G4 model family."""

from iqtree_inference import model_stage_main
from settings import CORRECT_MODEL


if __name__ == "__main__":
    model_stage_main(condition="correct", model=CORRECT_MODEL)
