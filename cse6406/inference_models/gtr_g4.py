"""Correctly specified reference: GTR substitution model plus four gamma rates."""

from .base import AnalysisModel

MODEL = AnalysisModel(
    name="gtr_g4",
    iqtree_model="GTR+G4",
    role="Correct-model reference retaining substitution complexity and rate heterogeneity.",
)
