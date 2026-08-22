"""Misspecification option: simpler HKY substitution model with gamma rates."""

from .base import AnalysisModel

MODEL = AnalysisModel(
    name="hky_g4",
    iqtree_model="HKY+G4",
    role="Simplifies substitution rates while retaining among-site rate heterogeneity.",
)
