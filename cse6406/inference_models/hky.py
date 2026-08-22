"""Misspecification option: HKY without gamma-distributed rate heterogeneity."""

from .base import AnalysisModel

MODEL = AnalysisModel(
    name="hky",
    iqtree_model="HKY",
    role="Combines a simpler substitution model with removal of rate heterogeneity.",
)
