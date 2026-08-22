"""Severe misspecification option: equal-rate, equal-frequency Jukes-Cantor."""

from .base import AnalysisModel

MODEL = AnalysisModel(
    name="jc",
    iqtree_model="JC",
    role="Severe stress test with equal substitution rates, frequencies, and no gamma.",
)
