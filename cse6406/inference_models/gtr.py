"""Misspecification option: GTR substitution model without gamma rates."""

from .base import AnalysisModel

MODEL = AnalysisModel(
    name="gtr",
    iqtree_model="GTR",
    role="Removes among-site rate heterogeneity while retaining GTR exchangeabilities.",
)
