from dataclasses import dataclass


@dataclass(frozen=True)
class AnalysisModel:
    """A named, fixed IQ-TREE analysis option."""

    name: str
    iqtree_model: str
    role: str
