from __future__ import annotations

from dataclasses import dataclass

from cse6406.trees.newick import informative_splits


# ``slots`` keeps the per-branch record small: the production design holds
# roughly 14 million of these in memory at once for the Stage 1 tables.
@dataclass(frozen=True, slots=True)
class BranchObservation:
    condition: str
    replicate: int
    locus: int
    ils: str
    sequence_length: int
    analysis_model: str
    support: float | None
    correct: bool
    split_size: int


class BranchScorer:
    """Score estimated internal branches against the locus's true gene tree."""

    def score(self, condition, locus: int, estimated, truth) -> list[BranchObservation]:
        true_splits = set(informative_splits(truth))
        return [
            BranchObservation(
                condition=condition.key, replicate=condition.replicate, locus=locus,
                ils=condition.ils, sequence_length=condition.sequence_length,
                analysis_model=condition.model.name, support=support,
                correct=split in true_splits, split_size=len(split),
            )
            for split, support in informative_splits(estimated, with_support=True)
        ]
