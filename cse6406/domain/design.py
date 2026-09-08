from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from cse6406.inference_models.base import AnalysisModel
from cse6406.inference_models.registry import ANALYSIS_MODELS

from .condition import ExperimentCondition


@dataclass(frozen=True)
class ILSLevel:
    name: str
    population_size: int


@dataclass(frozen=True)
class ExperimentalDesign:
    """The complete factorial design and reproducibility parameters."""

    output: Path = Path("work/full_experiment_r50")
    taxa: int = 51
    loci: int = 200
    replicates: int = 50
    sequence_lengths: tuple[int, ...] = (200, 800, 1600)
    ils_levels: tuple[ILSLevel, ...] = (
        ILSLevel("low", 100_000),
        ILSLevel("high", 1_000_000),
    )
    models: tuple[AnalysisModel, ...] = ANALYSIS_MODELS
    seed: int = 6_406_001
    threads: str = "1"
    birth_rate: float = 1e-7
    tree_height: int = 2_500_000
    substitution_rate_mean: float = 1e-7
    generating_model: str = "GTR{1,2,1,1,2}+F{0.30,0.20,0.20,0.30}+G4{0.20}"
    contracted_abayes_threshold: float = 0.90
    low_ils_max_nrf: float = 0.25
    high_ils_min_nrf: float = 0.55
    minimum_ils_nrf_gap: float = 0.30
    enforce_ils_gate: bool = True

    def validate(self) -> None:
        if self.taxa < 4 or self.loci < 1 or self.replicates < 1:
            raise ValueError("taxa >= 4, loci >= 1, and replicates >= 1 are required")
        if any(length < 50 for length in self.sequence_lengths):
            raise ValueError("every sequence length must be at least 50")
        if len({x.name for x in self.ils_levels}) != len(self.ils_levels):
            raise ValueError("ILS level names must be unique")
        if len({x.name for x in self.models}) != len(self.models):
            raise ValueError("analysis-model names must be unique")
        if self.substitution_rate_mean <= 0:
            raise ValueError("substitution_rate_mean must be positive")
        if not 1 / 3 <= self.contracted_abayes_threshold <= 1:
            raise ValueError("contracted_abayes_threshold must be between 1/3 and 1")
        if not 0 <= self.low_ils_max_nrf < self.high_ils_min_nrf <= 1:
            raise ValueError("ILS nRF thresholds must satisfy 0 <= low < high <= 1")
        if not 0 < self.minimum_ils_nrf_gap <= 1:
            raise ValueError("minimum_ils_nrf_gap must be in (0,1]")

    def conditions(self) -> list[ExperimentCondition]:
        return [
            ExperimentCondition(ils.name, replicate, length, model)
            for ils in self.ils_levels
            for replicate in range(1, self.replicates + 1)
            for length in self.sequence_lengths
            for model in self.models
        ]

    def ils_level(self, name: str) -> ILSLevel:
        return next(level for level in self.ils_levels if level.name == name)

    def replicate_dir(self, ils: str, replicate: int) -> Path:
        return self.output / f"ils_{ils}" / f"replicate_{replicate:02d}"

    def alignment_dir(self, ils: str, replicate: int, length: int) -> Path:
        return self.replicate_dir(ils, replicate) / "alignments" / f"L{length}"

    def model_dir(self, condition: ExperimentCondition) -> Path:
        return (
            self.replicate_dir(condition.ils, condition.replicate)
            / "estimated_gene_trees" / f"L{condition.sequence_length}"
            / condition.model.name
        )
