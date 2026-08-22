from dataclasses import dataclass

from cse6406.inference_models.base import AnalysisModel


@dataclass(frozen=True)
class ExperimentCondition:
    """A single ILS x replicate x sequence-length x model cell."""

    ils: str
    replicate: int
    sequence_length: int
    model: AnalysisModel

    @property
    def key(self) -> str:
        return (
            f"ils={self.ils}|rep={self.replicate:02d}|"
            f"length={self.sequence_length}|model={self.model.name}"
        )
