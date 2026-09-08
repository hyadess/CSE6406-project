from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import asdict

from .branch_scorer import BranchObservation


BINS = ((0, .5), (.5, .7), (.7, .9), (.9, .95), (.95, .99), (.99, 1.0000001))


def wilson(successes: int, total: int, z: float = 1.96) -> tuple[float, float]:
    if total == 0:
        return float("nan"), float("nan")
    proportion = successes / total
    denominator = 1 + z * z / total
    center = (proportion + z * z / (2 * total)) / denominator
    half = z * math.sqrt(
        proportion * (1 - proportion) / total + z * z / (4 * total * total)
    ) / denominator
    return max(0, center - half), min(1, center + half)


class CalibrationAnalyzer:
    """Produce Stage 1 branch, calibration, and high-support error tables."""

    def analyze(self, observations: list[BranchObservation]):
        supported = [row for row in observations if row.support is not None]
        calibration_groups = defaultdict(list)
        condition_groups = defaultdict(list)
        replicate_groups = defaultdict(list)
        for row in observations:
            condition_groups[self._pooled_key(row)].append(row)
            replicate_groups[row.condition].append(row)
        for row in supported:
            support_bin = next(
                f"[{low:g},{min(high, 1):g}{']' if high > 1 else ')'}"
                for low, high in BINS if low <= row.support < high
            )
            calibration_groups[(self._pooled_key(row), support_bin)].append(row)

        calibration = []
        for (condition, support_bin), rows in sorted(calibration_groups.items()):
            correct = sum(row.correct for row in rows)
            low, high = wilson(correct, len(rows))
            calibration.append({
                "condition": condition, "support_bin": support_bin,
                "n_branches": len(rows), "n_correct": correct,
                "mean_support": sum(float(row.support) for row in rows) / len(rows),
                "p_correct": correct / len(rows), "ci_low": low, "ci_high": high,
            })

        summary = self._summaries(condition_groups, pooled=True)
        replicate_summary = self._summaries(replicate_groups, pooled=False)
        return [asdict(row) for row in observations], calibration, summary, replicate_summary

    @staticmethod
    def _pooled_key(row: BranchObservation) -> str:
        return (
            f"ils={row.ils}|length={row.sequence_length}|"
            f"model={row.analysis_model}"
        )

    @staticmethod
    def _summaries(groups, *, pooled: bool):
        summary = []
        for condition, rows in sorted(groups.items()):
            wrong = [row for row in rows if not row.correct]
            high95 = [row for row in rows if row.support is not None and row.support >= .95]
            high99 = [row for row in rows if row.support is not None and row.support >= .99]
            supported = [row for row in rows if row.support is not None]
            brier = (
                sum((float(row.support) - float(row.correct)) ** 2 for row in supported)
                / len(supported) if supported else None
            )
            ece_parts = []
            for low, high in BINS:
                in_bin = [
                    row for row in supported
                    if low <= float(row.support) < high
                ]
                if in_bin:
                    mean_support = sum(float(row.support) for row in in_bin) / len(in_bin)
                    accuracy = sum(row.correct for row in in_bin) / len(in_bin)
                    ece_parts.append(len(in_bin) * abs(mean_support - accuracy))
            ece = sum(ece_parts) / len(supported) if supported else None
            exemplar = rows[0]
            summary.append({
                "condition": condition, "ils": exemplar.ils,
                "replicate": "pooled" if pooled else exemplar.replicate,
                "sequence_length": exemplar.sequence_length,
                "analysis_model": exemplar.analysis_model,
                "n_branches": len(rows),
                "n_supported_branches": sum(row.support is not None for row in rows),
                "n_missing_support": sum(row.support is None for row in rows),
                "branch_error_rate": len(wrong) / len(rows),
                "brier_score": brier,
                "expected_calibration_error": ece,
                "p_wrong_given_support_ge_0.95": (
                    sum(not row.correct for row in high95) / len(high95) if high95 else None
                ),
                "p_wrong_given_support_ge_0.99": (
                    sum(not row.correct for row in high99) / len(high99) if high99 else None
                ),
                "n_support_ge_0.95": len(high95), "n_support_ge_0.99": len(high99),
            })
        return summary
