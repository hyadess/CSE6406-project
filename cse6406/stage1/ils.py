from pathlib import Path

from cse6406.core.errors import PipelineError
from cse6406.trees.newick import read_tree, topology_distance


class ILSAnalyzer:
    """Verify intended ILS using true-gene versus true-species discordance."""

    def analyze(
        self, *, ils: str, replicate: int, species_tree: Path,
        gene_trees: list[Path], population_size: int,
    ) -> dict[str, object]:
        species = read_tree(species_tree)
        distances = [topology_distance(read_tree(path), species) for path in gene_trees]
        return {
            "ils": ils, "replicate": replicate, "population_size": population_size,
            "n_loci": len(distances),
            "mean_true_gene_species_nrf": sum(distances) / len(distances),
            "min_true_gene_species_nrf": min(distances),
            "max_true_gene_species_nrf": max(distances),
        }

    def summarize(
        self, rows: list[dict[str, object]], *,
        low_max: float = 0.25, high_min: float = 0.55, minimum_gap: float = 0.30,
    ) -> list[dict[str, object]]:
        output = []
        means = {}
        for ils in sorted({str(row["ils"]) for row in rows}):
            selected = [row for row in rows if row["ils"] == ils]
            total_loci = sum(int(row["n_loci"]) for row in selected)
            pooled = sum(
                float(row["mean_true_gene_species_nrf"]) * int(row["n_loci"])
                for row in selected
            ) / total_loci
            means[ils] = pooled
            output.append({
                "ils": ils, "replicates": len(selected), "n_loci": total_loci,
                "mean_true_gene_species_nrf": pooled,
            })
        gap = means.get("high", float("nan")) - means.get("low", float("nan"))
        acceptable = (
            "low" in means and "high" in means
            and means["low"] <= low_max
            and means["high"] >= high_min
            and gap >= minimum_gap
        )
        for row in output:
            row["low_max_nrf"] = low_max
            row["high_min_nrf"] = high_min
            row["minimum_nrf_gap"] = minimum_gap
            row["observed_high_minus_low_nrf"] = gap
            row["meets_preregistered_ils_criteria"] = acceptable
        return output

    @staticmethod
    def require_acceptable(summary: list[dict[str, object]]) -> None:
        if not summary or not bool(summary[0]["meets_preregistered_ils_criteria"]):
            raise PipelineError(
                "Simulated ILS treatments failed the preregistered nRF criteria. "
                "Inspect results/ils_summary.csv before running sequence inference."
            )
