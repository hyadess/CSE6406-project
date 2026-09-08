from pathlib import Path

from cse6406.trees.newick import read_tree, topology_distance


class SpeciesTreeEvaluator:
    """Compare four paired ASTRAL/wASTRAL outputs with the species truth."""

    def evaluate(
        self, *, condition, true_species_tree: Path,
        astral_tree: Path, contracted_astral_tree: Path,
        wastral_tree: Path, hybrid_wastral_tree: Path,
    ) -> dict[str, object]:
        truth = read_tree(true_species_tree)
        unweighted_error = topology_distance(read_tree(astral_tree), truth)
        contracted_error = topology_distance(read_tree(contracted_astral_tree), truth)
        weighted_error = topology_distance(read_tree(wastral_tree), truth)
        hybrid_error = topology_distance(read_tree(hybrid_wastral_tree), truth)
        return {
            "condition": condition.key, "ils": condition.ils,
            "replicate": condition.replicate,
            "sequence_length": condition.sequence_length,
            "analysis_model": condition.model.name,
            "unweighted_nrf": unweighted_error,
            "contracted_unweighted_nrf": contracted_error,
            "weighted_nrf": weighted_error,
            "hybrid_weighted_nrf": hybrid_error,
            "delta_weighted_minus_unweighted": weighted_error - unweighted_error,
            "delta_weighted_minus_contracted": weighted_error - contracted_error,
            "delta_hybrid_minus_unweighted": hybrid_error - unweighted_error,
            "delta_hybrid_minus_contracted": hybrid_error - contracted_error,
        }
