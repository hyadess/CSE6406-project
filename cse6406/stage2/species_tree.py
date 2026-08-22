from pathlib import Path

from cse6406.trees.newick import read_tree, topology_distance


class SpeciesTreeEvaluator:
    """Compare paired ASTRAL/wASTRAL outputs with the SimPhy species truth."""

    def evaluate(
        self, *, condition, true_species_tree: Path,
        astral_tree: Path, wastral_tree: Path,
    ) -> dict[str, object]:
        truth = read_tree(true_species_tree)
        unweighted_error = topology_distance(read_tree(astral_tree), truth)
        weighted_error = topology_distance(read_tree(wastral_tree), truth)
        return {
            "condition": condition.key, "ils": condition.ils,
            "replicate": condition.replicate,
            "sequence_length": condition.sequence_length,
            "analysis_model": condition.model.name,
            "unweighted_nrf": unweighted_error,
            "weighted_nrf": weighted_error,
            "delta_weighted_minus_unweighted": weighted_error - unweighted_error,
        }
