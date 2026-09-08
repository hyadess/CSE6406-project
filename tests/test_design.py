import unittest
from pathlib import Path

from cse6406.domain.design import ExperimentalDesign


class ExperimentalDesignTests(unittest.TestCase):
    def test_publication_defaults_use_fifty_replicates(self):
        design = ExperimentalDesign()
        self.assertEqual(design.replicates, 50)
        self.assertEqual(design.substitution_rate_mean, 1e-7)

    def test_full_grid_contains_every_factor(self):
        design = ExperimentalDesign(
            output=Path("unused"), replicates=2, sequence_lengths=(200, 800)
        )
        conditions = design.conditions()
        self.assertEqual(len(conditions), 2 * 2 * 2 * 5)
        self.assertEqual(
            {condition.model.iqtree_model for condition in conditions},
            {"GTR+G4", "GTR", "HKY+G4", "HKY", "JC"},
        )
        self.assertEqual({condition.ils for condition in conditions}, {"low", "high"})

    def test_model_directories_are_unique(self):
        design = ExperimentalDesign(output=Path("work/test"), replicates=1)
        paths = [design.model_dir(condition) for condition in design.conditions()]
        self.assertEqual(len(paths), len(set(paths)))


if __name__ == "__main__":
    unittest.main()
