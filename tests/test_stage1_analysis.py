import unittest

from cse6406.domain.condition import ExperimentCondition
from cse6406.inference_models.base import AnalysisModel
from cse6406.stage1.branch_scorer import BranchObservation, BranchScorer
from cse6406.stage1.calibration import CalibrationAnalyzer, wilson
from cse6406.trees.newick import read_tree_text


class Stage1AnalysisTests(unittest.TestCase):
    def setUp(self):
        self.condition = ExperimentCondition(
            "low", 1, 200, AnalysisModel("gtr_g4", "GTR+G4", "test")
        )

    def test_branch_correctness_uses_true_gene_tree(self):
        truth = read_tree_text("((A,B),(C,(D,E)));")
        estimate = read_tree_text("((A,B)/0.99,(C,(D,E)/0.80)/0.70);")
        rows = BranchScorer().score(self.condition, 1, estimate, truth)
        self.assertTrue(rows)
        self.assertTrue(all(row.correct for row in rows))
        self.assertIn(0.99, {row.support for row in rows})

    def test_summary_reports_conditional_high_support_error(self):
        common = dict(
            condition=self.condition.key, replicate=1, locus=1, ils="low",
            sequence_length=200, analysis_model="gtr_g4", split_size=2,
        )
        observations = [
            BranchObservation(**common, support=.99, correct=True),
            BranchObservation(**common, support=.98, correct=False),
            BranchObservation(**common, support=.50, correct=False),
            BranchObservation(**common, support=None, correct=False),
        ]
        branches, calibration, summary, replicate_summary = CalibrationAnalyzer().analyze(
            observations
        )
        self.assertEqual(len(branches), 4)
        self.assertTrue(calibration)
        self.assertEqual(summary[0]["n_missing_support"], 1)
        self.assertEqual(summary[0]["n_supported_branches"], 3)
        self.assertAlmostEqual(summary[0]["branch_error_rate"], 3 / 4)
        self.assertAlmostEqual(summary[0]["p_wrong_given_support_ge_0.95"], 0.5)
        self.assertIsNotNone(summary[0]["brier_score"])
        self.assertIsNotNone(summary[0]["expected_calibration_error"])
        self.assertIn("mean_support", calibration[0])
        self.assertEqual(replicate_summary[0]["replicate"], 1)

    def test_wilson_handles_boundary(self):
        low, high = wilson(10, 10)
        self.assertLess(low, 1)
        self.assertEqual(high, 1)


if __name__ == "__main__":
    unittest.main()
