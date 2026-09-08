import unittest

from cse6406.core.errors import PipelineError
from cse6406.stage1.ils import ILSAnalyzer


class ILSTreatmentTests(unittest.TestCase):
    def test_preregistered_ranges_are_required(self):
        rows = [
            {"ils": "low", "n_loci": 10, "mean_true_gene_species_nrf": 0.17},
            {"ils": "high", "n_loci": 10, "mean_true_gene_species_nrf": 0.69},
        ]
        summary = ILSAnalyzer().summarize(rows)
        self.assertTrue(summary[0]["meets_preregistered_ils_criteria"])
        ILSAnalyzer.require_acceptable(summary)

    def test_small_difference_is_rejected(self):
        rows = [
            {"ils": "low", "n_loci": 10, "mean_true_gene_species_nrf": 0.40},
            {"ils": "high", "n_loci": 10, "mean_true_gene_species_nrf": 0.45},
        ]
        summary = ILSAnalyzer().summarize(rows)
        self.assertFalse(summary[0]["meets_preregistered_ils_criteria"])
        with self.assertRaises(PipelineError):
            ILSAnalyzer.require_acceptable(summary)


if __name__ == "__main__":
    unittest.main()
