import tempfile
import unittest
from pathlib import Path

from cse6406.core.errors import PipelineError
from cse6406.domain.condition import ExperimentCondition
from cse6406.inference_models.base import AnalysisModel
from cse6406.stage2.species_tree import SpeciesTreeEvaluator
from cse6406.stage2.support_preparer import SupportTreePreparer
from cse6406.stage2.wastral import WeightedAstralRunner
from cse6406.trees.newick import read_tree, topology_distance


class RecordingRunner:
    def __init__(self):
        self.command = None

    def run(self, command, **_kwargs):
        self.command = [str(value) for value in command]
        Path(self.command[self.command.index("-o") + 1]).write_text("((A,B),(C,D));\n")
        return ""


class Stage2Tests(unittest.TestCase):
    def test_wastral_uses_support_mode_and_bayesian_scale(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inputs, output = root / "genes.tre", root / "species.tre"
            inputs.write_text("((A,B)/0.99,(C,D)/0.80);\n")
            recorder = RecordingRunner()
            WeightedAstralRunner(Path("/bin/wastral"), recorder).run(inputs, output)
            self.assertIn("--mode", recorder.command)
            self.assertIn("2", recorder.command)
            self.assertIn("-B", recorder.command)

    def test_wastral_rejects_missing_support(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "genes.tre"
            path.write_text("((A,B),(C,D));\n")
            with self.assertRaises(PipelineError):
                WeightedAstralRunner._validate_support(path)

    def test_support_preparer_imputes_documented_minimum(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source, ready = root / "source.tre", root / "ready.tre"
            source.write_text("((A,B),(C,(D,E)/0.9));\n")
            count = SupportTreePreparer().prepare(source, ready)
            self.assertGreater(count, 0)
            WeightedAstralRunner._validate_support(ready)
            self.assertIn("0.3333333333", ready.read_text())
            self.assertEqual(topology_distance(read_tree(source), read_tree(ready)), 0)

    def test_delta_sign_is_weighted_minus_unweighted(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            truth, unweighted, weighted = root / "t", root / "u", root / "w"
            truth.write_text("((A,B),(C,(D,E)));\n")
            unweighted.write_text("((A,C),(B,(D,E)));\n")
            weighted.write_text(truth.read_text())
            condition = ExperimentCondition(
                "high", 1, 200, AnalysisModel("jc", "JC", "test")
            )
            row = SpeciesTreeEvaluator().evaluate(
                condition=condition, true_species_tree=truth,
                astral_tree=unweighted, wastral_tree=weighted,
            )
            self.assertLess(row["delta_weighted_minus_unweighted"], 0)


if __name__ == "__main__":
    unittest.main()
