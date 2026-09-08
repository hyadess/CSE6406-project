import tempfile
import unittest
from pathlib import Path

from cse6406.core.errors import PipelineError
from cse6406.domain.condition import ExperimentCondition
from cse6406.inference_models.base import AnalysisModel
from cse6406.stage2.astral import AstralRunner
from cse6406.stage2.contractor import SupportTreeContractor
from cse6406.stage2.hybrid_wastral import HybridWeightedAstralRunner
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
    def test_astral_converts_auto_threads_to_an_integer(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inputs, output = root / "genes.tre", root / "species.tre"
            inputs.write_text("((A,B),(C,D));\n")
            recorder = RecordingRunner()
            AstralRunner(Path("/bin/astral4"), recorder).run(
                inputs, output, threads="AUTO"
            )
            value = recorder.command[recorder.command.index("-t") + 1]
            self.assertGreater(int(value), 0)

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

    def test_wastral_converts_auto_threads_to_an_integer(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inputs, output = root / "genes.tre", root / "species.tre"
            inputs.write_text("((A,B)/0.99,(C,D)/0.80);\n")
            recorder = RecordingRunner()
            WeightedAstralRunner(Path("/bin/wastral"), recorder).run(
                inputs, output, threads="AUTO"
            )
            value = recorder.command[recorder.command.index("-t") + 1]
            self.assertGreater(int(value), 0)

    def test_hybrid_wastral_uses_default_mode_and_bayesian_scale(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inputs, output = root / "genes.tre", root / "species.tre"
            inputs.write_text("((A,B)/0.99,(C,D)/0.80);\n")
            recorder = RecordingRunner()
            HybridWeightedAstralRunner(Path("/bin/wastral"), recorder).run(inputs, output)
            self.assertIn("-B", recorder.command)
            self.assertNotIn("--mode", recorder.command)

    def test_wastral_rejects_missing_support(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "genes.tre"
            path.write_text("((A,B),(C,D));\n")
            with self.assertRaises(PipelineError):
                WeightedAstralRunner._validate_support(path)

    def test_wastral_rejects_support_below_local_bayesian_minimum(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "genes.tre"
            path.write_text("((A,B)/0.20,(C,D)/0.90);\n")
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

    def test_support_preparer_normalizes_iqtree_rounded_minimum(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source, ready = root / "source.tre", root / "ready.tre"
            source.write_text("((A,B)/0.333,(C,D)/0.9);\n")
            count = SupportTreePreparer().prepare(source, ready)
            self.assertEqual(count, 0)
            WeightedAstralRunner._validate_support(ready)
            self.assertIn("0.3333333333", ready.read_text())
            self.assertIn(
                "rounded_minimum_adjustments\t1",
                ready.with_suffix(".support_imputation.txt").read_text(),
            )

    def test_contractor_collapses_low_and_missing_support(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source, contracted = root / "source.tre", root / "contracted.tre"
            source.write_text("(((A,B)/0.95,C),(D,(E,F)/0.80));\n")
            count = SupportTreeContractor().contract(source, contracted, threshold=0.90)
            tree = read_tree(contracted)
            self.assertGreaterEqual(count, 2)
            self.assertEqual(
                {node.taxon.label for node in tree.leaf_node_iter()},
                {"A", "B", "C", "D", "E", "F"},
            )
            self.assertLess(len(list(tree.preorder_node_iter())), 10)

    def test_delta_sign_is_weighted_minus_unweighted(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            truth = root / "t"
            unweighted, contracted = root / "u", root / "c"
            weighted, hybrid = root / "w", root / "h"
            truth.write_text("((A,B),(C,(D,E)));\n")
            unweighted.write_text("((A,C),(B,(D,E)));\n")
            contracted.write_text(unweighted.read_text())
            weighted.write_text(truth.read_text())
            hybrid.write_text(truth.read_text())
            condition = ExperimentCondition(
                "high", 1, 200, AnalysisModel("jc", "JC", "test")
            )
            row = SpeciesTreeEvaluator().evaluate(
                condition=condition, true_species_tree=truth,
                astral_tree=unweighted, contracted_astral_tree=contracted,
                wastral_tree=weighted, hybrid_wastral_tree=hybrid,
            )
            self.assertLess(row["delta_weighted_minus_unweighted"], 0)
            self.assertLess(row["delta_hybrid_minus_contracted"], 0)


if __name__ == "__main__":
    unittest.main()
