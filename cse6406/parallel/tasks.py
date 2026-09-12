"""Picklable units of work executed inside pool workers.

Every function here reuses the same classes as the sequential pipeline, so a
parallel run issues byte-identical external commands with byte-identical seeds
and therefore produces the same scientific artifacts.
"""
from __future__ import annotations

from pathlib import Path

from cse6406.core.command_runner import CommandRunner
from cse6406.core.toolchain import Toolchain
from cse6406.domain.condition import ExperimentCondition
from cse6406.domain.design import ExperimentalDesign
from cse6406.gene_trees.iqtree import IQTreeInferer
from cse6406.simulation.alisim import AliSimSimulator
from cse6406.simulation.simphy import SimPhySimulator
from cse6406.stage1.branch_scorer import BranchObservation, BranchScorer
from cse6406.stage1.ils import ILSAnalyzer
from cse6406.stage2.astral import AstralRunner
from cse6406.stage2.contractor import SupportTreeContractor
from cse6406.stage2.hybrid_wastral import HybridWeightedAstralRunner
from cse6406.stage2.species_tree import SpeciesTreeEvaluator
from cse6406.stage2.support_preparer import SupportTreePreparer
from cse6406.stage2.wastral import WeightedAstralRunner
from cse6406.trees.newick import combine_newick, read_tree


class WorkerContext:
    """Constant per-run state installed once in each worker process."""

    def __init__(self, tools: Toolchain, design: ExperimentalDesign, tool_threads: str):
        self.tools, self.design, self.tool_threads = tools, design, tool_threads
        self.runner = CommandRunner()
        self.scorer = BranchScorer()


CONTEXT: WorkerContext | None = None


def initialize(tools: Toolchain, design: ExperimentalDesign, tool_threads: str) -> None:
    global CONTEXT
    CONTEXT = WorkerContext(tools, design, tool_threads)


def _context() -> WorkerContext:
    if CONTEXT is None:
        raise RuntimeError("Worker context was never initialized")
    return CONTEXT


def simphy(ils: str, replicate: int) -> tuple[Path, list[Path]]:
    context = _context()
    assert context.tools.simphy is not None
    simulator = SimPhySimulator(context.tools.simphy, context.runner)
    return simulator.run(context.design, ils, replicate)


def verify_ils(
    ils: str, replicate: int, species_tree: Path,
    gene_trees: list[Path], population_size: int,
) -> dict[str, object]:
    return ILSAnalyzer().analyze(
        ils=ils, replicate=replicate, species_tree=species_tree,
        gene_trees=gene_trees, population_size=population_size,
    )


def alisim_locus(ils: str, replicate: int, length: int, locus: int, true_tree: Path) -> Path:
    context = _context()
    assert context.tools.iqtree is not None
    simulator = AliSimSimulator(context.tools.iqtree, context.runner)
    return simulator.run_locus(context.design, ils, replicate, length, locus, true_tree)


def gene_tree_locus(
    condition: ExperimentCondition, locus: int, alignment: Path, true_tree: Path,
) -> tuple[Path, list[BranchObservation]]:
    """Infer one gene tree and score it against its true gene tree.

    Scoring runs in the worker because parsing both Newick trees costs real CPU
    at 300,000 loci and is trivially parallel alongside the inference itself.
    """
    context = _context()
    assert context.tools.iqtree is not None
    inferer = IQTreeInferer(context.tools.iqtree, context.runner)
    tree = inferer.run_locus(
        context.design, condition, locus, alignment, threads=context.tool_threads
    )
    observations = context.scorer.score(
        condition, locus, read_tree(tree), read_tree(true_tree)
    )
    return tree, observations


def species_trees(
    condition: ExperimentCondition, species_tree: Path, estimates: list[Path],
) -> dict[str, object]:
    """Build the four Stage 2 species trees for one condition and score them."""
    context = _context()
    design, tools, runner = context.design, context.tools, context.runner
    assert tools.astral is not None and tools.wastral is not None
    threads = context.tool_threads

    output = design.model_dir(condition) / "species_trees"
    input_trees = output / "estimated_gene_trees.tre"
    combine_newick(estimates, input_trees)
    astral_tree = AstralRunner(tools.astral, runner).run(
        input_trees, output / "astral_unweighted.tre", threads=threads
    )
    threshold_label = f"{design.contracted_abayes_threshold:.2f}"
    contracted_input = output / f"estimated_gene_trees_contracted_abayes_{threshold_label}.tre"
    SupportTreeContractor().contract(
        input_trees, contracted_input, threshold=design.contracted_abayes_threshold,
    )
    contracted_astral_tree = AstralRunner(tools.astral, runner).run(
        contracted_input, output / f"astral_contracted_abayes_{threshold_label}.tre",
        threads=threads,
    )
    weighted_input = output / "estimated_gene_trees_wastral_ready.tre"
    SupportTreePreparer().prepare(input_trees, weighted_input)
    wastral_tree = WeightedAstralRunner(tools.wastral, runner).run(
        weighted_input, output / "wastral_support.tre", threads=threads
    )
    hybrid_wastral_tree = HybridWeightedAstralRunner(tools.wastral, runner).run(
        weighted_input, output / "wastral_hybrid.tre", threads=threads
    )
    return SpeciesTreeEvaluator().evaluate(
        condition=condition, true_species_tree=species_tree,
        astral_tree=astral_tree, contracted_astral_tree=contracted_astral_tree,
        wastral_tree=wastral_tree, hybrid_wastral_tree=hybrid_wastral_tree,
    )
