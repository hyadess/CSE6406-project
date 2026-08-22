from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

import dendropy


Split = frozenset[str]
_SUPPORT = re.compile(r"(?:^|/)([0-9]*\.?[0-9]+)$")
_SIMPHY_SUFFIX = re.compile(r"^(.+)_\d+_\d+$")


def normalize_label(label: str) -> str:
    match = _SIMPHY_SUFFIX.match(label)
    return match.group(1) if match else label


def clean_support(label: object | None) -> float | None:
    if label is None:
        return None
    match = _SUPPORT.search(str(label).strip())
    return float(match.group(1)) if match else None


def read_tree(path: Path):
    return read_tree_text(path.read_text().strip())


def read_tree_text(text: str):
    tree = dendropy.Tree.get(
        data=text, schema="newick", rooting="force-unrooted",
        preserve_underscores=True, suppress_internal_node_taxa=True,
    )
    tree.is_rooted = False
    tree.deroot()
    return tree


def canonical(side: Split, taxa: Split) -> Split:
    other = taxa - side
    if len(side) != len(other):
        return side if len(side) < len(other) else other
    return min(side, other, key=lambda value: tuple(sorted(value)))


def informative_splits(tree, *, with_support: bool = False):
    taxa = frozenset(normalize_label(x.taxon.label) for x in tree.leaf_node_iter())
    for node in tree.preorder_node_iter():
        if node.is_leaf() or node.parent_node is None:
            continue
        side = frozenset(normalize_label(x.taxon.label) for x in node.leaf_iter())
        if 1 < len(side) < len(taxa) - 1:
            split = canonical(side, taxa)
            yield (split, clean_support(node.label)) if with_support else split


def topology_distance(first, second) -> float:
    """Symmetric normalized RF distance after SimPhy label normalization."""
    first_splits = set(informative_splits(first))
    second_splits = set(informative_splits(second))
    denominator = len(first_splits) + len(second_splits)
    if denominator == 0:
        return 0.0
    return len(first_splits.symmetric_difference(second_splits)) / denominator


def combine_newick(paths: Iterable[Path], destination: Path) -> None:
    trees = []
    for path in paths:
        text = path.read_text().strip()
        if not text.endswith(";"):
            raise ValueError(f"Incomplete Newick tree: {path}")
        trees.append(text)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(trees) + "\n")
