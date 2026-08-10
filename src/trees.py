"""Newick I/O, bipartition extraction, and SimPhy label handling.

Verified against real files:
  * SimPhy true gene trees      leaves are  <species>_<locus>_<individual>, e.g. 6_0_0
  * SimPhy species trees        leaves are bare integers, e.g. 6
  * real S100 s_tree.trees      bare integers (checked against the file in the
                                official Weighted-ASTRAL_data repository)
  * IQ-TREE --abayes treefile   internal label written with a LEADING SLASH, /0.994
  * avian UCE gene trees        RAxML integer bootstrap labels, missing taxa per locus
"""
from __future__ import annotations

import re

import dendropy

SIMPHY_LEAF = re.compile(r"^(.+?)_\d+_\d+$")


def load_trees(path, tns=None, limit=None):
    """Load a multi-tree Newick file. Internal labels stay on node.label."""
    tl = dendropy.TreeList.get(
        path=str(path), schema="newick",
        suppress_internal_node_taxa=True,
        preserve_underscores=True,
        rooting="force-unrooted",
        taxon_namespace=tns,
    )
    if limit is not None:
        tl = dendropy.TreeList(tl[:limit], taxon_namespace=tl.taxon_namespace)
    for t in tl:
        t.is_rooted = False
        t.deroot()
    return tl


def load_tree(path, tns=None):
    return load_trees(path, tns=tns, limit=1)[0]


def strip_simphy_suffix(tree):
    """Rewrite `6_0_0` -> `6` in place. Returns the number of leaves renamed.

    Only strips when the label actually matches the SimPhy pattern, so calling
    it on a dataset that does not use that convention is a no-op.
    """
    n = 0
    for lf in tree.leaf_node_iter():
        m = SIMPHY_LEAF.match(lf.taxon.label)
        if m:
            lf.taxon.label = m.group(1)
            n += 1
    return n


def clean_support(label):
    """IQ-TREE may write '/0.994', '0.994' or 'a/b/0.994'; RAxML writes '95'.

    Returns the last slash-delimited numeric field, or None.
    """
    if label is None:
        return None
    parts = [p for p in str(label).strip().split("/") if p != ""]
    if not parts:
        return None
    try:
        return float(parts[-1])
    except ValueError:
        return None


def leafset(node):
    return frozenset(lf.taxon.label for lf in node.leaf_iter())


def canonical(side, all_labels):
    """A split is the same read from either side; pick one deterministically."""
    other = all_labels - side
    if len(side) != len(other):
        return side if len(side) < len(other) else other
    return min(side, other, key=lambda s: sorted(s))


def informative_branches(tree):
    """Yield (canonical_split, support) for every internal branch that induces
    a non-trivial split. Root and trivial splits are skipped."""
    all_lbl = frozenset(lf.taxon.label for lf in tree.leaf_node_iter())
    n = len(all_lbl)
    for nd in tree.preorder_node_iter():
        if nd.is_leaf() or nd.parent_node is None:
            continue
        side = leafset(nd)
        if not (1 < len(side) < n - 1):
            continue
        yield canonical(side, all_lbl), clean_support(nd.label)


def bipartitions(tree):
    return {s for s, _ in informative_branches(tree)}


def restrict_split(split, keep):
    """Project a split onto a subset of taxa. Returns None if the projection
    is trivial (fewer than 2 taxa on a side), which happens when a gene tree is
    missing taxa."""
    side = frozenset(split) & keep
    other = keep - side
    if len(side) < 2 or len(other) < 2:
        return None
    return canonical(side, keep)


def write_trees(trees, path):
    trees.write(path=str(path), schema="newick",
                suppress_rooting=True, unquoted_underscores=True,
                suppress_internal_node_labels=False,
                real_value_format_specifier=".10f")
