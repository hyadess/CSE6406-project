"""Unit tests for trees.py, using the label formats actually seen in the wild."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import dendropy
from trees import (clean_support, canonical, informative_branches, bipartitions,
                   restrict_split, strip_simphy_suffix)

ALL = frozenset("ABCDE")


def parse(nwk):
    t = dendropy.Tree.get(data=nwk, schema="newick",
                          suppress_internal_node_taxa=True,
                          preserve_underscores=True, rooting="force-unrooted")
    t.deroot()
    return t


def test_clean_support():
    cases = [("/0.994", 0.994),      # IQ-TREE --abayes, leading slash
             ("0.994", 0.994),       # plain
             ("95", 95.0),           # RAxML bootstrap
             ("88/0.97", 0.97),      # multiple support types
             ("", None), (None, None), ("NA", None)]
    for raw, want in cases:
        got = clean_support(raw)
        assert got == want, f"clean_support({raw!r})={got}, want {want}"
    print("ok  support labels: IQ-TREE slash, RAxML integer, multi-field")


def test_informative_branches_skips_trivial():
    t = parse("(A:0.1,B:0.1,((C:0.1,D:0.1)/0.9:0.1,E:0.1)/0.7:0.1);")
    got = {tuple(sorted(s)): v for s, v in informative_branches(t)}
    assert (("C", "D"), 0.9) in [(k, v) for k, v in got.items()], got
    for side in got:
        assert 1 < len(side) < len(ALL) - 1, f"trivial split kept: {side}"
    print(f"ok  informative branches only: {sorted(got)}")


def test_canonical_is_side_independent():
    assert canonical(frozenset("AB"), ALL) == canonical(frozenset("CDE"), ALL)
    print("ok  splits canonical regardless of which side is named")


def test_restrict_split_handles_missing_taxa():
    keep = frozenset("ABCD")
    assert restrict_split(frozenset("AB"), keep) == canonical(frozenset("AB"), keep)
    # projecting onto a set where one side has <2 taxa must be dropped, not
    # silently scored -- that is how missing data corrupts a comparison
    assert restrict_split(frozenset("AE"), keep) is None
    print("ok  splits that become trivial under missing data are dropped")


def test_strip_simphy_suffix():
    t = parse("(6_0_0:0.1,8_0_0:0.1,(10_0_0:0.1,13_0_0:0.1)/1:0.1);")
    n = strip_simphy_suffix(t)
    labels = sorted(lf.taxon.label for lf in t.leaf_node_iter())
    assert n == 4 and labels == ["10", "13", "6", "8"], (n, labels)
    t2 = parse("(TAEGU:0.1,GEOFO:0.1,(CORBR:0.1,MANVI:0.1)95:0.1);")
    assert strip_simphy_suffix(t2) == 0, "must be a no-op on non-SimPhy labels"
    print("ok  SimPhy 6_0_0 -> 6, and a no-op on avian-style labels")


def test_bipartitions_match_between_relabelled_trees():
    a = parse("(6_0_0:0.1,8_0_0:0.1,(10_0_0:0.1,13_0_0:0.1)/1:0.1);")
    b = parse("(6:0.2,8:0.2,(10:0.2,13:0.2)95:0.2);")
    strip_simphy_suffix(a)
    assert bipartitions(a) == bipartitions(b)
    print("ok  identical topologies match after label reconciliation")


if __name__ == "__main__":
    for fn in (test_clean_support, test_informative_branches_skips_trivial,
               test_canonical_is_side_independent,
               test_restrict_split_handles_missing_taxa,
               test_strip_simphy_suffix,
               test_bipartitions_match_between_relabelled_trees):
        fn()
    print("\nall trees.py unit tests passed")
