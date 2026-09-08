from pathlib import Path

from cse6406.trees.newick import clean_support, read_tree_text


LOCAL_BAYESIAN_MINIMUM = 1 / 3
IQTREE_THREE_DECIMAL_TOLERANCE = 0.0005


class SupportTreePreparer:
    """Create valid wASTRAL input and count conservative support imputations."""

    def prepare(
        self, source: Path, destination: Path, *, missing_support: float = 1 / 3,
    ) -> int:
        output, imputed, minimum_adjustments = [], 0, 0
        for text in source.read_text().splitlines():
            if not text.strip():
                continue
            tree = read_tree_text(text)
            for node in tree.preorder_node_iter():
                if node.is_leaf() or node.parent_node is None:
                    continue
                support = clean_support(node.label)
                if support is None:
                    support = missing_support
                    imputed += 1
                elif (
                    support < LOCAL_BAYESIAN_MINIMUM
                    and LOCAL_BAYESIAN_MINIMUM - support
                    <= IQTREE_THREE_DECIMAL_TOLERANCE
                ):
                    # IQ-TREE prints aBayes labels to three decimals, so its
                    # theoretical minimum of 1/3 can appear as 0.333.
                    support = LOCAL_BAYESIAN_MINIMUM
                    minimum_adjustments += 1
                node.label = f"{support:.17g}"
            output.append(tree.as_string(
                schema="newick", suppress_rooting=True,
                unquoted_underscores=True,
            ).strip())
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text("\n".join(output) + "\n")
        destination.with_suffix(".support_imputation.txt").write_text(
            f"imputed_internal_branch_labels\t{imputed}\n"
            f"imputed_value\t{missing_support:.17g}\n"
            f"rounded_minimum_adjustments\t{minimum_adjustments}\n"
            "reason\twASTRAL mode 2 requires numeric non-root internal labels; "
            "1/3 is the documented local-Bayesian minimum and receives minimum weight; "
            "IQ-TREE's rounded 0.333 is normalized to exactly 1/3.\n"
        )
        return imputed
