from pathlib import Path

from cse6406.trees.newick import clean_support, read_tree_text


class SupportTreePreparer:
    """Create valid wASTRAL input and count conservative support imputations."""

    def prepare(
        self, source: Path, destination: Path, *, missing_support: float = 1 / 3,
    ) -> int:
        output, imputed = [], 0
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
                node.label = f"{support:.10g}"
            output.append(tree.as_string(
                schema="newick", suppress_rooting=True,
                unquoted_underscores=True,
            ).strip())
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text("\n".join(output) + "\n")
        destination.with_suffix(".support_imputation.txt").write_text(
            f"imputed_internal_branch_labels\t{imputed}\n"
            f"imputed_value\t{missing_support:.10g}\n"
            "reason\twASTRAL mode 2 requires numeric non-root internal labels; "
            "1/3 is the documented local-Bayesian minimum and receives minimum weight.\n"
        )
        return imputed
