from pathlib import Path

from cse6406.trees.newick import clean_support, read_tree_text


class SupportTreeContractor:
    """Contract branches below an aBayes threshold for an ASTRAL baseline."""

    def contract(self, source: Path, destination: Path, *, threshold: float) -> int:
        output, contracted = [], 0
        for text in source.read_text().splitlines():
            if not text.strip():
                continue
            tree = read_tree_text(text)
            for node in list(tree.postorder_node_iter()):
                if node.is_leaf() or node.parent_node is None:
                    continue
                support = clean_support(node.label)
                if support is None or support < threshold:
                    node.edge.collapse(adjust_collapsed_head_children_edge_lengths=False)
                    contracted += 1
            output.append(tree.as_string(
                schema="newick", suppress_rooting=True,
                unquoted_underscores=True,
            ).strip())
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text("\n".join(output) + "\n")
        destination.with_suffix(".contraction.txt").write_text(
            f"contracted_internal_branches\t{contracted}\n"
            f"abayes_threshold\t{threshold:.10g}\n"
            "rule\tcontract support below threshold; missing support is contracted\n"
        )
        return contracted
