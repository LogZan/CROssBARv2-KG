"""Test script for UniprotKeywords adapter using shared implementation."""

from bccb.uniprot_keywords_adapter import (
    UniprotKeywords,
    UniprotKeywordNodeType,
    UniprotKeywordEdgeType,
)


def main() -> None:
    print("Testing UniprotKeywords Adapter")

    node_types = [
        UniprotKeywordNodeType.KEYWORD,
        UniprotKeywordNodeType.KEYWORD_CATEGORY,
        UniprotKeywordNodeType.GO_TERM,
    ]

    edge_types = [
        UniprotKeywordEdgeType.KEYWORD_TO_CATEGORY,
        UniprotKeywordEdgeType.KEYWORD_TO_PARENT,
        UniprotKeywordEdgeType.KEYWORD_TO_GO_TERM,
    ]

    adapter = UniprotKeywords(
        node_types=node_types,
        edge_types=edge_types,
        test_mode=True,
        scan_leaf_paths=False,
    )

    nodes = adapter.get_nodes()
    edges = adapter.get_edges()

    first_node = next(nodes)
    first_edge = next(edges)

    print("Sample node:", first_node)
    print("Sample edge:", first_edge)


if __name__ == "__main__":
    main()
