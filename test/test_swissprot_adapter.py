"""
Test script for UniprotSwissprot adapter with extended SwissProt parsing.
Tests the new node types (annotations, features, diseases) and edge types.
"""
import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from bccb.uniprot_swissprot_adapter import (
    UniprotSwissprot,
    UniprotNodeType,
    UniprotNodeField,
    UniprotEdgeType,
    UniprotIDField,
)

# Test mode flag
TEST_MODE = True

# Configuration - include new extended node and edge types
uniprot_node_types = [
    UniprotNodeType.PROTEIN,
    UniprotNodeType.GENE,
    UniprotNodeType.ORGANISM,
    # Extended types
    UniprotNodeType.FUNCTIONAL_ANNOTATION,
    UniprotNodeType.SEQUENCE_FEATURE,
    UniprotNodeType.UNIPROT_DISEASE,
]

# Use minimal fields for faster testing (skip embeddings)
uniprot_node_fields = [
    UniprotNodeField.PRIMARY_GENE_NAME,
    UniprotNodeField.LENGTH,
    UniprotNodeField.MASS,
    UniprotNodeField.ORGANISM,
    UniprotNodeField.ORGANISM_ID,
    UniprotNodeField.PROTEIN_NAMES,
    UniprotNodeField.PROTEIN_GENE_NAMES,
    UniprotNodeField.ENTREZ_GENE_IDS,
    UniprotNodeField.SEQUENCE,
    # Skip embeddings for faster testing
    # UniprotNodeField.PROTT5_EMBEDDING,
    # UniprotNodeField.ESM2_EMBEDDING,
    # UniprotNodeField.NT_EMBEDDING,
]

uniprot_edge_types = [
    UniprotEdgeType.PROTEIN_TO_ORGANISM,
    UniprotEdgeType.GENE_TO_PROTEIN,
    # Extended edge types
    UniprotEdgeType.PROTEIN_TO_KEYWORD,
    UniprotEdgeType.PROTEIN_TO_ANNOTATION,
    UniprotEdgeType.PROTEIN_TO_FEATURE,
    UniprotEdgeType.PROTEIN_TO_DISEASE,
    UniprotEdgeType.PROTEIN_INTERACTS_PROTEIN,
]

uniprot_id_type = [
    UniprotIDField.GENE_ENTREZ_ID,
]

print("=" * 60)
print("Testing UniprotSwissprot Adapter (Extended Version)")
print("=" * 60)

# Initialize adapter
print("\n1. Initializing adapter...")
adapter = UniprotSwissprot(
    json_path="/GenSIvePFS/users/data/UniProt/UniProtKB_SwissProt/uniprotkb_reviewed_true_2025_11_04.json",
    organism="*",
    node_types=uniprot_node_types,
    node_fields=uniprot_node_fields,
    edge_types=uniprot_edge_types,
    id_fields=uniprot_id_type,
    test_mode=TEST_MODE,
)

print(f"   Data source: {adapter.data_source}")
print(f"   Data version: {adapter.data_version}")
print(f"   Test mode: {adapter.test_mode} (limit: {adapter.test_limit} entries)")

# Load data
print("\n2. Loading data from JSON...")
adapter.download_uniprot_data()

print(f"   Loaded {len(adapter.uniprot_ids)} protein entries")

# Get nodes
print("\n3. Getting nodes (protein, gene, organism, extended)...")
nodes = list(adapter.get_nodes())
basic_node_labels = {"protein", "gene", "organism"}
basic_nodes = [n for n in nodes if n[1] in basic_node_labels]
extended_nodes = [n for n in nodes if n[1] not in basic_node_labels]
print(f"   Generated {len(basic_nodes)} basic nodes")
print(f"   Generated {len(extended_nodes)} extended nodes")

# Show sample nodes
print("\n   Sample basic nodes:")
node_counts = {}
for node_id, node_type, props in nodes[:6]:
    node_counts[node_type] = node_counts.get(node_type, 0) + 1
    if node_counts[node_type] <= 2:
        print(f"   - [{node_type}] {node_id}")
        print(f"     Props: {list(props.keys())[:5]}...")

# Show extended node types
print("\n4. Extended node types found:")
ext_node_types = {}
for node_id, node_type, props in extended_nodes:
    ext_node_types[node_type] = ext_node_types.get(node_type, 0) + 1
for ntype, count in sorted(ext_node_types.items()):
    print(f"   - {ntype}: {count}")

# Get edges
print("\n5. Getting edges...")
edges = list(adapter.get_edges())
basic_edge_labels = {"Protein_belongs_to_organism", "Gene_encodes_protein"}
basic_edges = [e for e in edges if e[3] in basic_edge_labels]
extended_edges = [e for e in edges if e[3] not in basic_edge_labels]
print(f"   Generated {len(basic_edges)} basic edges")
print(f"   Generated {len(extended_edges)} extended edges")

# Show sample basic edges
print("\n   Sample basic edges:")
for edge in basic_edges[:3]:
    _, source, target, label, _ = edge
    print(f"   - {source} --[{label}]--> {target}")

# Show extended edge types and counts
print("\n6. Extended edge types:")
ext_edge_types = {}
for _, source, target, label, props in extended_edges:
    ext_edge_types[label] = ext_edge_types.get(label, 0) + 1
for etype, count in sorted(ext_edge_types.items()):
    print(f"   - {etype}: {count}")

# Show sample extended edges with properties
print("\n   Sample extended edges with properties:")
for edge in extended_edges[:5]:
    _, source, target, label, props = edge
    # Show only non-base properties
    extra_props = {k: v for k, v in props.items() if k not in ('source', 'licence', 'version') and v}
    print(f"   - {source} --[{label}]--> {target}")
    if extra_props:
        # Truncate long text
        for k, v in extra_props.items():
            if isinstance(v, str) and len(v) > 50:
                extra_props[k] = v[:50] + "..."
        print(f"     Props: {extra_props}")

# Summary
annotation_node_types = {ntype for _, ntype, _ in extended_nodes if ntype.endswith("_annotation")}
feature_node_types = {ntype for _, ntype, _ in extended_nodes if ntype.endswith("_feature")}
disease_nodes = [n for n in extended_nodes if n[1] == "uniprot_disease"]
protein_with_keywords = {
    source for _, source, target, label, _ in extended_edges if label == "protein_has_keyword"
}
protein_with_interactions = {
    source for _, source, target, label, _ in extended_edges if label == "Protein_interacts_with_protein"
}

print("\n" + "=" * 60)
print("Test Summary")
print("=" * 60)
print(f"   Protein entries: {len(adapter.uniprot_ids)}")
print(f"   Basic nodes: {len(basic_nodes)}")
print(f"   Extended nodes: {len(extended_nodes)}")
print(f"   Basic edges: {len(basic_edges)}")
print(f"   Extended edges: {len(extended_edges)}")
print(f"   Annotation types: {len(annotation_node_types)}")
print(f"   Feature types: {len(feature_node_types)}")
print(f"   Disease nodes: {len(disease_nodes)}")
print(f"   Proteins with keywords: {len(protein_with_keywords)}")
print(f"   Proteins with interactions: {len(protein_with_interactions)}")
print("=" * 60)
print("Test completed successfully!")
print("=" * 60)
