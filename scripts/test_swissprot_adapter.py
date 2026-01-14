"""
Test script for UniprotSwissprot adapter.
Only runs the swissprot adapter in test mode.
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

from biocypher import BioCypher

# Output directory
output_dir = project_root / "biocypher-out"
output_dir.mkdir(exist_ok=True)

bc = BioCypher(
    biocypher_config_path=str(project_root / "config/biocypher_config.yaml"),
    schema_config_path=str(project_root / "config/schema_config.yaml"),
)

# Test mode flag
TEST_MODE = True

# Configuration - only include fields that don't require external embeddings for quick testing
uniprot_node_types = [
    UniprotNodeType.PROTEIN,
    UniprotNodeType.GENE,
    UniprotNodeType.ORGANISM,
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
    # Uncomment to include embeddings (requires h5 files)
    # UniprotNodeField.PROTT5_EMBEDDING,
    # UniprotNodeField.ESM2_EMBEDDING,
    # UniprotNodeField.NT_EMBEDDING,
]

uniprot_edge_types = [
    UniprotEdgeType.PROTEIN_TO_ORGANISM,
    UniprotEdgeType.GENE_TO_PROTEIN,
]

uniprot_id_type = [
    UniprotIDField.GENE_ENTREZ_ID,
]

print("=" * 60)
print("Testing UniprotSwissprot Adapter")
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
print("\n3. Getting nodes...")
nodes = list(adapter.get_nodes())
print(f"   Generated {len(nodes)} nodes")

# Show sample nodes
print("\n   Sample nodes:")
node_counts = {}
for node_id, node_type, props in nodes[:10]:
    node_counts[node_type] = node_counts.get(node_type, 0) + 1
    if node_counts[node_type] <= 2:
        print(f"   - [{node_type}] {node_id}")
        print(f"     Props: {list(props.keys())[:5]}...")

# Get edges
print("\n4. Getting edges...")
edges = adapter.get_edges()
if edges:
    edges = list(edges)
    print(f"   Generated {len(edges)} edges")
    
    # Show sample edges
    print("\n   Sample edges:")
    for edge in edges[:5]:
        _, source, target, label, _ = edge
        print(f"   - {source} --[{label}]--> {target}")
else:
    print("   No edges generated.")

# Write to BioCypher
print("\n5. Writing to BioCypher...")
bc.write_nodes(adapter.get_nodes())
if edges:
    bc.write_edges(adapter.get_edges())
bc.write_import_call()

print("\n" + "=" * 60)
print("Test completed!")
print(f"Output directory: {output_dir}")
print("=" * 60)
