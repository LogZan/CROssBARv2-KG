"""
Test script for UniprotSwissprot adapter to compare with expected output from original adapter.
Since the original adapter requires pypath (which has dependency issues), we compare against 
the expected output structure based on code analysis.
"""
import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

print("=" * 70)
print("UniprotSwissprot Adapter Output Comparison Test")
print("=" * 70)

# Expected properties from original Uniprot adapter (based on code analysis)
EXPECTED_PROTEIN_PROPS = {
    "length",
    "mass",
    "organism_name",
    "organism_id",
    "protein_names",          # mapped from protein_name
    "primary_protein_name",   # first element of protein_names
    "ec",
    "xref_proteomes",
    "sequence",
    "subcellular_location",
    "source",
    "licence",
    "version",
    # Embeddings (optional)
    # "prott5_embedding",
    # "esm2_embedding",
}

EXPECTED_GENE_PROPS = {
    "gene_symbol",            # mapped from gene_primary
    "gene_names",             # original field
    "xref_geneid",            # entrez gene ids (kept as is or renamed)
    "ensembl_transcript_ids", # mapped from xref_ensembl
    "ensembl_gene_ids",       # derived from ENST
    "kegg_ids",               # mapped from xref_kegg
    "source",
    "licence",
    "version",
    # "nt_embedding",
}

EXPECTED_ORGANISM_PROPS = {
    "organism_name",
    "source",
    "licence",
    "version",
}

# ==============================================================================
# Test SwissProt Adapter
# ==============================================================================
from bccb.uniprot_swissprot_adapter import (
    UniprotSwissprot,
    UniprotNodeType,
    UniprotNodeField,
    UniprotEdgeType,
    UniprotIDField,
)

# All non-embedding fields
node_fields = [
    UniprotNodeField.PRIMARY_GENE_NAME,
    UniprotNodeField.LENGTH,
    UniprotNodeField.MASS,
    UniprotNodeField.ORGANISM,
    UniprotNodeField.ORGANISM_ID,
    UniprotNodeField.PROTEIN_NAMES,
    UniprotNodeField.PROTEIN_GENE_NAMES,
    UniprotNodeField.ENTREZ_GENE_IDS,
    UniprotNodeField.ENSEMBL_TRANSCRIPT_IDS,
    UniprotNodeField.ENSEMBL_GENE_IDS,
    UniprotNodeField.KEGG_IDS,
    UniprotNodeField.PROTEOME,
    UniprotNodeField.SEQUENCE,
    UniprotNodeField.EC,
    UniprotNodeField.SUBCELLULAR_LOCATION,
]

print("\n1. Initializing UniprotSwissprot adapter...")
adapter = UniprotSwissprot(
    organism=9606,  # Human only
    node_types=[UniprotNodeType.PROTEIN, UniprotNodeType.GENE, UniprotNodeType.ORGANISM],
    node_fields=node_fields,
    edge_types=[UniprotEdgeType.PROTEIN_TO_ORGANISM, UniprotEdgeType.GENE_TO_PROTEIN],
    id_fields=[UniprotIDField.GENE_ENTREZ_ID],
    test_mode=True,
)

print("2. Loading data from local JSON...")
adapter.download_uniprot_data()
print(f"   Loaded {len(adapter.uniprot_ids)} proteins")

print("\n3. Analyzing output structure...")
nodes = list(adapter.get_nodes())
edges = adapter.get_edges()
if edges:
    edges = list(edges)
else:
    edges = []

# Collect actual properties
actual_protein_props = set()
actual_gene_props = set()
actual_organism_props = set()
node_counts = {}

sample_protein = None
sample_gene = None
sample_organism = None

for node_id, node_type, props in nodes:
    node_counts[node_type] = node_counts.get(node_type, 0) + 1
    if node_type == "protein":
        actual_protein_props.update(props.keys())
        if sample_protein is None:
            sample_protein = (node_id, props)
    elif node_type == "gene":
        actual_gene_props.update(props.keys())
        if sample_gene is None:
            sample_gene = (node_id, props)
    elif node_type == "organism":
        actual_organism_props.update(props.keys())
        if sample_organism is None:
            sample_organism = (node_id, props)

print(f"   Generated {len(nodes)} nodes, {len(edges)} edges")
print(f"   Node distribution: {node_counts}")

# ==============================================================================
# Comparison
# ==============================================================================
print("\n" + "=" * 70)
print("COMPARISON RESULTS")
print("=" * 70)

def compare_props(expected, actual, category):
    print(f"\n[{category}]")
    common = expected & actual
    missing = expected - actual
    extra = actual - expected
    
    print(f"  ✅ Common ({len(common)}): {sorted(common)}")
    if missing:
        print(f"  ❌ MISSING ({len(missing)}): {sorted(missing)}")
    else:
        print(f"  ✅ No missing properties")
    if extra:
        print(f"  ➕ Extra ({len(extra)}): {sorted(extra)}")
    
    return len(missing) == 0

protein_ok = compare_props(EXPECTED_PROTEIN_PROPS, actual_protein_props, "PROTEIN PROPERTIES")
gene_ok = compare_props(EXPECTED_GENE_PROPS, actual_gene_props, "GENE PROPERTIES")
organism_ok = compare_props(EXPECTED_ORGANISM_PROPS, actual_organism_props, "ORGANISM PROPERTIES")

# ==============================================================================
# Sample Data
# ==============================================================================
print("\n" + "=" * 70)
print("SAMPLE DATA (for verification)")
print("=" * 70)

if sample_protein:
    print(f"\n[Sample Protein] {sample_protein[0]}")
    for k, v in sample_protein[1].items():
        val_str = str(v)[:80] + "..." if len(str(v)) > 80 else str(v)
        print(f"  {k}: {val_str}")

if sample_gene:
    print(f"\n[Sample Gene] {sample_gene[0]}")
    for k, v in sample_gene[1].items():
        val_str = str(v)[:80] + "..." if len(str(v)) > 80 else str(v)
        print(f"  {k}: {val_str}")

if sample_organism:
    print(f"\n[Sample Organism] {sample_organism[0]}")
    for k, v in sample_organism[1].items():
        val_str = str(v)[:80] + "..." if len(str(v)) > 80 else str(v)
        print(f"  {k}: {val_str}")

# Edge samples
print("\n[Sample Edges]")
for edge in edges[:5]:
    _, source, target, label, _ = edge
    print(f"  {source} --[{label}]--> {target}")

print("\n" + "=" * 70)
if protein_ok and gene_ok and organism_ok:
    print("✅ ALL CHECKS PASSED!")
else:
    print("❌ SOME PROPERTIES ARE MISSING - NEEDS FIX")
print("=" * 70)
