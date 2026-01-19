import sys
import gc
import yaml
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
# Add parent directory to Python path
sys.path.insert(0, str(project_root))

# CRITICAL: Setup pypath cache directory FIRST before any pypath imports
# This ensures all temp files go to the large filesystem, not /tmp or /root
from bccb import cache_config
cache_config.setup_pypath_cache()

# Patch pypath settings for compatibility
try:
    from pypath.share import settings
    if not hasattr(settings, 'context') and hasattr(settings, 'settings'):
        settings.context = settings.settings.context
except ImportError:
    pass

from bccb.uniprot_swissprot_adapter import (
    UniprotSwissprot,
    UniprotNodeType,
    UniprotNodeField,
    UniprotEdgeType,
    UniprotIDField,
)

from bccb.ppi_adapter import (
    PPI
)

from bccb.interpro_adapter import (
    InterPro
)

from bccb.go_adapter import (
    GO
)

from bccb.drug_adapter import (
    Drug,
    DrugNodeField,
)

from bccb.compound_adapter import (
    Compound
)

from bccb.orthology_adapter import (
    Orthology
)

from bccb.disease_adapter import (
    Disease
)

from bccb.phenotype_adapter import (
    HPO
)

from bccb.pathway_adapter import (
    Pathway
)

from bccb.side_effect_adapter import (
    SideEffect
)

from bccb.ec_adapter import (
    EC
)

from bccb.tfgen_adapter import (
    TFGene
)

from biocypher import BioCypher


def update_schema_with_dynamic_types(schema_path: str, annotation_types: set, feature_types: set):
    """
    Dynamically update schema_config.yaml with discovered annotation and feature types.
    
    Args:
        schema_path: Path to schema_config.yaml
        annotation_types: Set of annotation type names (e.g., 'biophysicochemical_properties_annotation')
        feature_types: Set of feature type names (e.g., 'zinc_finger_feature')
    """
    # Load existing schema
    with open(schema_path, 'r', encoding='utf-8') as f:
        schema = yaml.safe_load(f)
    
    added_types = []
    
    # Add missing annotation types (inherit from functional annotation)
    for ann_type in annotation_types:
        # Convert to schema key format (underscore to space)
        schema_key = ann_type.replace('_', ' ')
        if schema_key not in schema:
            schema[schema_key] = {
                'is_a': 'functional annotation',
                'represented_as': 'node',
                'label_in_input': ann_type,
            }
            added_types.append(schema_key)
    
    # Add missing feature types (inherit from sequence feature)
    for feat_type in feature_types:
        # Convert to schema key format (underscore to space)
        schema_key = feat_type.replace('_', ' ')
        if schema_key not in schema:
            schema[schema_key] = {
                'is_a': 'sequence feature',
                'represented_as': 'node',
                'label_in_input': feat_type,
            }
            added_types.append(schema_key)
    
    if added_types:
        # Write updated schema
        with open(schema_path, 'w', encoding='utf-8') as f:
            yaml.dump(schema, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
        print(f"Schema updated with {len(added_types)} new types:")
        for t in sorted(added_types):
            print(f"  + {t}")
    else:
        print("No new schema types needed.")
    
    return added_types

# Whether to cache data by pypath for future usage
CACHE = True

# Flag for exporting node and edge files as csv format
export_as_csv = True

# Flag for test mode
TEST_MODE = True

# dirs
output_dir_path = "/GenSIvePFS/users/clzeng/workspace/CROssBARv2-KG/biocypher-out"

# user and passwd
drugbank_user = "zengchuanlong23@mails.ucas.ac.cn"
drugbank_passwd = "iHDTbZ3vpFaWzC"

# uniprot configuration
uniprot_node_types = [
    UniprotNodeType.PROTEIN,
    UniprotNodeType.GENE,
    UniprotNodeType.ORGANISM,
    # Extended types from SwissProt
    UniprotNodeType.FUNCTIONAL_ANNOTATION,
    UniprotNodeType.SEQUENCE_FEATURE,
    UniprotNodeType.UNIPROT_DISEASE,
]

uniprot_node_fields = [
    UniprotNodeField.PRIMARY_GENE_NAME,
    UniprotNodeField.LENGTH,
    UniprotNodeField.MASS,
    UniprotNodeField.ORGANISM,
    UniprotNodeField.ORGANISM_ID,
    UniprotNodeField.PROTEIN_NAMES,
    UniprotNodeField.PROTEIN_GENE_NAMES,
    UniprotNodeField.ENSEMBL_TRANSCRIPT_IDS,
    UniprotNodeField.ENSEMBL_GENE_IDS,
    UniprotNodeField.ENTREZ_GENE_IDS,
    UniprotNodeField.KEGG_IDS,
    UniprotNodeField.PROTEOME,
    UniprotNodeField.SEQUENCE,
    # UniprotNodeField.PROTT5_EMBEDDING,  # Disabled: slow to load
    # UniprotNodeField.ESM2_EMBEDDING,     # Disabled: slow to load
    # UniprotNodeField.NT_EMBEDDING,       # Disabled: slow to load
]

uniprot_edge_types = [
     UniprotEdgeType.PROTEIN_TO_ORGANISM,
     UniprotEdgeType.GENE_TO_PROTEIN,
     # Extended edge types from SwissProt
     UniprotEdgeType.PROTEIN_TO_KEYWORD,
     UniprotEdgeType.PROTEIN_TO_ANNOTATION,
     UniprotEdgeType.PROTEIN_TO_FEATURE,
     UniprotEdgeType.PROTEIN_TO_DISEASE,
]

uniprot_id_type = [
     UniprotIDField.GENE_ENTREZ_ID,
]

# ========================================
# Step 1: Load UniProt data FIRST to discover all annotation/feature types
# ========================================
print("=" * 60)
print("Step 1: Loading UniProt data to discover schema types...")
print("=" * 60)

uniprot_adapter = UniprotSwissprot(
        json_path="/GenSIvePFS/users/data/UniProt/UniProtKB_SwissProt/uniprotkb_reviewed_true_2025_11_04.json",
        organism="*",
        node_types=uniprot_node_types,
        node_fields=uniprot_node_fields,
        edge_types=uniprot_edge_types,
        id_fields=uniprot_id_type,
        test_mode=TEST_MODE,
    )

uniprot_adapter.download_uniprot_data(cache=CACHE, retries=6)

# ========================================
# Step 2: Dynamically update schema with discovered types
# ========================================
print("\n" + "=" * 60)
print("Step 2: Updating schema with discovered annotation/feature types...")
print("=" * 60)

schema_path = str(project_root / "config/schema_config.yaml")
annotation_types = set(uniprot_adapter.annotation_nodes.keys())
feature_types = set(uniprot_adapter.feature_nodes.keys())

print(f"Discovered {len(annotation_types)} annotation types and {len(feature_types)} feature types")

update_schema_with_dynamic_types(schema_path, annotation_types, feature_types)

# ========================================
# Step 3: Initialize BioCypher AFTER schema update
# ========================================
print("\n" + "=" * 60)
print("Step 3: Initializing BioCypher with updated schema...")
print("=" * 60)

bc = BioCypher(
    biocypher_config_path=str(project_root / "config/biocypher_config.yaml"),
    schema_config_path=str(project_root / "config/schema_config.yaml"),
)

# ========================================
# Step 4: Write nodes and edges
# ========================================
print("\n" + "=" * 60)
print("Step 4: Writing UniProt nodes and edges...")
print("=" * 60)

uniprot_nodes = uniprot_adapter.get_nodes()
uniprot_edges = uniprot_adapter.get_edges()

# Write basic nodes and edges
bc.write_nodes(uniprot_nodes)
bc.write_edges(uniprot_edges)

# Write extended nodes and edges from SwissProt
uniprot_extended_nodes = uniprot_adapter.get_all_extended_nodes()
uniprot_extended_edges = uniprot_adapter.get_all_extended_edges()

bc.write_nodes(uniprot_extended_nodes)
bc.write_edges(uniprot_extended_edges)

print(f"SwissProt extended nodes and edges written successfully.")
print(f"  Annotation types: {len(uniprot_adapter.annotation_nodes)}")
print(f"  Feature types: {len(uniprot_adapter.feature_nodes)}")
print(f"  Disease nodes: {len(uniprot_adapter.disease_nodes)}")
print(f"  Proteins with keywords: {len(uniprot_adapter.protein_keywords)}")

# PPI
ppi_adapter = PPI(organism=9606, 
                  output_dir=output_dir_path,
                  export_csv=export_as_csv,
                  test_mode=TEST_MODE)

ppi_adapter.download_ppi_data(cache=CACHE)

ppi_adapter.process_ppi_data()

bc.write_edges(ppi_adapter.get_ppi_edges())

# protein domain
# Note: pypath 0.16.28+ uses rescued data from OmniPath, currently only human (9606) available
try:
    interpro_adapter = InterPro(
        organism=9606,  # Use human for now, rescued data only supports 9606
        test_mode=TEST_MODE
    )

    interpro_adapter.download_interpro_data(cache=CACHE)

    if export_as_csv:
        interpro_adapter.export_as_csv(path=output_dir_path)

    bc.write_nodes(interpro_adapter.get_interpro_nodes())
    bc.write_edges(interpro_adapter.get_interpro_edges())
except Exception as e:
    print(f"WARNING: InterPro adapter failed: {e}")
    print("Skipping InterPro data, continuing with other adapters...")
finally:
    del interpro_adapter
    gc.collect()
    print("Memory cleaned up after InterPro adapter")

# gene ontology
try:
    go_adapter = GO(
        organism=9606, 
        test_mode=TEST_MODE
    )
    go_adapter.download_go_data(cache=CACHE)
    bc.write_nodes(go_adapter.get_go_nodes())
    bc.write_edges(go_adapter.get_go_edges())
    if export_as_csv:
        go_adapter.export_as_csv(path=output_dir_path)
except Exception as e:
    print(f"WARNING: GO adapter failed: {e}")
finally:
    del go_adapter
    gc.collect()
    print("Memory cleaned up after GO adapter")

# drug - with memory optimization
try:
    # Exclude SELFORMER_EMBEDDING (requires external file path)
    drug_node_fields = [f for f in DrugNodeField if f != DrugNodeField.SELFORMER_EMBEDDING]
    drug_adapter = Drug(
        drugbank_user=drugbank_user, 
        drugbank_passwd=drugbank_passwd,
        node_fields=drug_node_fields,
        export_csv=export_as_csv, 
        output_dir=output_dir_path,
        test_mode=TEST_MODE,
        low_memory_mode=True  # Enable memory optimization
    )
    drug_adapter.download_drug_data(cache=CACHE)
    drug_adapter.process_drug_data()
    bc.write_nodes(drug_adapter.get_drug_nodes())
    bc.write_edges(drug_adapter.get_edges())
except Exception as e:
    print(f"WARNING: Drug adapter failed: {e}")
    import traceback
    traceback.print_exc()
finally:
    if 'drug_adapter' in dir():
        del drug_adapter
    gc.collect()
    print("Memory cleaned up after Drug adapter")

# compound - with memory optimization
try:
    compound_adapter = Compound(
        stitch_organism=9606,
        export_csv=export_as_csv, 
        output_dir=output_dir_path,
        test_mode=TEST_MODE,
        low_memory_mode=True  # Enable memory optimization
    )
    compound_adapter.download_compound_data(cache=CACHE)
    compound_adapter.process_compound_data()
    bc.write_nodes(compound_adapter.get_compound_nodes())
    bc.write_edges(compound_adapter.get_cti_edges())
except Exception as e:
    print(f"WARNING: Compound adapter failed: {e}")
    import traceback
    traceback.print_exc()
finally:
    if 'compound_adapter' in dir():
        del compound_adapter
    gc.collect()
    print("Memory cleaned up after Compound adapter")

# orthology
try:
    orthology_adapter = Orthology(
        export_csv=export_as_csv, 
        output_dir=output_dir_path,
        test_mode=TEST_MODE
    )
    orthology_adapter.download_orthology_data(cache=CACHE)
    bc.write_edges(orthology_adapter.get_orthology_edges())
except Exception as e:
    print(f"WARNING: Orthology adapter failed: {e}")
finally:
    gc.collect()

# disease
try:
    disease_adapter = Disease(
        drugbank_user=drugbank_user, 
        drugbank_passwd=drugbank_passwd,
        export_csv=export_as_csv, 
        output_dir=output_dir_path,
        test_mode=TEST_MODE
    )
    disease_adapter.download_disease_data(cache=CACHE)
    bc.write_nodes(disease_adapter.get_nodes())
    bc.write_edges(disease_adapter.get_edges())
except Exception as e:
    print(f"WARNING: Disease adapter failed: {e}")
finally:
    gc.collect()

# phenotype
try:
    phenotype_adapter = HPO(
        export_csv=export_as_csv, 
        output_dir=output_dir_path,
        test_mode=TEST_MODE
    )
    phenotype_adapter.download_hpo_data(cache=CACHE)
    bc.write_nodes(phenotype_adapter.get_nodes())
    bc.write_edges(phenotype_adapter.get_edges())
except Exception as e:
    print(f"WARNING: Phenotype adapter failed: {e}")
finally:
    gc.collect()

# pathway
try:
    pathway_adapter = Pathway(
        drugbank_user=drugbank_user, 
        drugbank_passwd=drugbank_passwd,
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        test_mode=TEST_MODE,
    )
    pathway_adapter.download_pathway_data(cache=CACHE)
    bc.write_nodes(pathway_adapter.get_nodes())
    bc.write_edges(pathway_adapter.get_edges())
except Exception as e:
    print(f"WARNING: Pathway adapter failed: {e}")
    import traceback
    traceback.print_exc()
finally:
    if 'pathway_adapter' in dir():
        del pathway_adapter
    gc.collect()
    print("Memory cleaned up after Pathway adapter")


# side effect
try:
    side_effect_adapter = SideEffect(
        drugbank_user=drugbank_user,
        drugbank_passwd=drugbank_passwd,
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        test_mode=TEST_MODE
    )
    side_effect_adapter.download_side_effect_data(cache=CACHE)
    bc.write_nodes(side_effect_adapter.get_nodes())
    bc.write_edges(side_effect_adapter.get_edges())
except Exception as e:
    print(f"WARNING: Side effect adapter failed: {e}")
    import traceback
    traceback.print_exc()
finally:
    if 'side_effect_adapter' in dir():
        del side_effect_adapter
    gc.collect()
    print("Memory cleaned up after Side effect adapter")

# ec numbers
try:
    ec_adapter = EC(
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        test_mode=TEST_MODE
    )
    ec_adapter.download_ec_data(cache=CACHE)
    bc.write_nodes(ec_adapter.get_nodes())
    bc.write_edges(ec_adapter.get_edges())
except Exception as e:
    print(f"WARNING: EC adapter failed: {e}")
    import traceback
    traceback.print_exc()
finally:
    if 'ec_adapter' in dir():
        del ec_adapter
    gc.collect()
    print("Memory cleaned up after EC adapter")

# tf-gen
try:
    tfgene_adapter = TFGene(
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        test_mode=TEST_MODE
    )
    tfgene_adapter.download_tfgen_data(cache=CACHE)
    bc.write_edges(tfgene_adapter.get_edges())
except Exception as e:
    print(f"WARNING: TFGene adapter failed: {e}")
    import traceback
    traceback.print_exc()
finally:
    if 'tfgene_adapter' in dir():
        del tfgene_adapter
    gc.collect()
    print("Memory cleaned up after TFGene adapter")


# Write import call and other post-processing
bc.write_import_call()
bc.summary()

print("="*80)
print("Script completed successfully!")
print("="*80)
