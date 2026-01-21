import os
import sys
import gc
import yaml
import ctypes
from datetime import datetime, timezone, timedelta
from pathlib import Path

# Setup timezone (UTC+8 for Asia/Shanghai)
TZ = timezone(timedelta(hours=8))

project_root = Path(__file__).resolve().parent.parent
# Add parent directory to Python path
sys.path.insert(0, str(project_root))

# CRITICAL: Setup pypath cache directory FIRST before any pypath imports
# This ensures all temp files go to the large filesystem, not /tmp or /root
from bccb import cache_config
cache_config.setup_pypath_cache()


def aggressive_memory_cleanup(adapter_name: str = ""):
    """Perform aggressive memory cleanup between adapters."""
    # Force garbage collection multiple times
    gc.collect()
    gc.collect()
    gc.collect()
    
    # Try to release memory back to OS (Linux specific)
    try:
        libc = ctypes.CDLL("libc.so.6")
        libc.malloc_trim(0)
    except Exception:
        pass
    
    # Print memory status
    try:
        import resource
        mem_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024  # MB
        print(f"[{adapter_name}] Memory cleanup done. Peak RSS: {mem_usage:.0f} MB")
    except Exception:
        print(f"[{adapter_name}] Memory cleanup done.")


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

from bccb.uniprot_keywords_adapter import (
    UniprotKeywords,
    KeywordNodeType,
    KeywordEdgeType,
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

# Load configuration
with open(project_root / "config/crossbar_config.yaml", 'r') as f:
    config = yaml.safe_load(f)

# Extract config values
timestamp = datetime.now(TZ).strftime("%Y%m%d%H%M%S")
output_dir_path = str(project_root / "biocypher-out")
embeddings_dir = config['data_paths']['embeddings_dir']
malacards_dir_path = config['data_paths']['malacards_dir']
uniprot_json_path = config['data_paths']['uniprot_json']

# Embedding file paths
prott5_embedding_path = f"{embeddings_dir}/{config['embeddings']['prott5']}"
esm2_embedding_path = f"{embeddings_dir}/{config['embeddings']['esm2']}"
nt_embedding_path = f"{embeddings_dir}/{config['embeddings']['nucleotide_transformer']}"
selformer_drug_embedding_path = f"{embeddings_dir}/{config['embeddings']['selformer_drug']}"
selformer_compound_embedding_path = f"{embeddings_dir}/{config['embeddings']['selformer_compound']}"
doc2vec_disease_embedding_path = f"{embeddings_dir}/{config['embeddings']['doc2vec_disease']}"
biokeen_pathway_embedding_path = f"{embeddings_dir}/{config['embeddings']['biokeen_pathway']}"
rxnfp_ec_embedding_path = f"{embeddings_dir}/{config['embeddings']['rxnfp_ec']}"
anc2vec_go_embedding_path = f"{embeddings_dir}/{config['embeddings']['anc2vec_go']}"
cada_phenotype_embedding_path = f"{embeddings_dir}/{config['embeddings']['cada_phenotype']}"
dom2vec_domain_embedding_path = f"{embeddings_dir}/{config['embeddings']['dom2vec_domain']}"

# MalaCards file paths
malacards_json_path = f"{malacards_dir_path}/{config['malacards']['diseases']}"
malacards_related_diseases_json_path = f"{malacards_dir_path}/{config['malacards']['related_diseases']}"


# Organism parameter conversion helpers for pypath compatibility
# TODO: Unused, now setting organism=None in majority of uniprot.uniprot_data() calls
def organism_for_uniprot_data(organism):
    """Convert organism value for uniprot.uniprot_data() calls.

    In pypath v0.16+, uniprot_data() uses None for all organisms, not '*'.
    """
    return None if organism == "*" else organism

def organism_for_all_uniprots(organism):
    """Convert organism value for uniprot._all_uniprots() calls.

    _all_uniprots() uses '*' for all organisms.
    """
    return organism


# Helper for consistent logging
def log_adapter_boundary(adapter_name: str, phase: str):
    """Log adapter execution boundary."""
    separator = "=" * 60
    timestamp = datetime.now(TZ).strftime("%Y-%m-%d %H:%M:%S")
    if phase == "start":
        print(f"\n{separator}")
        print(f"[{timestamp}] [START] {adapter_name} adapter")
        print(f"{separator}")
    else:
        print(f"[{timestamp}] [END] {adapter_name} adapter")
        print(f"{separator}\n")

bc = BioCypher(biocypher_config_path= str(project_root / "config/biocypher_config.yaml"),
               schema_config_path= str(project_root / "config/schema_config.yaml")
)

# Load settings from config
CACHE = config['settings']['cache']
export_as_csv = config['settings']['export_csv']
TEST_MODE = config['settings']['test_mode']
UPDATE_SCHEMA_DYNAMICALLY = config['settings']['update_schema_dynamically']
ORGANISM = config['settings']['organism']
USE_EMBEDDINGS = config['settings']['use_embeddings']


def update_schema_with_dynamic_types(schema_path: str, annotation_types: set, feature_types: set):
    """
    Dynamically update schema_config.yaml with discovered annotation and feature types.
    Only needed on first run or when new annotation/feature types are discovered.
    
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
        # Write updated schema with blank lines between entries
        yaml_str = yaml.dump(schema, default_flow_style=False, allow_unicode=True, sort_keys=False)
        lines = yaml_str.split('\n')
        formatted_lines = []
        for i, line in enumerate(lines):
            if i > 0 and line and not line[0].isspace() and formatted_lines and formatted_lines[-1]:
                formatted_lines.append('')
            formatted_lines.append(line)

        with open(schema_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(formatted_lines))

        print(f"Schema updated with {len(added_types)} new types:")
        for t in sorted(added_types):
            print(f"  + {t}")
    else:
        print("No new schema types needed.")
    
    return added_types


# Load credentials and adapter configurations from config
drugbank_user = config['drugbank']['user']
drugbank_passwd = config['drugbank']['password']

# UniProt configuration
uniprot_node_types = [getattr(UniprotNodeType, nt) for nt in config['uniprot']['node_types']]
uniprot_node_fields = [getattr(UniprotNodeField, nf) for nf in config['uniprot']['node_fields']]
uniprot_edge_types = [getattr(UniprotEdgeType, et) for et in config['uniprot']['edge_types']]
uniprot_id_type = [getattr(UniprotIDField, idf) for idf in config['uniprot']['id_fields']]

log_adapter_boundary("UniProt SwissProt", "start")
try:
    uniprot_adapter = UniprotSwissprot(
            json_path=uniprot_json_path,
            organism=ORGANISM,
            node_types=uniprot_node_types,
            node_fields=uniprot_node_fields,
            edge_types=uniprot_edge_types,
            id_fields=uniprot_id_type,
            test_mode=TEST_MODE,
        )

    uniprot_adapter.download_uniprot_data(cache=CACHE,
            prott5_embedding_output_path=prott5_embedding_path if USE_EMBEDDINGS else None,
            esm2_embedding_path=esm2_embedding_path if USE_EMBEDDINGS else None,
            nucleotide_transformer_embedding_path=nt_embedding_path if USE_EMBEDDINGS else None)

    # Optionally update schema with dynamically discovered types
    if UPDATE_SCHEMA_DYNAMICALLY:
        print("Updating schema with dynamically discovered annotation/feature types...")
        schema_path = str(project_root / "config/schema_config.yaml")
        annotation_types = set(uniprot_adapter.annotation_nodes.keys())
        feature_types = set(uniprot_adapter.feature_nodes.keys())
        update_schema_with_dynamic_types(schema_path, annotation_types, feature_types)
        print("Schema updated. Please restart the script for changes to take effect.")
        sys.exit(0)

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
    if export_as_csv:
        uniprot_adapter.export_data_to_csv(path=output_dir_path,
                                            node_data=uniprot_nodes,
                                            edge_data=uniprot_edges)
    print(f"SwissProt data exported to CSV successfully.")
except Exception as e:
    print(f"WARNING: UniProt SwissProt adapter failed: {e}")
    import traceback
    traceback.print_exc()
finally:
    log_adapter_boundary("UniProt SwissProt", "end")
print(f"  Annotation types: {len(uniprot_adapter.annotation_nodes)}")
print(f"  Feature types: {len(uniprot_adapter.feature_nodes)}")
print(f"  Disease nodes: {len(uniprot_adapter.disease_nodes)}")
print(f"  Proteins with keywords: {len(uniprot_adapter.protein_keywords)}")

# UniProt Keywords (vocabulary with hierarchy and GO mappings)
log_adapter_boundary("UniProt Keywords", "start")
try:
    print("\n" + "="*60)
    print("Loading UniProt Keywords...")
    print("="*60)

    keywords_adapter = UniprotKeywords(
        json_path=config['data_paths']['uniprot_keywords_json'],
        node_types=[getattr(KeywordNodeType, nt) for nt in config['keywords']['node_types']],
        edge_types=[getattr(KeywordEdgeType, et) for et in config['keywords']['edge_types']],
        test_mode=TEST_MODE,
    )
    
    keywords_adapter.download_data(cache=CACHE)
    
    # Write keyword nodes
    keyword_nodes = keywords_adapter.get_nodes()
    bc.write_nodes(keyword_nodes)
    
    # Write keyword edges (hierarchy and GO mappings)
    keyword_edges = keywords_adapter.get_edges()
    bc.write_edges(keyword_edges)
    
    print(f"UniProt Keywords written successfully.")
    print(f"  Keywords loaded: {len(keywords_adapter.keywords_data)}")
    
except Exception as e:
    print(f"WARNING: UniProt Keywords adapter failed: {e}")
    import traceback
    traceback.print_exc()
finally:
    if 'keywords_adapter' in dir():
        del keywords_adapter
    gc.collect()
    aggressive_memory_cleanup("UniProt Keywords")
    log_adapter_boundary("UniProt Keywords", "end")

# PPI
log_adapter_boundary("PPI", "start")
try:
    ppi_adapter = PPI(organism=ORGANISM, 
                      output_dir=output_dir_path,
                      export_csv=export_as_csv,
                      test_mode=TEST_MODE)
    
    ppi_adapter.download_ppi_data(cache=CACHE)
    
    ppi_adapter.process_ppi_data()
    
    bc.write_edges(ppi_adapter.get_ppi_edges())
except Exception as e:
    print(f"WARNING: PPI adapter failed: {e}")
    import traceback
    traceback.print_exc()
finally:
    if 'ppi_adapter' in dir():
        del ppi_adapter
    aggressive_memory_cleanup("PPI")
    log_adapter_boundary("PPI", "end")

# protein domain
# Note: pypath 0.16.28+ uses rescued data from OmniPath, currently only human (9606) available
log_adapter_boundary("InterPro", "start")
try:
    interpro_adapter = InterPro(
        organism=9606,  # Use human for now, rescued data only supports 9606
        test_mode=TEST_MODE
    )

    interpro_adapter.download_interpro_data(cache=CACHE)
    interpro_adapter.download_domain_node_data(dom2vec_embedding_path=dom2vec_domain_embedding_path if USE_EMBEDDINGS else None, cache=CACHE)

    bc.write_nodes(interpro_adapter.get_interpro_nodes())
    bc.write_edges(interpro_adapter.get_interpro_edges())

    if export_as_csv:
        interpro_adapter.export_as_csv(path=output_dir_path)
except Exception as e:
    print(f"WARNING: InterPro adapter failed: {e}")
    print("Skipping InterPro data, continuing with other adapters...")
finally:
    if 'interpro_adapter' in dir():
        del interpro_adapter
    aggressive_memory_cleanup("InterPro")
    log_adapter_boundary("InterPro", "end")

# gene ontology
log_adapter_boundary("GO", "start")
try:
    go_adapter = GO(
        organism=ORGANISM,
        test_mode=TEST_MODE
    )
    go_adapter.download_go_data(cache=CACHE, anc2vec_embedding_path=anc2vec_go_embedding_path if USE_EMBEDDINGS else None)
    bc.write_nodes(go_adapter.get_go_nodes())
    bc.write_edges(go_adapter.get_go_edges())
    if export_as_csv:
        go_adapter.export_as_csv(path=output_dir_path)
except Exception as e:
    print(f"WARNING: GO adapter failed: {e}")
finally:
    if 'go_adapter' in dir():
        del go_adapter
    aggressive_memory_cleanup("GO")
    log_adapter_boundary("GO", "end")

# drug - with memory optimization
log_adapter_boundary("Drug", "start")
try:
    drug_adapter = Drug(
        drugbank_user=drugbank_user,
        drugbank_passwd=drugbank_passwd,
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        test_mode=TEST_MODE,
    )
    drug_adapter.download_drug_data(cache=CACHE, selformer_embedding_path=selformer_drug_embedding_path if USE_EMBEDDINGS else None)
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
    aggressive_memory_cleanup("Drug")
    log_adapter_boundary("Drug", "end")



log_adapter_boundary("Compound", "start")
try:
    compound_adapter = Compound(
        stitch_organism=ORGANISM,
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        test_mode=TEST_MODE,
    )
    compound_adapter.download_compound_data(cache=CACHE, selformer_embedding_path=selformer_compound_embedding_path if USE_EMBEDDINGS else None)
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
    aggressive_memory_cleanup("Compound")
    log_adapter_boundary("Compound", "end")



log_adapter_boundary("Orthology", "start")
try:
    orthology_adapter = Orthology(
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        test_mode=TEST_MODE
    )
    orthology_adapter.download_orthology_data(cache=CACHE)
    # Note: Orthology only has edges, no nodes
    bc.write_edges(orthology_adapter.get_orthology_edges())
except Exception as e:
    print(f"WARNING: Orthology adapter failed: {e}")
finally:
    if 'orthology_adapter' in dir():
        del orthology_adapter
    aggressive_memory_cleanup("Orthology")
    log_adapter_boundary("Orthology", "end")

# disease
log_adapter_boundary("Disease", "start")
try:
    disease_adapter = Disease(
        drugbank_user=drugbank_user,
        drugbank_passwd=drugbank_passwd,
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        test_mode=TEST_MODE
    )
    disease_adapter.download_disease_data(cache=CACHE,
        doc2vec_embedding_path=doc2vec_disease_embedding_path if USE_EMBEDDINGS else None,
        malacards_json_path=malacards_json_path,
        malacards_related_diseases_json_path=malacards_related_diseases_json_path
        )
    bc.write_nodes(disease_adapter.get_nodes())
    bc.write_edges(disease_adapter.get_edges())
except Exception as e:
    print(f"WARNING: Disease adapter failed: {e}")
finally:
    if 'disease_adapter' in dir():
        del disease_adapter
    aggressive_memory_cleanup("Disease")
    log_adapter_boundary("Disease", "end")

# phenotype
log_adapter_boundary("Phenotype", "start")
try:
    phenotype_adapter = HPO(
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        test_mode=TEST_MODE
    )
    phenotype_adapter.download_hpo_data(cache=CACHE, cada_embedding_path=cada_phenotype_embedding_path if USE_EMBEDDINGS else None)
    bc.write_nodes(phenotype_adapter.get_nodes())
    bc.write_edges(phenotype_adapter.get_edges())
except Exception as e:
    print(f"WARNING: Phenotype adapter failed: {e}")
finally:
    if 'phenotype_adapter' in dir():
        del phenotype_adapter
    aggressive_memory_cleanup("Phenotype")
    log_adapter_boundary("Phenotype", "end")

# pathway
log_adapter_boundary("Pathway", "start")
try:
    # Use config kegg_organism, or limit to human in test mode
    kegg_organism = config['pathway']['kegg_organism'] if not TEST_MODE else ["hsa"]
    pathway_adapter = Pathway(
        drugbank_user=drugbank_user,
        drugbank_passwd=drugbank_passwd,
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        test_mode=TEST_MODE,
        kegg_organism=kegg_organism,
    )
    pathway_adapter.download_pathway_data(cache=CACHE, biokeen_embedding_path=biokeen_pathway_embedding_path if USE_EMBEDDINGS else None)
    bc.write_nodes(pathway_adapter.get_nodes())
    bc.write_edges(pathway_adapter.get_edges())
except Exception as e:
    print(f"WARNING: Pathway adapter failed: {e}")
    import traceback
    traceback.print_exc()
finally:
    if 'pathway_adapter' in dir():
        del pathway_adapter
    aggressive_memory_cleanup("Pathway")
    log_adapter_boundary("Pathway", "end")


# side effect
log_adapter_boundary("Side Effect", "start")
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
    print("Memory cleaned up after Side effect adapter")
    log_adapter_boundary("Side Effect", "end")

# ec numbers
log_adapter_boundary("EC", "start")
try:
    ec_adapter = EC(
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        test_mode=TEST_MODE,
        organism=ORGANISM,
    )
    ec_adapter.download_ec_data(cache=CACHE, rxnfp_embedding_path=rxnfp_ec_embedding_path if USE_EMBEDDINGS else None)
    bc.write_nodes(ec_adapter.get_nodes())
    bc.write_edges(ec_adapter.get_edges())
except Exception as e:
    print(f"WARNING: EC adapter failed: {e}")
    import traceback
    traceback.print_exc()
finally:
    if 'ec_adapter' in dir():
        del ec_adapter
    aggressive_memory_cleanup("EC")
    log_adapter_boundary("EC", "end")

# tf-gen
log_adapter_boundary("TFGene", "start")
try:
    tfgene_adapter = TFGene(
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        organism=ORGANISM,
        test_mode=TEST_MODE,
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
    aggressive_memory_cleanup("TFGene")
    log_adapter_boundary("TFGene", "end")


# Write import call and other post-processing
bc.write_import_call()
bc.summary()

print("="*80)
print("Script completed successfully!")
print("="*80)
