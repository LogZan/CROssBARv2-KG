import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
# Add parent directory to Python path
sys.path.insert(0, str(project_root))

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
    Drug
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
import gc  # For garbage collection

bc = BioCypher(biocypher_config_path= str(project_root / "config/biocypher_config.yaml"),
               schema_config_path= str(project_root / "config/schema_config.yaml"),
)

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
]

uniprot_id_type = [
     UniprotIDField.GENE_ENTREZ_ID,
]

print("=" * 80)
print("Starting UniProt adapter...")
print("=" * 80)

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

uniprot_nodes = uniprot_adapter.get_nodes()
uniprot_edges = uniprot_adapter.get_edges()


bc.write_nodes(uniprot_nodes)
bc.write_edges(uniprot_edges)


if export_as_csv:
    uniprot_adapter.export_data_to_csv(path=output_dir_path,
                                    node_data=uniprot_nodes,
                                    edge_data=uniprot_edges)

# Clean up memory after UniProt
del uniprot_adapter, uniprot_nodes, uniprot_edges
gc.collect()

print("✓ UniProt adapter completed")

# PPI
# ppi_adapter = PPI(organism=None, 
#                   output_dir=output_dir_path,
#                   export_csv=export_as_csv,
#                   test_mode=TEST_MODE)

# ppi_adapter.download_ppi_data(cache=CACHE)

# ppi_adapter.process_ppi_data()

# bc.write_edges(ppi_adapter.get_ppi_edges())

# protein domain
# Note: SKIP InterPro due to memory issues - can be run separately
print("=" * 80)
print("SKIPPING InterPro adapter (causes OOM - Out Of Memory)")
print("To run InterPro separately, use a machine with more memory")
print("=" * 80)

# gene ontology
try:
    print("=" * 80)
    print("Starting GO adapter...")
    print("=" * 80)
    go_adapter = GO(
        organism="*", 
        test_mode=TEST_MODE
    )
    go_adapter.download_go_data(cache=CACHE)
    bc.write_nodes(go_adapter.get_go_nodes())
    bc.write_edges(go_adapter.get_go_edges())
    if export_as_csv:
        go_adapter.export_as_csv(path=output_dir_path)
    del go_adapter
    gc.collect()
    print("✓ GO adapter completed")
except Exception as e:
    print(f"WARNING: GO adapter failed: {e}")

# drug
try:
    print("=" * 80)
    print("Starting Drug adapter...")
    print("=" * 80)
    drug_adapter = Drug(
        drugbank_user=drugbank_user, 
        drugbank_passwd=drugbank_passwd,    
        export_csv=export_as_csv, 
        output_dir=output_dir_path,
        test_mode=TEST_MODE
    )
    drug_adapter.download_drug_data(cache=CACHE)
    drug_adapter.process_drug_data()
    bc.write_nodes(drug_adapter.get_drug_nodes())
    bc.write_edges(drug_adapter.get_edges())
    del drug_adapter
    gc.collect()
    print("✓ Drug adapter completed")
except Exception as e:
    print(f"WARNING: Drug adapter failed: {e}")

# compound
try:
    print("=" * 80)
    print("Starting Compound adapter...")
    print("=" * 80)
    compound_adapter = Compound(
        stitch_organism="*",
        export_csv=export_as_csv, 
        output_dir=output_dir_path
    )
    compound_adapter.download_compound_data(cache=CACHE)
    compound_adapter.process_compound_data()
    bc.write_nodes(compound_adapter.get_compound_nodes())
    bc.write_edges(compound_adapter.get_cti_edges())
    del compound_adapter
    gc.collect()
    print("✓ Compound adapter completed")
except Exception as e:
    print(f"WARNING: Compound adapter failed: {e}")

# orthology
try:
    print("=" * 80)
    print("Starting Orthology adapter...")
    print("=" * 80)
    orthology_adapter = Orthology(
        export_csv=export_as_csv, 
        output_dir=output_dir_path,
        test_mode=TEST_MODE
    )
    orthology_adapter.download_orthology_data(cache=CACHE)
    bc.write_edges(orthology_adapter.get_orthology_edges())
    del orthology_adapter
    gc.collect()
    print("✓ Orthology adapter completed")
except Exception as e:
    print(f"WARNING: Orthology adapter failed: {e}")

# disease
try:
    print("=" * 80)
    print("Starting Disease adapter...")
    print("=" * 80)
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
    del disease_adapter
    gc.collect()
    print("✓ Disease adapter completed")
except Exception as e:
    print(f"WARNING: Disease adapter failed: {e}")

# phenotype
try:
    print("=" * 80)
    print("Starting Phenotype adapter...")
    print("=" * 80)
    phenotype_adapter = HPO(
        export_csv=export_as_csv, 
        output_dir=output_dir_path,
        test_mode=TEST_MODE
    )
    phenotype_adapter.download_hpo_data(cache=CACHE)
    bc.write_nodes(phenotype_adapter.get_nodes())
    bc.write_edges(phenotype_adapter.get_edges())
    del phenotype_adapter
    gc.collect()
    print("✓ Phenotype adapter completed")
except Exception as e:
    print(f"WARNING: Phenotype adapter failed: {e}")

# pathway
try:
    print("=" * 80)
    print("Starting Pathway adapter...")
    print("=" * 80)
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
    del pathway_adapter
    gc.collect()
    print("✓ Pathway adapter completed")
except Exception as e:
    print(f"WARNING: Pathway adapter failed: {e}")


# side effect
try:
    print("=" * 80)
    print("Starting Side Effect adapter...")
    print("=" * 80)
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
    del side_effect_adapter
    gc.collect()
    print("✓ Side Effect adapter completed")
except Exception as e:
    print(f"WARNING: Side effect adapter failed: {e}")

# ec numbers
try:
    print("=" * 80)
    print("Starting EC adapter...")
    print("=" * 80)
    ec_adapter = EC(
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        test_mode=TEST_MODE
    )
    ec_adapter.download_ec_data(cache=CACHE)
    bc.write_nodes(ec_adapter.get_nodes())
    bc.write_edges(ec_adapter.get_edges())
    del ec_adapter
    gc.collect()
    print("✓ EC adapter completed")
except Exception as e:
    print(f"WARNING: EC adapter failed: {e}")

# tf-gen
try:
    print("=" * 80)
    print("Starting TFGene adapter...")
    print("=" * 80)
    tfgene_adapter = TFGene(
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        test_mode=TEST_MODE
    )
    tfgene_adapter.download_tfgen_data(cache=CACHE)
    bc.write_edges(tfgene_adapter.get_edges())
    del tfgene_adapter
    gc.collect()
    print("✓ TFGene adapter completed")
except Exception as e:
    print(f"WARNING: TFGene adapter failed: {e}")


# Write import call and other post-processing
print("=" * 80)
print("Writing import call and summary...")
print("=" * 80)
bc.write_import_call()
bc.summary()

print("=" * 80)
print("✓ Script completed successfully!")
print("=" * 80)
