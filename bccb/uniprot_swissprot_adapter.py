"""
UniProt SwissProt Adapter - reads from local JSON file instead of API.
Produces the same output format as uniprot_adapter.py.
Uses ijson for streaming large JSON files.
"""

from time import time
import collections
from typing import Optional, Union, Literal
from collections.abc import Generator
from enum import Enum, EnumMeta, auto
from functools import lru_cache
import pandas as pd
import numpy as np

import os
import h5py
import ijson
import json
import requests

from tqdm import tqdm
from biocypher._logger import logger
from bioregistry import normalize_curie

from pydantic import BaseModel, DirectoryPath, FilePath, HttpUrl, validate_call

# Default embedding directory
DEFAULT_EMBEDDING_DIR = "/GenSIvePFS/users/data/UniProt/UniProtKB_SwissProt/embeddings"

logger.debug(f"Loading module {__name__}.")


class UniprotEnumMeta(EnumMeta):
    def __contains__(cls, item):
        return item in cls.__members__.keys()


class UniprotNodeType(Enum, metaclass=UniprotEnumMeta):
    """
    Node types of the UniProt API represented in this adapter.
    """

    PROTEIN = auto()
    GENE = auto()
    ORGANISM = auto()
    CELLULAR_COMPARTMENT = auto()
    # Extended types from SwissProt
    FUNCTIONAL_ANNOTATION = auto()
    SEQUENCE_FEATURE = auto()
    UNIPROT_DISEASE = auto()


class UniprotNodeField(Enum, metaclass=UniprotEnumMeta):
    """
    Fields of nodes the UniProt API represented in this adapter. Overview of
    uniprot fields: https://www.uniprot.org/help/return_fields
    """

    # core attributes
    LENGTH = "length"
    SUBCELLULAR_LOCATION = "subcellular_location"
    MASS = "mass"
    ORGANISM = "organism_name"
    ORGANISM_ID = "organism_id"
    PROTEIN_NAMES = "protein_name"
    EC = "ec"
    PROTEIN_GENE_NAMES = "gene_names"
    PRIMARY_GENE_NAME = "gene_primary"
    SEQUENCE = "sequence"

    # xref attributes
    ENSEMBL_TRANSCRIPT_IDS = "xref_ensembl"
    PROTEOME = "xref_proteomes"
    ENTREZ_GENE_IDS = "xref_geneid"
    KEGG_IDS = "xref_kegg"

    # not from uniprot REST - computed via pypath
    ENSEMBL_GENE_IDS = "ensembl_gene_ids"

    # embeddings (external files)
    PROTT5_EMBEDDING = "prott5_embedding"
    ESM2_EMBEDDING = "esm2_embedding"
    NT_EMBEDDING = "nt_embedding"

    # === Extended protein properties ===
    # Basic info
    UNIPROT_KB_ID = "uniprot_kb_id"
    SECONDARY_ACCESSIONS = "secondary_accessions"
    PROTEIN_EXISTENCE = "protein_existence"
    ANNOTATION_SCORE = "annotation_score"

    # Version audit
    FIRST_PUBLIC_DATE = "first_public_date"
    LAST_ANNOTATION_UPDATE_DATE = "last_annotation_update_date"
    LAST_SEQUENCE_UPDATE_DATE = "last_sequence_update_date"
    ENTRY_VERSION = "entry_version"
    SEQUENCE_VERSION = "sequence_version"

    # Protein description extended
    RECOMMENDED_SHORT_NAMES = "recommended_short_names"
    EC_NUMBERS = "ec_numbers"
    ALTERNATIVE_NAMES_JSON = "alternative_names_json"
    PROTEIN_FLAG = "protein_flag"

    # Sequence checksums
    SEQUENCE_CRC64 = "sequence_crc64"
    SEQUENCE_MD5 = "sequence_md5"

    # Organism extended
    ORGANISM_COMMON_NAME = "organism_common_name"
    ORGANISM_SYNONYMS = "organism_synonyms"
    ORGANISM_LINEAGE = "organism_lineage"

    # Gene extended
    GENE_ORF_NAMES = "gene_orf_names"
    GENE_ORDERED_LOCUS_NAMES = "gene_ordered_locus_names"
    GENE_LOCATION = "gene_location"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None

    @classmethod
    def get_protein_properties(cls):
        return [
            cls.LENGTH.value,
            cls.MASS.value,
            cls.PROTEIN_NAMES.value,
            cls.PROTEOME.value,
            cls.EC.value,
            cls.ORGANISM.value,
            cls.ORGANISM_ID.value,
            cls.SEQUENCE.value,
            cls.SUBCELLULAR_LOCATION.value,
            cls.PROTT5_EMBEDDING.value,
            cls.ESM2_EMBEDDING.value,
            # Extended properties
            cls.UNIPROT_KB_ID.value,
            cls.SECONDARY_ACCESSIONS.value,
            cls.PROTEIN_EXISTENCE.value,
            cls.ANNOTATION_SCORE.value,
            cls.FIRST_PUBLIC_DATE.value,
            cls.LAST_ANNOTATION_UPDATE_DATE.value,
            cls.ENTRY_VERSION.value,
            cls.SEQUENCE_VERSION.value,
            cls.RECOMMENDED_SHORT_NAMES.value,
            cls.EC_NUMBERS.value,
            cls.ALTERNATIVE_NAMES_JSON.value,
            cls.PROTEIN_FLAG.value,
            cls.SEQUENCE_CRC64.value,
            cls.SEQUENCE_MD5.value,
        ]

    @classmethod
    def get_gene_properties(cls) -> list:
        return [
            cls.PROTEIN_GENE_NAMES.value,
            cls.ENTREZ_GENE_IDS.value,
            cls.KEGG_IDS.value,
            cls.ENSEMBL_TRANSCRIPT_IDS.value,
            cls.ENSEMBL_GENE_IDS.value,
            cls.PRIMARY_GENE_NAME.value,
            cls.NT_EMBEDDING.value,
            cls.GENE_ORF_NAMES.value,
            cls.GENE_ORDERED_LOCUS_NAMES.value,
            cls.GENE_LOCATION.value,
        ]

    @classmethod
    def get_organism_properties(cls) -> list:
        return [
            cls.ORGANISM.value,
            cls.ORGANISM_COMMON_NAME.value,
            cls.ORGANISM_SYNONYMS.value,
            cls.ORGANISM_LINEAGE.value,
        ]

    @classmethod
    def get_split_fields(cls) -> list:
        return [
            cls.PROTEOME.value,
            cls.PROTEIN_GENE_NAMES.value,
            cls.EC.value,
            cls.ENTREZ_GENE_IDS.value,
            cls.ENSEMBL_TRANSCRIPT_IDS.value,
            cls.KEGG_IDS.value,
        ]


class UniprotEdgeType(Enum, metaclass=UniprotEnumMeta):
    """
    Edge types of the UniProt API represented in this adapter.
    """

    PROTEIN_TO_ORGANISM = auto()
    GENE_TO_PROTEIN = auto()
    # Extended edge types from SwissProt
    PROTEIN_TO_KEYWORD = auto()
    PROTEIN_TO_ANNOTATION = auto()
    PROTEIN_TO_FEATURE = auto()
    PROTEIN_TO_DISEASE = auto()
    PROTEIN_INTERACTS_PROTEIN = auto()


class UniprotIDField(Enum, metaclass=UniprotEnumMeta):
    """
    Fields of edges of the UniProt API represented in this adapter. Used to
    assign source and target identifiers in `get_edges()`.
    """

    # default
    PROTEIN_UNIPROT_ACCESSION = auto()
    GENE_ENTREZ_ID = auto()
    ORGANISM_NCBI_TAXONOMY_ID = auto()

    # optional
    GENE_ENSEMBL_GENE_ID = auto()


class UniProtSwissProtModel(BaseModel):
    json_path: str = "/GenSIvePFS/users/data/UniProt/UniProtKB_SwissProt/uniprotkb_reviewed_true_2025_11_04.json"
    organism: Literal["*"] | int | None = "*"
    node_types: Optional[Union[list[UniprotNodeType], None]] = None
    node_fields: Optional[Union[list[UniprotNodeField], None]] = None
    edge_types: Optional[Union[list[UniprotEdgeType], None]] = None
    id_fields: Optional[Union[list[UniprotIDField], None]] = None
    add_prefix: bool = True
    test_mode: bool = False


class UniprotSwissprot:
    """
    Adapter that reads UniProt SwissProt data from a local JSON file.
    Produces the same output format as the original Uniprot adapter.

    Args:
        json_path: path to the SwissProt JSON file.
        organism: organism code in NCBI taxid format, e.g. 9606 for human. "*" for all.
        node_types: `UniprotNodeType` fields that will be included in graph, if it is None, select all fields.
        node_fields: `UniprotNodeField` fields that will be included in graph, if it is None, select all fields.
        edge_types: `UniprotEdgeType` fields that will be included in graph, if it is None, select all fields.
        id_fields: `UniprotIDField` field that will be included in graph as node identifier, if it is None, selects first 3 fields.
        add_prefix: if True, add prefix to database identifiers.
        test_mode: if True, limits amount of output data.
    """

    def __init__(
        self,
        json_path: str = "/GenSIvePFS/users/data/UniProt/UniProtKB_SwissProt/uniprotkb_reviewed_true_2025_11_04.json",
        organism: Optional[Literal["*"] | int | None] = "*",
        node_types: Optional[Union[list[UniprotNodeType], None]] = None,
        node_fields: Optional[Union[list[UniprotNodeField], None]] = None,
        edge_types: Optional[Union[list[UniprotEdgeType], None]] = None,
        id_fields: Optional[Union[list[UniprotIDField], None]] = None,
        add_prefix: Optional[bool] = True,
        test_mode: Optional[bool] = False,
    ):
        model = UniProtSwissProtModel(
            json_path=json_path,
            organism=organism,
            node_types=node_types,
            node_fields=node_fields,
            edge_types=edge_types,
            id_fields=id_fields,
            add_prefix=add_prefix,
            test_mode=test_mode,
        ).model_dump()

        # params
        self.json_path = model["json_path"]
        self.organism = model["organism"]
        self.add_prefix = model["add_prefix"]
        self.test_mode = model["test_mode"]
        self.test_limit = 100

        # provenance
        self.data_source = "uniprot"
        self.data_version = "2025_11"
        self.data_licence = "CC BY 4.0"

        self._configure_fields()

        self._set_node_and_edge_fields(
            node_types=model["node_types"],
            node_fields=model["node_fields"],
            edge_types=model["edge_types"],
        )

        self.set_id_fields(id_fields=model["id_fields"])

        # loading of ligands and receptors sets
        self.ligands = self._read_ligands_set()
        self.receptors = self._read_receptors_set()

        # loading of subcellular locations set
        self.locations = set()

        # gene property name mappings that will be used gene node properties in KG
        self.gene_property_name_mappings = {
            "gene_primary": "gene_symbol",
            "xref_ensembl": "ensembl_transcript_ids",
            "xref_kegg": "kegg_ids",
        }
        # protein property name mappings that will be used protein node properties in KG
        self.protein_property_name_mappings = {"protein_name": "protein_names"}

        # Data storage
        self.data = {}
        self.uniprot_ids = set()
        self.prott5_embedding_df = None
        self.esm2_embedding_df = None
        self.entrez_id_to_nucleotide_transformer_embedding = {}
        
        # Extended data storage for SwissProt
        self.raw_entries = {}  # protein_id -> raw entry for extended processing
        self.annotation_nodes = {}  # annotation_id -> (node_label, props)
        self.feature_nodes = {}  # feature_id -> (node_label, props)
        self.disease_nodes = {}  # disease_id -> props
        self.protein_keywords = {}  # protein_id -> [keyword_ids]
        self.protein_annotations = {}  # protein_id -> [(annotation_id, edge_props)]
        self.protein_features = {}  # protein_id -> [(feature_id, edge_props)]
        self.protein_diseases = {}  # protein_id -> [(disease_id, edge_props)]
        self.protein_interactions = {}  # protein_id -> [(target_protein_id, edge_props)]

    def _read_ligands_set(self) -> set:
        if not os.path.isfile("data/ligands_curated.csv"):
            return set()
        ligand_file = pd.read_csv("data/ligands_curated.csv", header=None)
        return set(ligand_file[0])

    def _read_receptors_set(self) -> set:
        if not os.path.isfile("data/receptors_curated.csv"):
            return set()
        receptor_file = pd.read_csv("data/receptors_curated.csv", header=None)
        return set(receptor_file[0])

    def _stream_entries(self) -> Generator[dict, None, None]:
        """Stream entries from JSON file using ijson."""
        count = 0
        with open(self.json_path, 'rb') as f:
            for entry in ijson.items(f, 'results.item'):
                # Filter by organism if specified
                if self.organism not in ("*", None):
                    taxon_id = entry.get("organism", {}).get("taxonId")
                    if taxon_id != self.organism:
                        continue

                if self.test_mode and count >= self.test_limit:
                    break

                yield entry
                count += 1

    @validate_call
    def download_uniprot_data(
        self,
        cache: bool = False,
        debug: bool = False,
        retries: int = 3,
        prott5_embedding_output_path: FilePath | None = None,
        esm2_embedding_path: FilePath | None = None,
        nucleotide_transformer_embedding_path: FilePath | None = None,
    ):
        """
        Load uniprot data from local JSON file.

        Args:
            cache: unused, kept for API compatibility
            debug: unused, kept for API compatibility
            retries: unused, kept for API compatibility
            prott5_embedding_output_path: path to ProtT5 embeddings h5 file
            esm2_embedding_path: path to ESM2 embeddings h5 file
            nucleotide_transformer_embedding_path: path to NT embeddings h5 file
        """
        # Set adapter-specific cache directory for pypath utilities (e.g. mapping)
        from . import cache_config
        cache_config.set_adapter_cache('uniprot')
        
        self._load_json_data()

        # Load embeddings if requested
        if UniprotNodeField.PROTT5_EMBEDDING.value in self.node_fields and not self.test_mode:
            self.data[UniprotNodeField.PROTT5_EMBEDDING.value] = {}
            self.download_prott5_embeddings(
                prott5_embedding_output_path=prott5_embedding_output_path
            )

        if UniprotNodeField.ESM2_EMBEDDING.value in self.node_fields and not self.test_mode:
            self.data[UniprotNodeField.ESM2_EMBEDDING.value] = {}
            self.retrieve_esm2_embeddings(esm2_embedding_path)

        if UniprotNodeField.NT_EMBEDDING.value in self.node_fields and not self.test_mode:
            self.data[UniprotNodeField.NT_EMBEDDING.value] = {}
            self.retrieve_nucleotide_transformer_embeddings(
                nucleotide_transformer_embedding_path
            )

    def _load_json_data(self):
        """Load and parse JSON data from file."""
        logger.info(f"Loading SwissProt data from {self.json_path}...")
        t0 = time()

        # Initialize data storage for each field
        for field in self.node_fields:
            if field not in [
                UniprotNodeField.PROTT5_EMBEDDING.value,
                UniprotNodeField.ESM2_EMBEDDING.value,
                UniprotNodeField.NT_EMBEDDING.value,
            ]:
                self.data[field] = {}

        # Add ensembl gene ids storage
        self.data[UniprotNodeField.ENSEMBL_GENE_IDS.value] = {}

        # Stream and process entries
        for entry in tqdm(self._stream_entries(), desc="Loading JSON entries"):
            self._process_entry(entry)

        t1 = time()
        msg = f"Loaded {len(self.uniprot_ids)} UniProt entries in {round((t1-t0) / 60, 2)} mins."
        logger.info(msg)

    def _process_entry(self, entry: dict):
        """Process a single JSON entry and extract all fields."""
        protein_id = entry.get("primaryAccession")
        if not protein_id:
            return

        self.uniprot_ids.add(protein_id)

        # Extract sequence info
        seq = entry.get("sequence", {})
        if UniprotNodeField.LENGTH.value in self.node_fields:
            self.data[UniprotNodeField.LENGTH.value][protein_id] = seq.get("length")

        if UniprotNodeField.MASS.value in self.node_fields:
            self.data[UniprotNodeField.MASS.value][protein_id] = seq.get("molWeight")

        if UniprotNodeField.SEQUENCE.value in self.node_fields:
            self.data[UniprotNodeField.SEQUENCE.value][protein_id] = seq.get("value")

        if UniprotNodeField.SEQUENCE_CRC64.value in self.node_fields:
            self.data[UniprotNodeField.SEQUENCE_CRC64.value][protein_id] = seq.get("crc64")

        if UniprotNodeField.SEQUENCE_MD5.value in self.node_fields:
            self.data[UniprotNodeField.SEQUENCE_MD5.value][protein_id] = seq.get("md5")

        # Extract organism info
        org = entry.get("organism", {})
        if UniprotNodeField.ORGANISM.value in self.node_fields:
            self.data[UniprotNodeField.ORGANISM.value][protein_id] = org.get("scientificName")

        if UniprotNodeField.ORGANISM_ID.value in self.node_fields:
            self.data[UniprotNodeField.ORGANISM_ID.value][protein_id] = org.get("taxonId")

        if UniprotNodeField.ORGANISM_COMMON_NAME.value in self.node_fields:
            self.data[UniprotNodeField.ORGANISM_COMMON_NAME.value][protein_id] = org.get("commonName")

        if UniprotNodeField.ORGANISM_SYNONYMS.value in self.node_fields:
            synonyms = org.get("synonyms", [])
            if synonyms:
                self.data[UniprotNodeField.ORGANISM_SYNONYMS.value][protein_id] = synonyms

        if UniprotNodeField.ORGANISM_LINEAGE.value in self.node_fields:
            lineage = org.get("lineage", [])
            if lineage:
                self.data[UniprotNodeField.ORGANISM_LINEAGE.value][protein_id] = lineage

        # Extract protein description
        desc = entry.get("proteinDescription", {})
        rec_name = desc.get("recommendedName", {})

        if UniprotNodeField.PROTEIN_NAMES.value in self.node_fields:
            full_name = rec_name.get("fullName", {}).get("value")
            # Also collect alternative names
            alt_names = []
            for alt in desc.get("alternativeNames", []):
                alt_full = alt.get("fullName", {}).get("value")
                if alt_full:
                    alt_names.append(alt_full)
            
            if full_name:
                if alt_names:
                    self.data[UniprotNodeField.PROTEIN_NAMES.value][protein_id] = [full_name] + alt_names
                else:
                    self.data[UniprotNodeField.PROTEIN_NAMES.value][protein_id] = full_name

        if UniprotNodeField.EC.value in self.node_fields:
            ec_nums = rec_name.get("ecNumbers", [])
            ec_values = [e.get("value") for e in ec_nums if e.get("value")]
            # Also get EC from alternative names
            for alt in desc.get("alternativeNames", []):
                for ec in alt.get("ecNumbers", []):
                    if ec.get("value"):
                        ec_values.append(ec.get("value"))
            if ec_values:
                self.data[UniprotNodeField.EC.value][protein_id] = ec_values if len(ec_values) > 1 else ec_values[0]

        # Extract gene info
        genes = entry.get("genes", [])
        if genes:
            gene_names = []
            primary_gene = None
            orf_names = []
            ordered_locus_names = []

            for gene_data in genes:
                if gene_data.get("geneName"):
                    gene_name = gene_data["geneName"].get("value")
                    if gene_name:
                        gene_names.append(gene_name)
                        if not primary_gene:
                            primary_gene = gene_name
                # Also add synonyms
                for syn in gene_data.get("synonyms", []):
                    if syn.get("value"):
                        gene_names.append(syn.get("value"))
                # Extract orfNames
                for orf in gene_data.get("orfNames", []):
                    if orf.get("value"):
                        orf_names.append(orf.get("value"))
                # Extract orderedLocusNames
                for oln in gene_data.get("orderedLocusNames", []):
                    if oln.get("value"):
                        ordered_locus_names.append(oln.get("value"))

            if UniprotNodeField.PROTEIN_GENE_NAMES.value in self.node_fields and gene_names:
                self.data[UniprotNodeField.PROTEIN_GENE_NAMES.value][protein_id] = gene_names if len(gene_names) > 1 else gene_names[0]

            if UniprotNodeField.PRIMARY_GENE_NAME.value in self.node_fields and primary_gene:
                self.data[UniprotNodeField.PRIMARY_GENE_NAME.value][protein_id] = primary_gene

            if UniprotNodeField.GENE_ORF_NAMES.value in self.node_fields and orf_names:
                self.data[UniprotNodeField.GENE_ORF_NAMES.value][protein_id] = orf_names if len(orf_names) > 1 else orf_names[0]

            if UniprotNodeField.GENE_ORDERED_LOCUS_NAMES.value in self.node_fields and ordered_locus_names:
                self.data[UniprotNodeField.GENE_ORDERED_LOCUS_NAMES.value][protein_id] = ordered_locus_names if len(ordered_locus_names) > 1 else ordered_locus_names[0]

        # Extract gene locations
        if UniprotNodeField.GENE_LOCATION.value in self.node_fields:
            gene_locs = entry.get("geneLocations", [])
            if gene_locs:
                loc_data = []
                for loc in gene_locs:
                    loc_info = {
                        "value": loc.get("value"),
                        "type": loc.get("geneEncodingType")
                    }
                    loc_data.append(loc_info)
                self.data[UniprotNodeField.GENE_LOCATION.value][protein_id] = loc_data if len(loc_data) > 1 else loc_data[0]

        # === Extended protein properties extraction ===
        
        # Basic info
        if UniprotNodeField.UNIPROT_KB_ID.value in self.node_fields:
            self.data[UniprotNodeField.UNIPROT_KB_ID.value][protein_id] = entry.get("uniProtkbId")
        
        if UniprotNodeField.SECONDARY_ACCESSIONS.value in self.node_fields:
            sec_acc = entry.get("secondaryAccessions", [])
            if sec_acc:
                self.data[UniprotNodeField.SECONDARY_ACCESSIONS.value][protein_id] = sec_acc
        
        if UniprotNodeField.PROTEIN_EXISTENCE.value in self.node_fields:
            self.data[UniprotNodeField.PROTEIN_EXISTENCE.value][protein_id] = entry.get("proteinExistence")
        
        if UniprotNodeField.ANNOTATION_SCORE.value in self.node_fields:
            self.data[UniprotNodeField.ANNOTATION_SCORE.value][protein_id] = entry.get("annotationScore")
        
        # Version audit info
        entry_audit = entry.get("entryAudit", {})
        if UniprotNodeField.FIRST_PUBLIC_DATE.value in self.node_fields:
            self.data[UniprotNodeField.FIRST_PUBLIC_DATE.value][protein_id] = entry_audit.get("firstPublicDate")

        if UniprotNodeField.LAST_ANNOTATION_UPDATE_DATE.value in self.node_fields:
            self.data[UniprotNodeField.LAST_ANNOTATION_UPDATE_DATE.value][protein_id] = entry_audit.get("lastAnnotationUpdateDate")

        if UniprotNodeField.LAST_SEQUENCE_UPDATE_DATE.value in self.node_fields:
            self.data[UniprotNodeField.LAST_SEQUENCE_UPDATE_DATE.value][protein_id] = entry_audit.get("lastSequenceUpdateDate")

        if UniprotNodeField.ENTRY_VERSION.value in self.node_fields:
            self.data[UniprotNodeField.ENTRY_VERSION.value][protein_id] = entry_audit.get("entryVersion")

        if UniprotNodeField.SEQUENCE_VERSION.value in self.node_fields:
            self.data[UniprotNodeField.SEQUENCE_VERSION.value][protein_id] = entry_audit.get("sequenceVersion")
        
        # Extended protein description
        if UniprotNodeField.RECOMMENDED_SHORT_NAMES.value in self.node_fields:
            short_names = []
            rec_name = desc.get("recommendedName", {})
            for sn in rec_name.get("shortNames", []):
                if sn.get("value"):
                    short_names.append(sn.get("value"))
            if short_names:
                self.data[UniprotNodeField.RECOMMENDED_SHORT_NAMES.value][protein_id] = short_names
        
        if UniprotNodeField.EC_NUMBERS.value in self.node_fields:
            ec_nums = []
            rec_name = desc.get("recommendedName", {})
            for ec in rec_name.get("ecNumbers", []):
                if ec.get("value"):
                    ec_nums.append(ec.get("value"))
            # Also get from alternative names
            for alt in desc.get("alternativeNames", []):
                for ec in alt.get("ecNumbers", []):
                    if ec.get("value"):
                        ec_nums.append(ec.get("value"))
            if ec_nums:
                self.data[UniprotNodeField.EC_NUMBERS.value][protein_id] = ec_nums
        
        if UniprotNodeField.ALTERNATIVE_NAMES_JSON.value in self.node_fields:
            alt_names = desc.get("alternativeNames", [])
            if alt_names:
                import json
                self.data[UniprotNodeField.ALTERNATIVE_NAMES_JSON.value][protein_id] = json.dumps(alt_names)
        
        if UniprotNodeField.PROTEIN_FLAG.value in self.node_fields:
            flag = desc.get("flag")
            if flag:
                self.data[UniprotNodeField.PROTEIN_FLAG.value][protein_id] = flag

        # Extract cross-references
        xrefs = entry.get("uniProtKBCrossReferences", [])
        self._extract_cross_refs(protein_id, xrefs)

        # Extract subcellular location from comments
        if UniprotNodeField.SUBCELLULAR_LOCATION.value in self.node_fields:
            locations = self._extract_subcellular_locations(entry.get("comments", []))
            if locations:
                self.data[UniprotNodeField.SUBCELLULAR_LOCATION.value][protein_id] = locations
                for loc in locations:
                    self.locations.add(loc)

        # Extended data extraction for SwissProt
        self._extract_extended_data(protein_id, entry)

    def _extract_cross_refs(self, protein_id: str, xrefs: list):
        """Extract cross-references by database type."""
        ensembl_ids = []
        proteome_ids = []
        entrez_ids = []
        kegg_ids = []

        for xref in xrefs:
            db = xref.get("database")
            xref_id = xref.get("id")

            if db == "Ensembl" and UniprotNodeField.ENSEMBL_TRANSCRIPT_IDS.value in self.node_fields:
                ensembl_ids.append(xref_id)

            elif db == "Proteomes" and UniprotNodeField.PROTEOME.value in self.node_fields:
                proteome_ids.append(xref_id)

            elif db == "GeneID" and UniprotNodeField.ENTREZ_GENE_IDS.value in self.node_fields:
                entrez_ids.append(xref_id)

            elif db == "KEGG" and UniprotNodeField.KEGG_IDS.value in self.node_fields:
                kegg_ids.append(xref_id)

        if ensembl_ids:
            self.data[UniprotNodeField.ENSEMBL_TRANSCRIPT_IDS.value][protein_id] = ensembl_ids if len(ensembl_ids) > 1 else ensembl_ids[0]
            # Try to extract ENSG IDs
            ensg_ids = self._find_ensg_from_enst(ensembl_ids)
            if ensg_ids:
                self.data[UniprotNodeField.ENSEMBL_GENE_IDS.value][protein_id] = ensg_ids

        if proteome_ids:
            self.data[UniprotNodeField.PROTEOME.value][protein_id] = proteome_ids if len(proteome_ids) > 1 else proteome_ids[0]

        if entrez_ids:
            # Take first entrez ID as per original adapter
            self.data[UniprotNodeField.ENTREZ_GENE_IDS.value][protein_id] = entrez_ids[0]

        if kegg_ids:
            # Remove prefix from KEGG IDs (e.g., "hsa:1234" -> "1234")
            clean_kegg = [kid.split(":")[-1] if ":" in kid else kid for kid in kegg_ids]
            self.data[UniprotNodeField.KEGG_IDS.value][protein_id] = clean_kegg if len(clean_kegg) > 1 else clean_kegg[0]

    def _extract_subcellular_locations(self, comments: list) -> list:
        """Extract subcellular locations from comments."""
        locations = []
        for comment in comments:
            if comment.get("commentType") == "SUBCELLULAR LOCATION":
                for loc_data in comment.get("subcellularLocations", []):
                    loc = loc_data.get("location", {}).get("value")
                    if loc:
                        locations.append(loc.replace("'", "").strip())
        return locations

    def _find_ensg_from_enst(self, enst_list):
        """
        Take ensembl transcript ids, return ensembl gene ids.
        Falls back to extracting from the IDs if pypath mapping not available.
        """
        enst_list = self._ensure_iterable(enst_list)
        enst_list = [enst.split(" [")[0].split(".")[0] for enst in enst_list]

        ensg_ids = set()
        try:
            from pypath.utils import mapping
            for enst_id in enst_list:
                ensg_id = list(
                    mapping.map_name(
                        enst_id, "enst_biomart", "ensg_biomart"
                    )
                )
                ensg_id = ensg_id[0] if ensg_id else None
                if ensg_id:
                    ensg_ids.add(ensg_id)
        except Exception:
            # pypath not available or has dependency issues, skip ENSG mapping
            pass

        ensg_ids = list(ensg_ids)
        if len(ensg_ids) == 1:
            ensg_ids = ensg_ids[0]
        return ensg_ids if ensg_ids else None

    @validate_call
    def download_prott5_embeddings(
        self,
        prott5_embedding_output_path: FilePath | str | None = None
    ):
        """
        Downloads ProtT5 embedding from uniprot website.
        If the file exists in the defined file path, then directly reads it.

        Args:
            prott5_embedding_output_path: Path to the h5 file. Defaults to embedding dir.
        """
        url: HttpUrl = (
            "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/embeddings/uniprot_sprot/per-protein.h5"
        )
        
        if prott5_embedding_output_path:
            full_path = str(prott5_embedding_output_path)
        else:
            # Create default embedding directory if it doesn't exist
            os.makedirs(DEFAULT_EMBEDDING_DIR, exist_ok=True)
            full_path = os.path.join(DEFAULT_EMBEDDING_DIR, "prott5_protein_embeddings.h5")

        if not os.path.isfile(full_path):
            logger.info(f"Downloading ProtT5 embeddings to {full_path}...")
            try:
                with requests.get(url, stream=True) as response:
                    response.raise_for_status()
                    with open(full_path, "wb") as f:
                        for chunk in response.iter_content(512 * 1024):
                            if chunk:
                                f.write(chunk)
                logger.info("ProtT5 embeddings downloaded successfully.")
            except Exception as e:
                logger.error(f"Failed to download ProtT5 embeddings: {e}")
                return
        else:
            logger.info(f"ProtT5 Embedding file exists. Reading from {full_path}...")

        df_list = []
        with h5py.File(full_path, "r") as file:
            for uniprot_id, embedding in file.items():
                embedding = np.array(embedding).astype(np.float16)
                count = np.count_nonzero(~np.isnan(embedding))
                if (
                    self.organism not in ("*", None)
                    and uniprot_id in self.uniprot_ids
                    and count == 1024
                ):
                    df_list.append((uniprot_id, embedding))
                elif count == 1024:
                    df_list.append((uniprot_id, embedding))

        self.prott5_embedding_df = pd.DataFrame(df_list, columns=['uniprot_id', 'embedding'])
        del df_list

    @validate_call
    def retrieve_esm2_embeddings(
        self,
        esm2_embedding_path: FilePath | str | None = None
    ) -> None:
        """
        Load ESM2 embeddings from h5 file.
        Note: ESM2 embeddings are not available for download from UniProt,
        so this method requires a local file path.
        """
        if esm2_embedding_path:
            full_path = str(esm2_embedding_path)
        else:
            # Default path for ESM2 embeddings
            full_path = os.path.join(DEFAULT_EMBEDDING_DIR, "esm2_embeddings.h5")

        if not os.path.isfile(full_path):
            logger.warning(f"ESM2 embedding file not found at {full_path}. Skipping ESM2 embeddings.")
            return

        logger.info(f"Loading ESM2 embeddings from {full_path}...")

        df_list = []
        with h5py.File(full_path, "r") as file:
            for uniprot_id, embedding in file.items():
                embedding = np.array(embedding).astype(np.float16)
                count = np.count_nonzero(~np.isnan(embedding))
                if (
                    self.organism not in ("*", None)
                    and uniprot_id in self.uniprot_ids
                    and count == 1280
                ):
                    df_list.append((uniprot_id, embedding))
                elif count == 1280:
                    df_list.append((uniprot_id, embedding))

        self.esm2_embedding_df = pd.DataFrame(df_list, columns=['uniprot_id', 'embedding'])
        del df_list

    def retrieve_nucleotide_transformer_embeddings(
        self,
        nucleotide_transformer_embedding_path: FilePath | str | None = None
    ) -> None:
        """
        Load Nucleotide Transformer embeddings from h5 file.
        Note: NT embeddings require a local file path.
        """
        if nucleotide_transformer_embedding_path:
            full_path = str(nucleotide_transformer_embedding_path)
        else:
            # Default path for NT embeddings
            full_path = os.path.join(DEFAULT_EMBEDDING_DIR, "nt_embeddings.h5")

        if not os.path.isfile(full_path):
            logger.warning(f"Nucleotide Transformer embedding file not found at {full_path}. Skipping NT embeddings.")
            return

        logger.info(f"Loading Nucleotide Transformer embeddings from {full_path}...")

        self.entrez_id_to_nucleotide_transformer_embedding = {}
        with h5py.File(full_path, "r") as file:
            for entrez_id, embedding in file.items():
                embedding = np.array(embedding).astype(np.float16)
                self.entrez_id_to_nucleotide_transformer_embedding[entrez_id] = embedding

    @validate_call
    def _get_ligand_or_receptor(self, uniprot_id: str):
        """Tell if UniProt protein node is a L, R or nothing."""
        uniprot_id = uniprot_id[8:] if uniprot_id.startswith("uniprot:") else uniprot_id

        if uniprot_id in self.ligands:
            return "ligand"
        return "receptor" if uniprot_id in self.receptors else "protein"

    @validate_call
    def get_nodes(
        self,
        ligand_or_receptor: bool = False,
        protein_label: str = "protein",
        gene_label: str = "gene",
        organism_label: str = "organism"
    ) -> Generator[tuple[str, str, dict]]:
        """
        Yield nodes (protein, gene, organism) from SwissProt data.
        """
        if ligand_or_receptor and (not self.ligands or not self.receptors):
            raise ValueError(
                "No ligands or receptors found in the 'data' directory. "
                "Please set ligand_or_receptor to False or add the files."
            )

        logger.info(
            "Preparing UniProt nodes of the types "
            f"{[type.name for type in self.node_types]}."
        )

        seen_organisms = set()
        seen_genes = set()

        for protein_id, all_props in self._reformat_and_filter_proteins():

            if UniprotNodeType.PROTEIN in self.node_types:
                protein_props = self._get_protein_properties(all_props, protein_id.split(":")[-1])

                if ligand_or_receptor:
                    lor = self._get_ligand_or_receptor(protein_id)
                    yield (protein_id, lor, protein_props)
                else:
                    yield (protein_id, protein_label, protein_props)

            if UniprotNodeType.GENE in self.node_types:
                gene_list = self._get_gene(all_props)
                for gene_id, gene_props in gene_list:
                    if gene_id and gene_id not in seen_genes:
                        seen_genes.add(gene_id)
                        yield (gene_id, gene_label, gene_props)

            if UniprotNodeType.ORGANISM in self.node_types:
                organism_id, organism_props = self._get_organism(all_props)
                if organism_id and organism_id not in seen_organisms:
                    seen_organisms.add(organism_id)
                    yield (organism_id, organism_label, organism_props)

    @validate_call
    def get_edges(
        self,
        gene_to_protein_label: str = "Gene_encodes_protein",
        protein_to_organism_label: str = "Protein_belongs_to_organism"
    ) -> Generator[tuple[None, str, str, str, dict]]:
        """Get edges from SwissProt data."""
        logger.info(
            "Preparing UniProt edges of the types "
            f"{[type.name for type in self.edge_types]}."
        )

        edge_list = []
        properties = {
            "source": self.data_source,
            "licence": self.data_licence,
            "version": self.data_version,
        }

        for protein in tqdm(self.uniprot_ids, desc="Getting edges"):
            protein_id = self.add_prefix_to_id("uniprot", protein)

            if UniprotEdgeType.GENE_TO_PROTEIN in self.edge_types:
                type_dict = {
                    UniprotNodeField.ENTREZ_GENE_IDS.value: "ncbigene",
                    UniprotNodeField.ENSEMBL_GENE_IDS.value: "ensembl",
                }

                # Find preferred identifier for gene
                if UniprotIDField.GENE_ENTREZ_ID in self.id_fields:
                    id_type = UniprotNodeField.ENTREZ_GENE_IDS.value
                elif UniprotIDField.GENE_ENSEMBL_GENE_ID in self.id_fields:
                    id_type = UniprotNodeField.ENSEMBL_GENE_IDS.value
                else:
                    id_type = UniprotNodeField.ENTREZ_GENE_IDS.value

                genes = self.data.get(id_type, {}).get(protein)
                if genes:
                    genes = self._ensure_iterable(genes)
                    for gene in genes:
                        if not gene:
                            continue
                        gene_id = self.add_prefix_to_id(type_dict[id_type], str(gene))
                        edge_list.append((
                            None,
                            gene_id,
                            protein_id,
                            gene_to_protein_label,
                            properties,
                        ))

            if UniprotEdgeType.PROTEIN_TO_ORGANISM in self.edge_types:
                organism_id = self.data.get(UniprotNodeField.ORGANISM_ID.value, {}).get(protein)
                if organism_id:
                    organism_id = self.add_prefix_to_id("ncbitaxon", str(organism_id))
                    edge_list.append((
                        None,
                        protein_id,
                        organism_id,
                        protein_to_organism_label,
                        properties,
                    ))

        if edge_list:
            return edge_list

    def _reformat_and_filter_proteins(self):
        """
        For each uniprot id, select desired fields and reformat to give a tuple
        containing id and properties. Yield a tuple for each protein.
        """
        for protein in tqdm(self.uniprot_ids, desc="Processing proteins"):
            protein_id = self.add_prefix_to_id("uniprot", protein)
            _props = {arg: self.data.get(arg, {}).get(protein) for arg in self.node_fields}
            yield protein_id, _props

    @validate_call
    def _get_gene(self, all_props: dict) -> list:
        """
        Get gene node representation from UniProt data per protein.
        """
        if not (
            UniprotNodeField.PROTEIN_GENE_NAMES.value in all_props.keys()
            and UniprotNodeField.ENTREZ_GENE_IDS.value in all_props.keys()
        ):
            return []

        # Find preferred identifier for gene
        if UniprotIDField.GENE_ENTREZ_ID in self.id_fields:
            id_type = UniprotNodeField.ENTREZ_GENE_IDS.value
        elif UniprotIDField.GENE_ENSEMBL_GENE_ID in self.id_fields:
            id_type = UniprotNodeField.ENSEMBL_GENE_IDS.value
        else:
            id_type = UniprotNodeField.ENTREZ_GENE_IDS.value

        gene_raw = all_props.get(id_type)
        if not gene_raw:
            return []

        type_dict = {
            UniprotNodeField.ENTREZ_GENE_IDS.value: "ncbigene",
            UniprotNodeField.ENSEMBL_GENE_IDS.value: "ensembl",
        }

        genes = self._ensure_iterable(gene_raw)
        gene_props = {}

        for k in all_props.keys():
            if k not in self.gene_properties:
                continue

            gene_props[
                (
                    self.gene_property_name_mappings[k]
                    if self.gene_property_name_mappings.get(k)
                    else k.replace(" ", "_").replace("-", "_").lower()
                )
            ] = all_props[k]

            if k == UniprotNodeField.NT_EMBEDDING.value:
                if genes and self.entrez_id_to_nucleotide_transformer_embedding.get(str(genes[0])) is not None:
                    gene_props[k] = [str(emb) for emb in self.entrez_id_to_nucleotide_transformer_embedding[str(genes[0])]]

        gene_props["source"] = self.data_source
        gene_props["licence"] = self.data_licence
        gene_props["version"] = self.data_version

        gene_list = []
        for gene in genes:
            gene_id = self.add_prefix_to_id(type_dict[id_type], str(gene))
            gene_list.append((gene_id, gene_props.copy()))

        return gene_list

    @validate_call
    def _get_organism(self, all_props: dict):
        """Get organism node representation."""
        organism_props = {}

        org_id = all_props.get(UniprotNodeField.ORGANISM_ID.value)
        if not org_id:
            return None, None

        organism_id = self.add_prefix_to_id("ncbitaxon", str(org_id))

        for k in all_props.keys():
            if k in self.organism_properties:
                organism_props[k] = all_props[k]

        organism_props["source"] = self.data_source
        organism_props["licence"] = self.data_licence
        organism_props["version"] = self.data_version

        return organism_id, organism_props

    @validate_call
    def _get_protein_properties(self, all_props: dict, protein_id: str) -> dict:
        """Extract protein properties from all_props."""
        protein_props = {}

        for k in all_props.keys():
            if k not in self.protein_properties:
                continue

            elif k == UniprotNodeField.PROTEIN_NAMES.value:
                names = self._ensure_iterable(all_props[k]) if all_props[k] else None
                protein_props["primary_protein_name"] = names[0] if names else None
                protein_props[
                    (
                        self.protein_property_name_mappings[k]
                        if self.protein_property_name_mappings.get(k)
                        else k.replace(" ", "_").replace("-", "_")
                    )
                ] = all_props[k]

            elif k == UniprotNodeField.PROTT5_EMBEDDING.value:
                if self.prott5_embedding_df is not None:
                    res = self.prott5_embedding_df[self.prott5_embedding_df["uniprot_id"] == protein_id]["embedding"]
                    if not res.empty:
                        protein_props[k.replace(" ", "_").replace("-", "_")] = [str(emb) for emb in res.values[0]]

            elif k == UniprotNodeField.ESM2_EMBEDDING.value:
                if self.esm2_embedding_df is not None:
                    res = self.esm2_embedding_df[self.esm2_embedding_df["uniprot_id"] == protein_id]["embedding"]
                    if not res.empty:
                        protein_props[k.replace(" ", "_").replace("-", "_")] = [str(emb) for emb in res.values[0]]

            else:
                protein_props[
                    (
                        self.protein_property_name_mappings[k]
                        if self.protein_property_name_mappings.get(k)
                        else k.replace(" ", "_").replace("-", "_")
                    )
                ] = all_props[k]

        protein_props["source"] = self.data_source
        protein_props["licence"] = self.data_licence
        protein_props["version"] = self.data_version

        return protein_props

    @lru_cache
    @validate_call
    def add_prefix_to_id(
        self, prefix: str = None, identifier: str = None, sep: str = ":"
    ) -> str:
        """Adds prefix to database id."""
        if self.add_prefix and identifier:
            return normalize_curie(prefix + sep + identifier)
        return identifier

    def _configure_fields(self):
        """Configure field lists."""
        self.split_fields = UniprotNodeField.get_split_fields()
        self.protein_properties = UniprotNodeField.get_protein_properties()
        self.gene_properties = UniprotNodeField.get_gene_properties()
        self.organism_properties = UniprotNodeField.get_organism_properties()

    def _set_node_and_edge_fields(
        self,
        node_types,
        node_fields,
        edge_types,
    ):
        """Set node and edge fields."""
        # Ensure computation of ENSGs
        if (
            node_fields
            and UniprotNodeField.ENSEMBL_GENE_IDS in node_fields
            and UniprotNodeField.ENSEMBL_TRANSCRIPT_IDS not in node_fields
        ):
            node_fields.append(UniprotNodeField.ENSEMBL_TRANSCRIPT_IDS)

        self.node_types = node_types or list(UniprotNodeType)
        self.node_fields = [field.value for field in node_fields] if node_fields else [field.value for field in UniprotNodeField]
        self.edge_types = edge_types or list(UniprotEdgeType)

    def set_id_fields(self, id_fields):
        """Set ID fields."""
        self.id_fields = id_fields or list(UniprotIDField)[:3]

    def _ensure_iterable(self, value):
        """Ensure value is iterable (list)."""
        if value is None:
            return []
        return [value] if isinstance(value, (str, int)) else list(value)

    @validate_call
    def export_data_to_csv(
        self,
        node_data: Generator[tuple[str, str, dict]] = None,
        edge_data: Generator[tuple[None, str, str, str, dict]] = None,
        path: DirectoryPath | None = None,
    ) -> None:
        """
        Save node and edge data to csv.
            node_data: output of `get_nodes()` function
            edge_data: output of `get_edges()` function
            path: Directory to save the output csv file
        """
        if node_data:
            logger.debug("Saving uniprot node data as csv")
            node_types_dict = collections.defaultdict(list)
            for _id, _type, props in node_data:
                _dict = {"id": _id} | props
                node_types_dict[_type].append(_dict)

            for _type, values in node_types_dict.items():
                df = pd.DataFrame.from_records(values)
                if path:
                    full_path = os.path.join(path, f"{_type.capitalize()}.csv")
                else:
                    full_path = os.path.join(os.getcwd(), f"{_type.capitalize()}.csv")

                df.to_csv(full_path, index=False)
                logger.info(f"{_type.capitalize()} data is written: {full_path}")

        if edge_data:
            logger.debug("Saving uniprot edge data as csv")
            edge_types_dict = collections.defaultdict(list)
            for _, source_id, target_id, _type, props in edge_data:
                _dict = {"source_id": source_id, "target_id": target_id} | props
                edge_types_dict[_type].append(_dict)

            for _type, values in edge_types_dict.items():
                df = pd.DataFrame.from_records(values)
                if path:
                    full_path = os.path.join(path, f"{_type.capitalize()}.csv")
                else:
                    full_path = os.path.join(os.getcwd(), f"{_type.capitalize()}.csv")

                df.to_csv(full_path, index=False)
                logger.info(f"{_type.capitalize()} data is written: {full_path}")

    # ========================================
    # Extended SwissProt Data Extraction
    # ========================================

    def _extract_extended_data(self, protein_id: str, entry: dict):
        """Extract extended data from SwissProt entry (comments, features, keywords, diseases)."""
        
        # Extract keywords
        if UniprotEdgeType.PROTEIN_TO_KEYWORD in self.edge_types:
            self._extract_keywords(protein_id, entry)
        
        # Extract functional annotations (comments)
        if UniprotEdgeType.PROTEIN_TO_ANNOTATION in self.edge_types:
            self._extract_annotations(protein_id, entry)
        
        # Extract sequence features
        if UniprotEdgeType.PROTEIN_TO_FEATURE in self.edge_types:
            self._extract_features(protein_id, entry)
        
        # Extract diseases
        if UniprotEdgeType.PROTEIN_TO_DISEASE in self.edge_types:
            self._extract_diseases(protein_id, entry)
        
        # Extract protein-protein interactions
        if UniprotEdgeType.PROTEIN_INTERACTS_PROTEIN in self.edge_types:
            self._extract_interactions(protein_id, entry)

    def _extract_keywords(self, protein_id: str, entry: dict):
        """Extract keywords from entry with full data (id, category, name)."""
        keywords = entry.get("keywords", [])
        keyword_data = []
        for kw in keywords:
            kw_id = kw.get("id")
            if kw_id:
                keyword_data.append({
                    "id": kw_id,
                    "category": kw.get("category"),
                    "name": kw.get("name")
                })
        if keyword_data:
            self.protein_keywords[protein_id] = keyword_data

    def _extract_annotations(self, protein_id: str, entry: dict):
        """Dynamically extract all comment types as functional annotations."""
        comments = entry.get("comments", [])
        annotations = []
        
        for idx, comment in enumerate(comments):
            comment_type = comment.get("commentType", "UNKNOWN")
            
            # Skip INTERACTION (handled separately) and SUBCELLULAR LOCATION (handled in basic extraction)
            if comment_type in ("INTERACTION",):
                continue
            
            # Generate annotation ID
            annotation_id = f"{protein_id}_{comment_type}_{idx}"
            
            # Determine node label from comment type
            node_label = comment_type.lower().replace(" ", "_").replace("-", "_") + "_annotation"
            
            # Extract text content
            text = self._get_comment_text(comment)
            
            # Edge properties
            edge_props = {
                "text": text,
            }
            
            # Handle specific comment types for additional properties
            if comment_type == "CATALYTIC ACTIVITY":
                reaction = comment.get("reaction", {})
                edge_props["reaction_name"] = reaction.get("name")
                ec_entries = reaction.get("ecNumber")
                if ec_entries:
                    edge_props["rhea_id"] = None
                # Extract cross-references
                xrefs = reaction.get("reactionCrossReferences", [])
                rhea_ids = [x.get("id") for x in xrefs if x.get("database") == "Rhea"]
                chebi_ids = [x.get("id") for x in xrefs if x.get("database") == "ChEBI"]
                if rhea_ids:
                    edge_props["rhea_id"] = rhea_ids[0]
                if chebi_ids:
                    edge_props["chebi_ids"] = chebi_ids
            
            elif comment_type == "COFACTOR":
                cofactors = comment.get("cofactors", [])
                if cofactors:
                    edge_props["cofactor_name"] = cofactors[0].get("name")
                    cf_xrefs = cofactors[0].get("cofactorCrossReference", {})
                    if cf_xrefs.get("database") == "ChEBI":
                        edge_props["cofactor_chebi_id"] = cf_xrefs.get("id")
            
            elif comment_type == "SUBCELLULAR LOCATION":
                locs = comment.get("subcellularLocations", [])
                if locs:
                    loc_info = locs[0].get("location", {})
                    edge_props["location_id"] = loc_info.get("id")
            
            # Store annotation node (just the label, no content)
            if node_label not in self.annotation_nodes:
                self.annotation_nodes[node_label] = {}
            
            # Store protein -> annotation relationship
            annotations.append((node_label, edge_props))
        
        if annotations:
            self.protein_annotations[protein_id] = annotations

    def _get_comment_text(self, comment: dict) -> str:
        """Extract text from various comment structures."""
        comment_type = comment.get("commentType", "")
        
        # Standard text field
        texts = comment.get("texts", [])
        if texts:
            return texts[0].get("value", "")
        
        # Subcellular location
        if comment_type == "SUBCELLULAR LOCATION":
            locs = comment.get("subcellularLocations", [])
            if locs:
                loc_values = []
                for loc_data in locs:
                    loc_name = loc_data.get("location", {}).get("value")
                    if loc_name:
                        loc_values.append(loc_name)
                return "; ".join(loc_values)
        
        # Catalytic activity
        if comment_type == "CATALYTIC ACTIVITY":
            reaction = comment.get("reaction", {})
            return reaction.get("name", "")
        
        # Cofactor
        if comment_type == "COFACTOR":
            cofactors = comment.get("cofactors", [])
            if cofactors:
                names = [c.get("name", "") for c in cofactors if c.get("name")]
                return "; ".join(names)
        
        return ""

    def _extract_features(self, protein_id: str, entry: dict):
        """Dynamically extract all feature types as sequence features."""
        features = entry.get("features", [])
        feature_data = []
        
        for idx, feature in enumerate(features):
            feature_type = feature.get("type", "UNKNOWN")
            
            # Generate feature ID
            feature_id = f"{protein_id}_{feature_type}_{idx}"
            
            # Determine node label from feature type
            node_label = feature_type.lower().replace(" ", "_").replace("-", "_") + "_feature"
            
            # Extract location with modifiers
            location = feature.get("location", {})
            start_info = location.get("start", {})
            end_info = location.get("end", {})
            start = start_info.get("value")
            end = end_info.get("value")
            
            # Edge properties - basic info
            edge_props = {
                "start_position": start,
                "end_position": end,
                "start_modifier": start_info.get("modifier"),
                "end_modifier": end_info.get("modifier"),
                "description": feature.get("description", ""),
                "feature_id": feature.get("featureId"),
            }
            
            # Extract evidences (PubMed, PDB, PROSITE references)
            evidences = feature.get("evidences", [])
            if evidences:
                evidence_codes = [e.get("evidenceCode") for e in evidences if e.get("evidenceCode")]
                pubmed_ids = [e.get("id") for e in evidences if e.get("source") == "PubMed" and e.get("id")]
                pdb_ids = [e.get("id") for e in evidences if e.get("source") == "PDB" and e.get("id")]
                prosite_ids = [e.get("id") for e in evidences if e.get("source") == "PROSITE-ProRule" and e.get("id")]
                
                if evidence_codes:
                    edge_props["evidence_codes"] = evidence_codes
                if pubmed_ids:
                    edge_props["pubmed_ids"] = pubmed_ids
                if pdb_ids:
                    edge_props["pdb_ids"] = pdb_ids
                if prosite_ids:
                    edge_props["prosite_ids"] = prosite_ids
            
            # Extract feature cross-references (ChEBI, etc.)
            xrefs = feature.get("featureCrossReferences", [])
            if xrefs:
                chebi_ids = [x.get("id") for x in xrefs if x.get("database") == "ChEBI" and x.get("id")]
                if chebi_ids:
                    edge_props["chebi_ids"] = chebi_ids
            
            # Handle alternative sequence (variants, mutagenesis)
            alt_seq = feature.get("alternativeSequence", {})
            if alt_seq:
                edge_props["original_sequence"] = alt_seq.get("originalSequence")
                alt_seqs = alt_seq.get("alternativeSequences", [])
                if alt_seqs:
                    edge_props["alternative_sequence"] = alt_seqs[0]
            
            # Handle ligand (binding sites)
            ligand = feature.get("ligand", {})
            if ligand:
                edge_props["ligand_name"] = ligand.get("name")
                edge_props["ligand_id"] = ligand.get("id")
            
            # Store feature node (just the label)
            if node_label not in self.feature_nodes:
                self.feature_nodes[node_label] = {}
            
            # Store protein -> feature relationship
            feature_data.append((node_label, edge_props))
        
        if feature_data:
            self.protein_features[protein_id] = feature_data

    def _extract_diseases(self, protein_id: str, entry: dict):
        """Extract disease associations from DISEASE comments."""
        comments = entry.get("comments", [])
        disease_data = []
        
        for comment in comments:
            if comment.get("commentType") != "DISEASE":
                continue
            
            disease = comment.get("disease", {})
            if not disease:
                continue
            
            disease_id = disease.get("diseaseAccession")
            if not disease_id:
                continue
            
            # Store disease node
            if disease_id not in self.disease_nodes:
                self.disease_nodes[disease_id] = {
                    "disease_id": disease.get("diseaseId"),
                    "acronym": disease.get("acronym"),
                    "description": disease.get("description"),
                }
                # Extract OMIM ID from cross-reference
                xref = disease.get("diseaseCrossReference", {})
                if xref.get("database") == "MIM":
                    self.disease_nodes[disease_id]["omim_id"] = xref.get("id")
            
            # Edge properties
            edge_props = {
                "note": comment.get("note", {}).get("texts", [{}])[0].get("value") if comment.get("note") else None
            }
            
            disease_data.append((disease_id, edge_props))
        
        if disease_data:
            self.protein_diseases[protein_id] = disease_data

    def _extract_interactions(self, protein_id: str, entry: dict):
        """Extract protein-protein interactions from INTERACTION comments."""
        comments = entry.get("comments", [])
        interactions = []
        
        for comment in comments:
            if comment.get("commentType") != "INTERACTION":
                continue
            
            for interaction in comment.get("interactions", []):
                interactant1 = interaction.get("interactantOne", {})
                interactant2 = interaction.get("interactantTwo", {})
                
                # Get target protein ID
                target_id = interactant2.get("uniProtKBAccession")
                if not target_id:
                    continue
                
                # Edge properties
                edge_props = {
                    "number_of_experiments": interaction.get("numberOfExperiments"),
                    "organism_differ": interaction.get("organismDiffer", False),
                    "intact_id_source": interactant1.get("intActId"),
                    "intact_id_target": interactant2.get("intActId"),
                }
                
                interactions.append((target_id, edge_props))
        
        if interactions:
            self.protein_interactions[protein_id] = interactions

    # ========================================
    # Extended Node Generators
    # ========================================

    def get_annotation_nodes(
        self,
        label_prefix: str = ""
    ) -> Generator[tuple[str, str, dict], None, None]:
        """
        Yield Functional Annotation nodes.
        """
        if UniprotNodeType.FUNCTIONAL_ANNOTATION not in self.node_types:
            return

        logger.info("Preparing Functional Annotation nodes...")
        
        base_props = {
            "source": self.data_source,
            "licence": self.data_licence,
            "version": self.data_version,
        }
        
        # Yield one node per annotation type
        for node_label in self.annotation_nodes.keys():
            node_id = node_label
            yield (node_id, node_label, base_props.copy())

    def get_feature_nodes(
        self,
        label_prefix: str = ""
    ) -> Generator[tuple[str, str, dict], None, None]:
        """
        Yield Sequence Feature nodes.
        """
        if UniprotNodeType.SEQUENCE_FEATURE not in self.node_types:
            return

        logger.info("Preparing Sequence Feature nodes...")
        
        base_props = {
            "source": self.data_source,
            "licence": self.data_licence,
            "version": self.data_version,
        }
        
        # Yield one node per feature type
        for node_label in self.feature_nodes.keys():
            node_id = node_label
            yield (node_id, node_label, base_props.copy())

    def get_disease_nodes(
        self,
        disease_label: str = "uniprot_disease"
    ) -> Generator[tuple[str, str, dict], None, None]:
        """
        Yield UniProt Disease nodes.
        """
        if UniprotNodeType.UNIPROT_DISEASE not in self.node_types:
            return

        logger.info("Preparing UniProt Disease nodes...")
        
        for disease_id, props in self.disease_nodes.items():
            node_id = self.add_prefix_to_id("uniprot.disease", disease_id)
            node_props = {
                "source": self.data_source,
                "licence": self.data_licence,
                "version": self.data_version,
                "disease_accession": disease_id,
                **props
            }
            yield (node_id, disease_label, node_props)

    # ========================================
    # Extended Edge Generators
    # ========================================

    def get_keyword_edges(
        self,
        edge_label: str = "protein_has_keyword"
    ) -> Generator[tuple[None, str, str, str, dict], None, None]:
        """
        Yield Protein -> Keyword edges with category and name properties.
        """
        if UniprotEdgeType.PROTEIN_TO_KEYWORD not in self.edge_types:
            return

        logger.info("Preparing Protein -> Keyword edges...")

        base_props = {
            "source": self.data_source,
            "licence": self.data_licence,
            "version": self.data_version,
        }

        for protein_id, keyword_data in self.protein_keywords.items():
            source_id = self.add_prefix_to_id("uniprot", protein_id)
            for kw in keyword_data:
                target_id = self.add_prefix_to_id("uniprot.keyword", kw["id"])
                props = base_props.copy()
                props["keyword_category"] = kw.get("category")
                props["keyword_name"] = kw.get("name")
                yield (None, source_id, target_id, edge_label, props)

    def get_annotation_edges(
        self,
        edge_label: str = "protein_has_annotation"
    ) -> Generator[tuple[None, str, str, str, dict], None, None]:
        """
        Yield Protein -> Functional Annotation edges with text properties.
        """
        if UniprotEdgeType.PROTEIN_TO_ANNOTATION not in self.edge_types:
            return

        logger.info("Preparing Protein -> Annotation edges...")
        
        base_props = {
            "source": self.data_source,
            "licence": self.data_licence,
            "version": self.data_version,
        }
        
        for protein_id, annotations in self.protein_annotations.items():
            source_id = self.add_prefix_to_id("uniprot", protein_id)
            for node_label, edge_props in annotations:
                target_id = node_label  # Node ID is the label itself
                props = base_props.copy()
                props.update(edge_props)
                yield (None, source_id, target_id, edge_label, props)

    def get_feature_edges(
        self,
        edge_label: str = "protein_has_feature"
    ) -> Generator[tuple[None, str, str, str, dict], None, None]:
        """
        Yield Protein -> Sequence Feature edges with position properties.
        """
        if UniprotEdgeType.PROTEIN_TO_FEATURE not in self.edge_types:
            return

        logger.info("Preparing Protein -> Feature edges...")
        
        base_props = {
            "source": self.data_source,
            "licence": self.data_licence,
            "version": self.data_version,
        }
        
        for protein_id, features in self.protein_features.items():
            source_id = self.add_prefix_to_id("uniprot", protein_id)
            for node_label, edge_props in features:
                target_id = node_label  # Node ID is the label itself
                props = base_props.copy()
                props.update(edge_props)
                yield (None, source_id, target_id, edge_label, props)

    def get_disease_edges(
        self,
        edge_label: str = "protein_associated_with_uniprot_disease"
    ) -> Generator[tuple[None, str, str, str, dict], None, None]:
        """
        Yield Protein -> Disease edges.
        """
        if UniprotEdgeType.PROTEIN_TO_DISEASE not in self.edge_types:
            return

        logger.info("Preparing Protein -> Disease edges...")
        
        base_props = {
            "source": self.data_source,
            "licence": self.data_licence,
            "version": self.data_version,
        }
        
        for protein_id, diseases in self.protein_diseases.items():
            source_id = self.add_prefix_to_id("uniprot", protein_id)
            for disease_id, edge_props in diseases:
                target_id = self.add_prefix_to_id("uniprot.disease", disease_id)
                props = base_props.copy()
                props.update(edge_props)
                yield (None, source_id, target_id, edge_label, props)

    def get_interaction_edges(
        self,
        edge_label: str = "uniprot_protein_interacts_with_protein"
    ) -> Generator[tuple[None, str, str, str, dict], None, None]:
        """
        Yield Protein -> Protein interaction edges from UniProt.
        """
        if UniprotEdgeType.PROTEIN_INTERACTS_PROTEIN not in self.edge_types:
            return

        logger.info("Preparing UniProt Protein Interaction edges...")
        
        base_props = {
            "source": self.data_source,
            "licence": self.data_licence,
            "version": self.data_version,
        }
        
        for protein_id, interactions in self.protein_interactions.items():
            source_id = self.add_prefix_to_id("uniprot", protein_id)
            for target_protein_id, edge_props in interactions:
                target_id = self.add_prefix_to_id("uniprot", target_protein_id)
                props = base_props.copy()
                props.update(edge_props)
                yield (None, source_id, target_id, edge_label, props)

    def get_all_extended_nodes(self) -> Generator[tuple[str, str, dict], None, None]:
        """Yield all extended nodes (annotations, features, diseases)."""
        yield from self.get_annotation_nodes()
        yield from self.get_feature_nodes()
        yield from self.get_disease_nodes()

    def get_all_extended_edges(self) -> Generator[tuple[None, str, str, str, dict], None, None]:
        """Yield all extended edges (keywords, annotations, features, diseases, interactions)."""
        yield from self.get_keyword_edges()
        yield from self.get_annotation_edges()
        yield from self.get_feature_edges()
        yield from self.get_disease_edges()
        yield from self.get_interaction_edges()

