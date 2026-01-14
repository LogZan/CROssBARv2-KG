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

from tqdm import tqdm
from biocypher._logger import logger
from bioregistry import normalize_curie

from pydantic import BaseModel, DirectoryPath, FilePath, validate_call

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

    # not from uniprot REST
    # we provide these by mapping ENSTs via pypath
    ENSEMBL_GENE_IDS = "ensembl_gene_ids"

    # not from uniprot REST
    # we provide these by downloading the ProtT5 embeddings from uniprot
    PROTT5_EMBEDDING = "prott5_embedding"

    # not from uniprot REST
    # we provide these by getting embeddings from ESM2 650M model
    ESM2_EMBEDDING = "esm2_embedding"

    # not from uniprot REST
    # we provide these by getting embeddings from Nucletide tranformer 2 model
    NT_EMBEDDING = "nt_embedding"

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
            cls.PROTT5_EMBEDDING.value,
            cls.ESM2_EMBEDDING.value,
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
        ]

    @classmethod
    def get_organism_properties(cls) -> list:
        return [cls.ORGANISM.value]

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
        self._load_json_data()

        # Load embeddings if requested
        if UniprotNodeField.PROTT5_EMBEDDING.value in self.node_fields:
            self.data[UniprotNodeField.PROTT5_EMBEDDING.value] = {}
            if prott5_embedding_output_path:
                self.download_prott5_embeddings(
                    prott5_embedding_output_path=prott5_embedding_output_path
                )

        if UniprotNodeField.ESM2_EMBEDDING.value in self.node_fields:
            self.data[UniprotNodeField.ESM2_EMBEDDING.value] = {}
            if esm2_embedding_path:
                self.retrieve_esm2_embeddings(esm2_embedding_path)

        if UniprotNodeField.NT_EMBEDDING.value in self.node_fields:
            self.data[UniprotNodeField.NT_EMBEDDING.value] = {}
            if nucleotide_transformer_embedding_path:
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

        # Extract organism info
        org = entry.get("organism", {})
        if UniprotNodeField.ORGANISM.value in self.node_fields:
            self.data[UniprotNodeField.ORGANISM.value][protein_id] = org.get("scientificName")

        if UniprotNodeField.ORGANISM_ID.value in self.node_fields:
            self.data[UniprotNodeField.ORGANISM_ID.value][protein_id] = org.get("taxonId")

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

            if UniprotNodeField.PROTEIN_GENE_NAMES.value in self.node_fields and gene_names:
                self.data[UniprotNodeField.PROTEIN_GENE_NAMES.value][protein_id] = gene_names if len(gene_names) > 1 else gene_names[0]

            if UniprotNodeField.PRIMARY_GENE_NAME.value in self.node_fields and primary_gene:
                self.data[UniprotNodeField.PRIMARY_GENE_NAME.value][protein_id] = primary_gene

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
        except ImportError:
            # pypath not available, skip ENSG mapping
            pass

        ensg_ids = list(ensg_ids)
        if len(ensg_ids) == 1:
            ensg_ids = ensg_ids[0]
        return ensg_ids if ensg_ids else None

    @validate_call
    def download_prott5_embeddings(
        self,
        prott5_embedding_output_path: FilePath | None = None
    ):
        """
        Load ProtT5 embeddings from h5 file.
        """
        if not prott5_embedding_output_path or not os.path.isfile(prott5_embedding_output_path):
            logger.warning("ProtT5 embedding file not found.")
            return

        logger.info("Loading ProtT5 embeddings...")

        df_list = []
        with h5py.File(prott5_embedding_output_path, "r") as file:
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
        esm2_embedding_path: FilePath | None = None
    ) -> None:
        """Load ESM2 embeddings from h5 file."""
        if not esm2_embedding_path or not os.path.isfile(esm2_embedding_path):
            logger.warning("ESM2 embedding file not found.")
            return

        logger.info("Loading ESM2 embeddings...")

        df_list = []
        with h5py.File(esm2_embedding_path, "r") as file:
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
        nucleotide_transformer_embedding_path: FilePath | None = None
    ) -> None:
        """Load Nucleotide Transformer embeddings from h5 file."""
        if not nucleotide_transformer_embedding_path or not os.path.isfile(nucleotide_transformer_embedding_path):
            logger.warning("Nucleotide Transformer embedding file not found.")
            return

        logger.info("Loading Nucleotide Transformer embeddings...")

        self.entrez_id_to_nucleotide_transformer_embedding = {}
        with h5py.File(nucleotide_transformer_embedding_path, "r") as file:
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
