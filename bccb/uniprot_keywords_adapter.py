"""
UniProt Keywords Adapter - reads from local keywords JSON file.
Produces Keyword nodes and Keyword hierarchy/GO mapping edges.
"""

from time import time
from typing import Optional, Union
from collections.abc import Generator
from enum import Enum, EnumMeta, auto
import json
import os

from tqdm import tqdm
from biocypher._logger import logger
from bioregistry import normalize_curie

from pydantic import BaseModel, validate_call

logger.debug(f"Loading module {__name__}.")


class KeywordEnumMeta(EnumMeta):
    def __contains__(cls, item):
        return item in cls.__members__.keys()


class KeywordNodeType(Enum, metaclass=KeywordEnumMeta):
    """
    Node types of the UniProt Keywords represented in this adapter.
    """
    KEYWORD = auto()


class KeywordNodeField(Enum, metaclass=KeywordEnumMeta):
    """
    Fields of Keyword nodes represented in this adapter.
    """
    NAME = "name"
    DEFINITION = "definition"
    CATEGORY = "category"
    SYNONYMS = "synonyms"
    LINKS = "links"
    REVIEWED_PROTEIN_COUNT = "reviewed_protein_count"
    UNREVIEWED_PROTEIN_COUNT = "unreviewed_protein_count"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class KeywordEdgeType(Enum, metaclass=KeywordEnumMeta):
    """
    Edge types of the UniProt Keywords represented in this adapter.
    """
    KEYWORD_HIERARCHY = auto()  # Keyword -> Keyword (IS-A / CHILD_OF)
    KEYWORD_TO_GO = auto()      # Keyword -> GO Term


class KeywordEdgeField(Enum, metaclass=KeywordEnumMeta):
    """
    Fields of Keyword edges represented in this adapter.
    """
    GO_NAME = "go_name"


class UniprotKeywordsModel(BaseModel):
    json_path: str = "/GenSIvePFS/users/data/UniProt/keywords_all_2025_12_25.json"
    node_types: Optional[Union[list[KeywordNodeType], None]] = None
    node_fields: Optional[Union[list[KeywordNodeField], None]] = None
    edge_types: Optional[Union[list[KeywordEdgeType], None]] = None
    add_prefix: bool = True
    test_mode: bool = False


class UniprotKeywords:
    """
    Adapter that reads UniProt Keywords from a local JSON file.
    
    The keywords vocabulary provides:
    - 1201 keyword entries with hierarchical relationships
    - Mappings to Gene Ontology terms
    - Category classifications (11 categories)
    
    Args:
        json_path: path to the keywords JSON file.
        node_types: `KeywordNodeType` to include. Default: all.
        node_fields: `KeywordNodeField` to include. Default: all.
        edge_types: `KeywordEdgeType` to include. Default: all.
        add_prefix: if True, add prefix to identifiers.
        test_mode: if True, limits output data for testing.
    """

    def __init__(
        self,
        json_path: str = "/GenSIvePFS/users/data/UniProt/keywords_all_2025_12_25.json",
        node_types: Optional[Union[list[KeywordNodeType], None]] = None,
        node_fields: Optional[Union[list[KeywordNodeField], None]] = None,
        edge_types: Optional[Union[list[KeywordEdgeType], None]] = None,
        add_prefix: Optional[bool] = True,
        test_mode: Optional[bool] = False,
    ):
        model = UniprotKeywordsModel(
            json_path=json_path,
            node_types=node_types,
            node_fields=node_fields,
            edge_types=edge_types,
            add_prefix=add_prefix,
            test_mode=test_mode,
        ).model_dump()

        # params
        self.json_path = model["json_path"]
        self.add_prefix = model["add_prefix"]
        self.test_mode = model["test_mode"]
        self.test_limit = 100

        # provenance
        self.data_source = "uniprot"
        self.data_version = "2025_12"
        self.data_licence = "CC BY 4.0"

        # Set node and edge types/fields
        self.node_types = model["node_types"] or list(KeywordNodeType)
        self.node_fields = [f.value for f in (model["node_fields"] or list(KeywordNodeField))]
        self.edge_types = model["edge_types"] or list(KeywordEdgeType)

        # Data storage
        self.keywords_data = []
        self.keyword_ids = set()

    def download_data(self, cache: bool = False):
        """
        Load keywords data from local JSON file.
        """
        self._load_json_data()

    def _load_json_data(self):
        """Load and parse JSON data from file."""
        logger.info(f"Loading Keywords data from {self.json_path}...")
        t0 = time()

        if not os.path.exists(self.json_path):
            raise FileNotFoundError(f"Keywords file not found: {self.json_path}")

        with open(self.json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)

        results = data.get("results", [])
        
        count = 0
        for entry in results:
            if self.test_mode and count >= self.test_limit:
                break
            self.keywords_data.append(entry)
            kw_id = entry.get("keyword", {}).get("id")
            if kw_id:
                self.keyword_ids.add(kw_id)
            count += 1

        t1 = time()
        msg = f"Loaded {len(self.keywords_data)} keywords in {round((t1-t0), 2)} seconds."
        logger.info(msg)

    @validate_call
    def get_nodes(
        self,
        keyword_label: str = "uniprot_keyword"
    ) -> Generator[tuple[str, str, dict], None, None]:
        """
        Yield Keyword nodes.
        
        Returns:
            Generator of (node_id, node_label, properties) tuples.
        """
        if KeywordNodeType.KEYWORD not in self.node_types:
            return

        logger.info("Preparing Keyword nodes...")

        for entry in tqdm(self.keywords_data, desc="Processing keywords"):
            kw_info = entry.get("keyword", {})
            kw_id = kw_info.get("id")
            kw_name = kw_info.get("name")

            if not kw_id:
                continue

            node_id = self.add_prefix_to_id("uniprot.keyword", kw_id)
            
            props = {
                "source": self.data_source,
                "licence": self.data_licence,
                "version": self.data_version,
            }

            # Add properties based on node_fields
            if KeywordNodeField.NAME.value in self.node_fields:
                props["name"] = kw_name

            if KeywordNodeField.DEFINITION.value in self.node_fields:
                props["definition"] = entry.get("definition")

            if KeywordNodeField.CATEGORY.value in self.node_fields:
                category = entry.get("category", {})
                props["category"] = category.get("name")

            if KeywordNodeField.SYNONYMS.value in self.node_fields:
                props["synonyms"] = entry.get("synonyms", [])

            if KeywordNodeField.LINKS.value in self.node_fields:
                props["links"] = entry.get("links", [])

            if KeywordNodeField.REVIEWED_PROTEIN_COUNT.value in self.node_fields:
                stats = entry.get("statistics", {})
                props["reviewed_protein_count"] = stats.get("reviewedProteinCount")

            if KeywordNodeField.UNREVIEWED_PROTEIN_COUNT.value in self.node_fields:
                stats = entry.get("statistics", {})
                props["unreviewed_protein_count"] = stats.get("unreviewedProteinCount")

            yield (node_id, keyword_label, props)

    @validate_call
    def get_edges(
        self,
        hierarchy_label: str = "keyword_is_a_keyword",
        go_mapping_labels: dict = None
    ) -> Generator[tuple[None, str, str, str, dict], None, None]:
        """
        Yield Keyword edges (hierarchy and GO mappings).
        
        Returns:
            Generator of (None, source_id, target_id, edge_label, properties) tuples.
        """
        if go_mapping_labels is None:
            go_mapping_labels = {
                "biological process": "keyword_maps_to_biological_process",
                "molecular function": "keyword_maps_to_molecular_function",
                "cellular component": "keyword_maps_to_cellular_component",
            }

        logger.info("Preparing Keyword edges...")

        base_props = {
            "source": self.data_source,
            "licence": self.data_licence,
            "version": self.data_version,
        }

        for entry in tqdm(self.keywords_data, desc="Processing keyword edges"):
            kw_info = entry.get("keyword", {})
            kw_id = kw_info.get("id")

            if not kw_id:
                continue

            source_id = self.add_prefix_to_id("uniprot.keyword", kw_id)

            # Hierarchy edges (Keyword -> Parent Keyword)
            if KeywordEdgeType.KEYWORD_HIERARCHY in self.edge_types:
                for parent in entry.get("parents", []):
                    parent_kw = parent.get("keyword", {})
                    parent_id = parent_kw.get("id")
                    if parent_id:
                        target_id = self.add_prefix_to_id("uniprot.keyword", parent_id)
                        yield (None, source_id, target_id, hierarchy_label, base_props.copy())

            # # GO mapping edges (Keyword -> GO Term)
            # if KeywordEdgeType.KEYWORD_TO_GO in self.edge_types:
            #     for go_entry in entry.get("geneOntologies", []):
            #         go_id = go_entry.get("goId")
            #         go_name = go_entry.get("name")
            #         if go_id:
            #             target_id = self.add_prefix_to_id("go", go_id)
                        
            #             # Determine edge label based on GO category (would need lookup)
            #             # For now, use a generic label or derive from category
            #             category = entry.get("category", {}).get("name", "").lower()
            #             edge_label = go_mapping_labels.get(category, "keyword_maps_to_go")
                        
            #             edge_props = base_props.copy()
            #             edge_props["go_name"] = go_name
                        
            #             yield (None, source_id, target_id, edge_label, edge_props)

    def add_prefix_to_id(
        self, prefix: str = None, identifier: str = None, sep: str = ":"
    ) -> str:
        """Adds prefix to database id."""
        if self.add_prefix and identifier:
            return normalize_curie(prefix + sep + identifier)
        return identifier
