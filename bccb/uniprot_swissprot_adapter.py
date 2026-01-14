"""
UniProt SwissProt Adapter - reads from local JSON file instead of API.
Produces the same output format as uniprot_adapter.py.
Uses ijson for streaming large JSON files.
"""

import ijson
from typing import Optional, Generator
from enum import Enum, EnumMeta, auto
from functools import lru_cache

from tqdm import tqdm
import logging

logger = logging.getLogger(__name__)

try:
    from bioregistry import normalize_curie
except ImportError:
    def normalize_curie(curie):
        return curie


class UniprotEnumMeta(EnumMeta):
    def __contains__(cls, item):
        return item in cls.__members__.keys()


class UniprotNodeType(Enum, metaclass=UniprotEnumMeta):
    PROTEIN = auto()
    GENE = auto()
    ORGANISM = auto()


class UniprotEdgeType(Enum, metaclass=UniprotEnumMeta):
    PROTEIN_TO_ORGANISM = auto()
    GENE_TO_PROTEIN = auto()


class UniprotSwissprot:
    """
    Adapter that reads UniProt SwissProt data from a local JSON file.
    Produces the same output format as the original Uniprot adapter.
    """

    def __init__(
        self,
        json_path: str = "/GenSIvePFS/users/data/UniProt/UniProtKB_SwissProt/uniprotkb_reviewed_true_2025_11_04.json",
        node_types: Optional[list[UniprotNodeType]] = None,
        edge_types: Optional[list[UniprotEdgeType]] = None,
        add_prefix: bool = True,
        test_mode: bool = False,
    ):
        self.json_path = json_path
        self.add_prefix = add_prefix
        self.test_mode = test_mode
        self.test_limit = 100

        self.data_source = "uniprot"
        self.data_version = "2025_11"
        self.data_licence = "CC BY 4.0"

        self.node_types = node_types or list(UniprotNodeType)
        self.edge_types = edge_types or list(UniprotEdgeType)

    def _stream_entries(self) -> Generator[dict, None, None]:
        """Stream entries from JSON file using ijson."""
        count = 0
        with open(self.json_path, 'rb') as f:
            for entry in ijson.items(f, 'results.item'):
                if self.test_mode and count >= self.test_limit:
                    break
                yield entry
                count += 1

    @lru_cache
    def _add_prefix(self, prefix: str, identifier: str) -> str:
        if self.add_prefix and identifier:
            return normalize_curie(f"{prefix}:{identifier}")
        return identifier

    def _get_gene_id(self, gene_data: dict) -> Optional[str]:
        """Extract gene ID with priority: geneName > orfNames > orderedLocusNames"""
        if gene_data.get("geneName"):
            return gene_data["geneName"].get("value")
        if gene_data.get("orfNames"):
            return gene_data["orfNames"][0].get("value")
        if gene_data.get("orderedLocusNames"):
            return gene_data["orderedLocusNames"][0].get("value")
        return None

    def _extract_protein_props(self, entry: dict) -> dict:
        """Extract protein properties from entry."""
        props = {
            "source": self.data_source,
            "licence": self.data_licence,
            "version": self.data_version,
        }

        seq = entry.get("sequence", {})
        props["length"] = seq.get("length")
        props["mass"] = seq.get("molWeight")
        props["sequence"] = seq.get("value")

        desc = entry.get("proteinDescription", {})
        rec_name = desc.get("recommendedName", {})
        full_name = rec_name.get("fullName", {})
        props["primary_protein_name"] = full_name.get("value")
        props["protein_names"] = full_name.get("value")

        ec_nums = rec_name.get("ecNumbers", [])
        if ec_nums:
            props["ec"] = [e.get("value") for e in ec_nums]
            if len(props["ec"]) == 1:
                props["ec"] = props["ec"][0]

        org = entry.get("organism", {})
        props["organism_name"] = org.get("scientificName")
        props["organism_id"] = org.get("taxonId")

        return props

    def _extract_gene_props(self, gene_data: dict) -> dict:
        """Extract gene properties."""
        props = {
            "source": self.data_source,
            "licence": self.data_licence,
            "version": self.data_version,
        }

        if gene_data.get("geneName"):
            props["gene_symbol"] = gene_data["geneName"].get("value")

        synonyms = gene_data.get("synonyms", [])
        if synonyms:
            props["gene_names"] = [s.get("value") for s in synonyms]

        orf_names = gene_data.get("orfNames", [])
        if orf_names:
            props["orf_names"] = [o.get("value") for o in orf_names]

        return props

    def get_nodes(
        self,
        protein_label: str = "protein",
        gene_label: str = "gene",
        organism_label: str = "organism",
    ) -> Generator[tuple[str, str, dict], None, None]:
        """Yield nodes from SwissProt data (streaming)."""

        logger.info(f"Generating nodes for types: {[t.name for t in self.node_types]}")

        seen_organisms = set()
        seen_genes = set()

        for entry in tqdm(self._stream_entries(), desc="Processing entries"):
            protein_id = entry.get("primaryAccession")
            if not protein_id:
                continue

            protein_id_prefixed = self._add_prefix("uniprot", protein_id)

            if UniprotNodeType.PROTEIN in self.node_types:
                props = self._extract_protein_props(entry)
                yield (protein_id_prefixed, protein_label, props)

            if UniprotNodeType.GENE in self.node_types:
                for gene_data in entry.get("genes", []):
                    gene_id = self._get_gene_id(gene_data)
                    if gene_id and gene_id not in seen_genes:
                        seen_genes.add(gene_id)
                        gene_id_prefixed = self._add_prefix("gene", gene_id)
                        gene_props = self._extract_gene_props(gene_data)
                        yield (gene_id_prefixed, gene_label, gene_props)

            if UniprotNodeType.ORGANISM in self.node_types:
                org = entry.get("organism", {})
                taxon_id = org.get("taxonId")
                if taxon_id and taxon_id not in seen_organisms:
                    seen_organisms.add(taxon_id)
                    org_id = self._add_prefix("ncbitaxon", str(taxon_id))
                    org_props = {
                        "organism_name": org.get("scientificName"),
                        "source": self.data_source,
                        "licence": self.data_licence,
                        "version": self.data_version,
                    }
                    yield (org_id, organism_label, org_props)

    def get_edges(
        self,
        gene_to_protein_label: str = "Gene_encodes_protein",
        protein_to_organism_label: str = "Protein_belongs_to_organism",
    ) -> Generator[tuple[None, str, str, str, dict], None, None]:
        """Yield edges from SwissProt data (streaming)."""

        logger.info(f"Generating edges for types: {[t.name for t in self.edge_types]}")

        props = {
            "source": self.data_source,
            "licence": self.data_licence,
            "version": self.data_version,
        }

        for entry in tqdm(self._stream_entries(), desc="Processing edges"):
            protein_id = entry.get("primaryAccession")
            if not protein_id:
                continue

            protein_id_prefixed = self._add_prefix("uniprot", protein_id)

            if UniprotEdgeType.GENE_TO_PROTEIN in self.edge_types:
                for gene_data in entry.get("genes", []):
                    gene_id = self._get_gene_id(gene_data)
                    if gene_id:
                        gene_id_prefixed = self._add_prefix("gene", gene_id)
                        yield (None, gene_id_prefixed, protein_id_prefixed, gene_to_protein_label, props)

            if UniprotEdgeType.PROTEIN_TO_ORGANISM in self.edge_types:
                org = entry.get("organism", {})
                taxon_id = org.get("taxonId")
                if taxon_id:
                    org_id = self._add_prefix("ncbitaxon", str(taxon_id))
                    yield (None, protein_id_prefixed, org_id, protein_to_organism_label, props)
