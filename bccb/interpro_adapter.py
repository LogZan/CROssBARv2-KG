from __future__ import annotations

import os
import json
import hashlib
import collections
import h5py
import pandas as pd
import numpy as np
import gc
import requests

from pypath.share import curl, settings
from pypath.inputs import interpro

from contextlib import ExitStack
from bioregistry import normalize_curie
from tqdm import tqdm

from time import time
from biocypher._logger import logger

from typing import Literal, Union, Optional, Generator
from pydantic import BaseModel, DirectoryPath, FilePath, validate_call

from enum import Enum, EnumMeta

from . import cache_config
cache_config.setup_pypath_cache()

logger.debug(f"Loading module {__name__}.")


class InterProEnumMeta(EnumMeta):
    def __contains__(cls, item):
        return item in cls.__members__.keys()


class InterProNodeField(Enum, metaclass=InterProEnumMeta):
    """
    Domain node fields in InterPro
    """

    # primary attributes
    PROTEIN_COUNT = "protein_count"
    NAME = "name"
    TYPE = "type"
    PARENT_LIST = "parent_list"
    CHILD_LIST = "child_list"

    # member list attributes
    PFAM = "PFAM"

    # external attributes
    EC = "EC"

    # structural attributes
    PDB = "PDB"

    # embedding
    DOM2VEC_EMBEDDING = "dom2vec_embedding"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None

    @classmethod
    def get_primary_attributes(cls):
        """
        Returns primary InterPro attributes
        """
        return [
            cls.PROTEIN_COUNT.value,
            cls.NAME.value,
            cls.TYPE.value,
            cls.PARENT_LIST.value,
            cls.CHILD_LIST.value,
        ]

    @classmethod
    def get_member_list_attributes(cls):
        """
        Returns external InterPro attributes
        """
        return [cls.PFAM.value]

    @classmethod
    def get_external_attributes(cls):
        """
        Returns external InterPro attributes
        """
        return [cls.EC.value]

    @classmethod
    def get_structural_attributes(cls):
        """
        Returns structural InterPro attributes
        """
        return [cls.PDB.value]


class InterProEdgeField(Enum, metaclass=InterProEnumMeta):
    """
    Domain edge fields in InterPro
    """

    START = "start"
    END = "end"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class InterProModel(BaseModel):
    page_size: int = 150
    organism: int | Literal["*"] | None = None
    add_prefix: bool = True
    node_fields: Union[list[InterProNodeField], None] = None
    edge_fields: Union[list[InterProEdgeField], None] = None
    test_mode: bool = False
    data_dir: str | None = None


class InterPro:
    """
    Class that downloads InterPro data using pypath and reformats it to be ready
    for import into a BioCypher database.
    """

    def __init__(
        self,
        page_size: Optional[int] = 150,
        organism: Optional[int | Literal["*"] | None] = None,
        add_prefix: Optional[bool] = True,
        node_fields: Optional[Union[list[InterProNodeField], None]] = None,
        edge_fields: Optional[Union[list[InterProEdgeField], None]] = None,
        test_mode: Optional[bool] = False,
        data_dir: Optional[str] = None,
    ):
        """
        Args:
            page_size: page size of downloaded annotation data
            organism: organism code in NCBI taxid format, e.g. "9606" for human. If it is None or "*", downloads all organism data.
            add_prefix: if True, add prefix to database identifiers
            node_fields: `InterProNodeField` fields to be used in the graph, if it is None, select all fields.
            edge_fields: `InterProEdgeField` fields to be used in the graph, if it is None, select all fields.
            test_mode: limits amount of data for testing
            data_dir: directory containing InterPro data files (protein2ipr.dat, entry.list)
        """

        model = InterProModel(
            page_size=page_size,
            organism=organism,
            add_prefix=add_prefix,
            node_fields=node_fields,
            edge_fields=edge_fields,
            test_mode=test_mode,
            data_dir=data_dir,
        ).model_dump()

        self.page_size = model["page_size"]
        self.organism = (
            None if model["organism"] in ("*", None) else model["organism"]
        )
        self.add_prefix = model["add_prefix"]
        self.test_mode = model["test_mode"]
        self.data_dir = model["data_dir"]

        # set node and edge fields
        self.set_node_and_edge_fields(
            node_fields=model["node_fields"], edge_fields=model["edge_fields"]
        )

        self.early_stopping = None
        if model["test_mode"]:
            self.early_stopping = 100
            
    @validate_call
    def download_interpro_data(self, cache: bool = False,
                               debug: bool = False,
                               retries: int = 3,
                               dom2vec_embedding_path: FilePath | None = None) -> None:
        """
        Wrapper function to download InterPro data using pypath; used to access
        settings.
        Args:
            cache: if True, it uses the cached version of the data, otherwise
            forces download.
            debug: if True, turns on debug mode in pypath.
            retries: number of retries in case of download error.
        """
        # Set adapter-specific cache directory
        cache_config.set_adapter_cache('interpro')
        
        # stack pypath context managers
        with ExitStack() as stack:

            stack.enter_context(settings.settings.context(curl_retries=retries))

            if debug:
                stack.enter_context(curl.debug_on())

            if not cache:
                stack.enter_context(curl.cache_off())

            self.download_domain_node_data(
                dom2vec_embedding_path=dom2vec_embedding_path,
                cache=cache,
            )
            self.download_domain_edge_data(cache=cache)

    @validate_call
    def download_domain_node_data(self, dom2vec_embedding_path: FilePath | None = None,
                                  cache: bool = False) -> None:
        """
        Downloads domain node data from InterPro API (reviewed proteins).

        Args:
            dom2vec_embedding_path: Path to dom2vec embedding file
            cache: If True, uses cached data
        """

        logger.info("Downloading InterPro domain data from InterPro API...")
        t0 = time()

        # Initialize data structures
        self.interpro_entries = {}  # entry_id -> entry metadata
        self.protein_domains = {}   # protein_id -> list of domains with locations
        # No xrefs available from the protein/entry API; keep empty dicts
        self.interpro_external_xrefs = {}
        self.interpro_structural_xrefs = {}

        # If cache=True and data_dir is provided, use local files
        if cache and self.data_dir:
            protein2ipr_path = os.path.join(self.data_dir, "protein2ipr.dat")
            entry_list_path = os.path.join(self.data_dir, "entry.list")
            if os.path.exists(protein2ipr_path) and os.path.exists(entry_list_path):
                logger.info(f"Using local InterPro files from {self.data_dir}")
                self._load_interpro_from_local_files(
                    protein2ipr_path=protein2ipr_path,
                    entry_list_path=entry_list_path,
                )

                # Load embeddings if requested
                if (
                    InterProNodeField.DOM2VEC_EMBEDDING.value in self.node_fields
                    and dom2vec_embedding_path is not None
                ):
                    self.retrieve_dom2vec_embeddings(
                        dom2vec_embedding_path=dom2vec_embedding_path
                    )

                t1 = time()
                logger.info(
                    f"InterPro domain data is retrieved in {round((t1-t0) / 60, 2)} mins"
                )
                return

        # Step 1: Get list of reviewed proteins
        logger.info("Fetching reviewed protein list...")
        protein_list = []
        page_size = 100
        protein_url = (
            "https://www.ebi.ac.uk/interpro/api/protein/reviewed/taxonomy/uniprot/"
            f"?page_size={page_size}"
        )

        # Set protein limit based on test mode
        max_proteins = 100 if self.test_mode else None

        cache_dir = cache_config.get_cache_dir("interpro")
        tax_label = (
            str(self.organism)
            if self.organism not in (None, "*")
            else "all"
        )

        with tqdm(desc="Fetching protein list", unit="page") as pbar:
            while protein_url and (max_proteins is None or len(protein_list) < max_proteins):
                try:
                    data = self._fetch_json(
                        protein_url,
                        cache=cache,
                        cache_dir=cache_dir,
                        tax_label=tax_label,
                    )

                    for protein in data.get('results', []):
                        protein_acc = protein.get('metadata', {}).get('accession')
                        if protein_acc:
                            protein_list.append(protein_acc.upper())
                            if max_proteins and len(protein_list) >= max_proteins:
                                break

                    protein_url = data.get('next')
                    pbar.update(1)
                    pbar.set_postfix({'proteins': len(protein_list)})

                except Exception as e:
                    logger.error(f"Error fetching protein list: {e}")
                    break

        logger.info(f"Fetched {len(protein_list)} proteins")

        # Step 2: For each protein, get its InterPro entries
        logger.info("Fetching domain mappings for proteins...")

        with tqdm(total=len(protein_list), desc="Fetching protein domains", unit="protein") as pbar:
            for protein_id in protein_list:
                try:
                    url = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/reviewed/{protein_id}/"
                    data = self._fetch_json(
                        url,
                        cache=cache,
                        cache_dir=cache_dir,
                        tax_label=tax_label,
                    )

                    # Process entries for this protein
                    for entry in data.get('results', []):
                        entry_metadata = entry.get('metadata', {})
                        entry_id = entry_metadata.get('accession')

                        if not entry_id:
                            continue

                        # Store entry metadata if not already stored
                        if entry_id not in self.interpro_entries:
                            self.interpro_entries[entry_id] = {
                                'name': entry_metadata.get('name'),
                                'type': entry_metadata.get('type'),
                                'member_databases': entry_metadata.get('member_databases', {})
                            }

                        # Get domain locations for this protein
                        for protein_data in entry.get('proteins', []):
                            if protein_data.get('accession', '').upper() == protein_id:
                                for location in protein_data.get('entry_protein_locations', []):
                                    for fragment in location.get('fragments', []):
                                        domain_info = {
                                            'entry_id': entry_id,
                                            'entry_name': entry_metadata.get('name'),
                                            'entry_type': entry_metadata.get('type'),
                                            'start': fragment.get('start'),
                                            'end': fragment.get('end')
                                        }

                                        if protein_id not in self.protein_domains:
                                            self.protein_domains[protein_id] = []
                                        self.protein_domains[protein_id].append(domain_info)

                    pbar.update(1)

                except Exception as e:
                    logger.warning(f"Error fetching domains for {protein_id}: {e}")
                    pbar.update(1)
                    continue

        logger.info(f"Downloaded {len(self.interpro_entries)} InterPro entries")
        logger.info(f"Mapped {len(self.protein_domains)} proteins to domains")

        # Load embeddings if requested
        if InterProNodeField.DOM2VEC_EMBEDDING.value in self.node_fields and dom2vec_embedding_path is not None:
            self.retrieve_dom2vec_embeddings(dom2vec_embedding_path=dom2vec_embedding_path)

        t1 = time()
        action = "retrieved" if cache else "downloaded"
        logger.info(
            f"InterPro domain data is {action} in {round((t1-t0) / 60, 2)} mins"
        )

    def _load_interpro_from_local_files(
        self,
        protein2ipr_path: str,
        entry_list_path: str,
    ) -> None:
        """Save file paths for streaming access - no data loaded into memory."""
        self.protein2ipr_path = protein2ipr_path
        self.entry_list_path = entry_list_path
        logger.info(f"InterPro file paths saved for streaming access: {protein2ipr_path}")

    def _fetch_json(self, url: str, cache: bool, cache_dir: str, tax_label: str | None = None) -> dict:
        """Fetch JSON with optional file cache under the InterPro cache dir."""
        cache_root = os.path.join(cache_dir, "interpro_api")
        os.makedirs(cache_root, exist_ok=True)
        cache_key = hashlib.md5(url.encode("utf-8")).hexdigest()
        label = tax_label if tax_label else "all"
        cache_path = os.path.join(cache_root, f"{cache_key}-{label}-interpro.json")

        if cache and os.path.exists(cache_path):
            with open(cache_path, "r", encoding="utf-8") as fh:
                return json.load(fh)

        response = requests.get(url, timeout=30)
        response.raise_for_status()
        data = response.json()

        if cache:
            with open(cache_path, "w", encoding="utf-8") as fh:
                json.dump(data, fh)

        return data

    def download_domain_edge_data(self, cache: bool = False) -> None:
        """
        Builds protein-domain annotations from InterPro API results.
        """

        logger.debug("Started downloading InterPro annotation data")
        t0 = time()

        if not hasattr(self, "protein_domains") or not self.protein_domains:
            logger.warning("No protein domain data found; cannot build edges.")
            self.interpro_annotations = {}
        else:
            InterproAnnotation = collections.namedtuple(
                "InterproAnnotation", ("interpro_id", "start", "end")
            )
            annotations = {}
            for protein_id, domains in self.protein_domains.items():
                ann_list = []
                for domain in domains:
                    ann_list.append(
                        InterproAnnotation(
                            interpro_id=domain.get("entry_id"),
                            start=domain.get("start"),
                            end=domain.get("end"),
                        )
                    )
                if ann_list:
                    annotations[protein_id] = ann_list
            self.interpro_annotations = annotations

        t1 = time()
        action = "retrieved" if cache else "downloaded"
        logger.info(
            f"InterPro annotation data is {action} in {round((t1-t0) / 60, 2)} mins"
        )

    def retrieve_dom2vec_embeddings(self, 
                                    dom2vec_embedding_path: FilePath | None = None) -> None:

        logger.info("Retrieving dom2vec domain embeddings.")

        self.interpro_id_to_dom2vec_embedding = {}
        with h5py.File(dom2vec_embedding_path, "r") as f:
            for interpro_id, embedding in f.items():
                self.interpro_id_to_dom2vec_embedding[interpro_id] = np.array(embedding).astype(np.float16)
                    
    @validate_call
    def get_interpro_nodes(self, node_label: str = "domain"):
        """
        Stream InterPro domain nodes for BioCypher using generator
        Args:
            node_label : label of interpro nodes
        """
        # Check if we have file paths (local file mode)
        if hasattr(self, 'entry_list_path') and self.entry_list_path:
            # Use streaming from local files
            logger.info("Streaming InterPro nodes from local file")
            counter = 0

            file_size = os.path.getsize(self.entry_list_path)

            with open(self.entry_list_path, "r", encoding="utf-8") as fh:
                with tqdm(total=file_size, unit='B', unit_scale=True, unit_divisor=1024,
                         desc="Processing InterPro nodes") as pbar:
                    bytes_since_update = 0
                    header = True
                    for line in fh:
                        line_bytes = len(line.encode('utf-8'))
                        bytes_since_update += line_bytes

                        line = line.rstrip("\n")
                        if not line:
                            continue
                        if header:
                            header = False
                            continue
                        parts = line.split("\t")
                        if len(parts) < 3:
                            continue

                        entry_id, entry_type, entry_name = parts[0], parts[1], parts[2]
                        domain_id = self.add_prefix_to_id("interpro", entry_id)

                        props = {}
                        if InterProNodeField.NAME.value in self.node_fields:
                            props["name"] = entry_name.replace("'", "^").replace("|", ",") if entry_name else None
                        if InterProNodeField.TYPE.value in self.node_fields:
                            props["type"] = entry_type

                        yield (domain_id, node_label, props)

                        counter += 1

                        # Update progress bar every 1000 lines
                        if counter % 1000 == 0:
                            pbar.update(bytes_since_update)
                            bytes_since_update = 0

                        if self.early_stopping and counter >= self.early_stopping:
                            break

                    # Final update for remaining bytes
                    if bytes_since_update > 0:
                        pbar.update(bytes_since_update)

            logger.info(f"Streamed {counter} InterPro nodes")
        else:
            # Fallback to API mode (load data into memory)
            if not getattr(self, "interpro_entries", None):
                self.download_domain_node_data()

            logger.info("Generating InterPro nodes from API data")
            counter = 0

            for entry_id, meta in getattr(self, "interpro_entries", {}).items():
                domain_id = self.add_prefix_to_id("interpro", entry_id)

                props = {}
                if InterProNodeField.NAME.value in self.node_fields:
                    props["name"] = meta.get("name", "").replace("'", "^").replace("|", ",") if meta.get("name") else None
                if InterProNodeField.TYPE.value in self.node_fields:
                    props["type"] = meta.get("type")

                yield (domain_id, node_label, props)

                counter += 1
                if self.early_stopping and counter >= self.early_stopping:
                    break

            logger.info(f"Generated {counter} InterPro nodes")

    def get_interpro_edges(
        self, edge_label: str = "protein_has_domain"
    ):
        """
        Stream Protein-Domain edges for BioCypher using generator
        Args:
            edge_label: label of protein-domain edge
        """
        # Check if we have file paths (local file mode)
        if hasattr(self, 'protein2ipr_path') and self.protein2ipr_path:
            # Use streaming from local files
            logger.info("Streaming InterPro edges from local file")
            counter = 0
            max_proteins = 100 if self.test_mode else None
            seen_proteins = set()

            file_size = os.path.getsize(self.protein2ipr_path)

            # Pre-compute prefix settings for hot loop
            use_prefix = self.add_prefix
            check_start = "start" in self.edge_fields
            check_end = "end" in self.edge_fields

            with open(self.protein2ipr_path, "rt", encoding="utf-8") as fh:
                with tqdm(total=file_size, unit='B', unit_scale=True, unit_divisor=1024,
                         desc="Processing protein-domain edges") as pbar:
                    bytes_since_update = 0
                    for line in fh:
                        line_bytes = len(line.encode('utf-8'))
                        bytes_since_update += line_bytes

                        line = line.rstrip("\n")
                        if not line:
                            continue
                        parts = line.split("\t")
                        if len(parts) < 6:
                            continue

                        protein_id, entry_id, entry_name, _, start, end = parts[:6]
                        protein_id = protein_id.upper()

                        # Early stopping for test mode
                        if max_proteins is not None:
                            if protein_id not in seen_proteins:
                                if len(seen_proteins) >= max_proteins:
                                    break
                                seen_proteins.add(protein_id)

                        # Inline prefix logic (avoid function call overhead)
                        if use_prefix:
                            interpro_id = f"interpro:{entry_id}"
                            uniprot_id = f"uniprot:{protein_id}"
                        else:
                            interpro_id = entry_id
                            uniprot_id = protein_id

                        props = {}
                        if check_start:
                            props["start"] = int(start) if start.isdigit() else start
                        if check_end:
                            props["end"] = int(end) if end.isdigit() else end

                        yield (None, uniprot_id, interpro_id, edge_label, props)

                        counter += 1

                        # Update progress bar every 100k lines (reduced overhead)
                        if counter % 100000 == 0:
                            pbar.update(bytes_since_update)
                            bytes_since_update = 0

                        if self.early_stopping and counter >= self.early_stopping:
                            break

                    # Final update for remaining bytes
                    if bytes_since_update > 0:
                        pbar.update(bytes_since_update)

            logger.info(f"Streamed {counter} InterPro edges")
        else:
            # Fallback to API mode (load data into memory)
            if not hasattr(self, "interpro_annotations"):
                self.download_domain_edge_data()

            logger.info("Generating InterPro edges from API data")
            counter = 0

            for protein_id, annotations in getattr(self, "interpro_annotations", {}).items():
                for annotation in annotations:
                    interpro_id = self.add_prefix_to_id("interpro", annotation.interpro_id)
                    uniprot_id = self.add_prefix_to_id("uniprot", protein_id)

                    props = {}
                    if hasattr(annotation, 'start') and "start" in self.edge_fields:
                        props["start"] = annotation.start
                    if hasattr(annotation, 'end') and "end" in self.edge_fields:
                        props["end"] = annotation.end

                    yield (None, uniprot_id, interpro_id, edge_label, props)

                    counter += 1
                    if self.early_stopping and counter >= self.early_stopping:
                        break

                if self.early_stopping and counter >= self.early_stopping:
                    break

            logger.info(f"Generated {counter} InterPro edges")

    @validate_call
    def check_length(self, element: str | list | int | float) -> str | list | int | float:
        """
        If the type of given entry is a list and has just one element returns this one element
        """
        if isinstance(element, list) and len(element) == 1:
            return element[0]
        else:
            return element

    @validate_call
    def add_prefix_to_id(
        self, prefix: str = None, identifier: str = None, sep: str = ":"
    ) -> str:
        """
        Adds prefix to ids
        """
        if self.add_prefix:
            return prefix + sep + identifier

        return identifier

    def set_node_and_edge_fields(self, node_fields, edge_fields) -> None:
        """
        Sets Interpro node and edge fields
        """

        if node_fields:
            self.node_fields = [field.value for field in node_fields]
        else:
            self.node_fields = [field.value for field in InterProNodeField]

        if edge_fields:
            self.edge_fields = [field.value for field in edge_fields]
        else:
            self.edge_fields = [field.value for field in InterProEdgeField]

    def export_as_csv(self, path: DirectoryPath | None = None):
        import csv
        
        if path:
            node_full_path = os.path.join(path, "Domain.csv")
            edge_full_path = os.path.join(path, "Protein_has_domain.csv")
        else:
            node_full_path = os.path.join(os.getcwd(), "Domain.csv")
            edge_full_path = os.path.join(os.getcwd(), "Protein_has_domain.csv")

        # write nodes (small dataset, can use pandas)
        nodes = self.get_interpro_nodes()
        node_df_list = []
        for n in nodes:
            props = {"id": n[0]}
            props |= n[2]
            node_df_list.append(props)

        nodes_df = pd.DataFrame.from_records(node_df_list)
        nodes_df.to_csv(node_full_path, index=False)
        logger.info(f"Domain node data is written: {node_full_path}")

        # write edges (large dataset, stream to CSV to avoid OOM)
        edges = self.get_interpro_edges()
        edge_count = 0
        
        with open(edge_full_path, 'w', newline='', encoding='utf-8') as f:
            writer = None
            for e in edges:
                props = {"source_id": e[1], "target_id": e[2]}
                props |= e[4]
                
                if writer is None:
                    writer = csv.DictWriter(f, fieldnames=list(props.keys()))
                    writer.writeheader()
                
                writer.writerow(props)
                edge_count += 1
        
        logger.info(f"Domain edge data is written: {edge_full_path} ({edge_count:,} edges)")
