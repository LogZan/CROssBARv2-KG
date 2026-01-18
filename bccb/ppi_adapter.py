import os
import pandas as pd
import numpy as np
import time
import collections
from time import time

from pypath.inputs import intact
from pypath.inputs import string
from pypath.inputs import biogrid
from pypath.share import curl, settings
from pypath.inputs import uniprot

from tqdm import tqdm  # progress bar

from biocypher._logger import logger
from pypath.resources import urls
from contextlib import ExitStack
from pydantic import BaseModel, DirectoryPath, validate_call
from typing import Literal, Union, Optional

from bioregistry import normalize_curie

from enum import Enum, EnumMeta

from . import cache_config
cache_config.setup_pypath_cache()


class PPIEnumMeta(EnumMeta):
    def __contains__(cls, item):
        return item in cls.__members__.keys()


class IntactEdgeField(Enum, metaclass=PPIEnumMeta):
    SOURCE = "source"
    PUBMED_IDS = "pubmeds"
    INTACT_SCORE = "mi_score"
    METHODS = "methods"
    INTERACTION_TYPES = "interaction_types"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class BiogridEdgeField(Enum, metaclass=PPIEnumMeta):
    SOURCE = "source"
    PUBMED_IDS = "pmid"
    EXPERIMENTAL_SYSTEM = "experimental_system"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class StringEdgeField(Enum, metaclass=PPIEnumMeta):
    SOURCE = "source"
    COMBINED_SCORE = "combined_score"
    PHYSICAL_COMBINED_SCORE = "physical_combined_score"

    @classmethod
    def _missing_(cls, value: str):
        value = value.lower()
        for member in cls.__members__.values():
            if member.value.lower() == value:
                return member
        return None


class PPIModel(BaseModel):
    output_dir: Optional[Union[DirectoryPath, None]] = None
    export_csv: Optional[bool] = False
    organism: Optional[Union[int, Literal["*"], None]] = None
    intact_fields: Optional[Union[list[IntactEdgeField], None]] = None
    biogrid_fields: Optional[Union[list[BiogridEdgeField], None]] = None
    string_fields: Optional[Union[list[StringEdgeField], None]] = None
    add_prefix: Optional[bool] = True
    test_mode: Optional[bool] = False


class PPI:
    def __init__(
        self,
        output_dir: Optional[Union[DirectoryPath, None]] = None,
        export_csv: Optional[bool] = False,
        organism: Optional[Union[int, Literal["*"], None]] = None,
        intact_fields: Optional[Union[list[IntactEdgeField], None]] = None,
        biogrid_fields: Optional[Union[list[BiogridEdgeField], None]] = None,
        string_fields: Optional[Union[list[StringEdgeField], None]] = None,
        add_prefix: Optional[bool] = True,
        test_mode: Optional[bool] = False,
    ):
        """
        Downloads and processes PPI data

            Args:
                output_dir: Directory to save the output csv file
                export_csv: Flag for whether or not create csvs of outputs of databases
                cache: if True, it uses the cached version of the data, otherwise
                forces download.
                debug: if True, turns on debug mode in pypath.
                retries: number of retries in case of download error.
                organism: organism code in NCBI taxid format, e.g. "9606" for human. If it is None or "*", downloads all organism data.
                intact_fields: `IntactEdgeField` fields to be used in the graph, if it is None, select all fields.
                biogrid_fields: `BiogridEdgeField` fields to be used in the graph, if it is None, select all fields.
                string_fields: `StringEdgeField` fields to be used in the graph, if it is None, select all fields.
                add_prefix: if True, add prefix to uniprot ids.
                test_mode: if True, take small sample from data for testing.

        """
        model = PPIModel(
            output_dir=output_dir,
            export_csv=export_csv,
            organism=organism,
            intact_fields=intact_fields,
            biogrid_fields=biogrid_fields,
            string_fields=string_fields,
            add_prefix=add_prefix,
            test_mode=test_mode,
        ).model_dump()

        self.export_csv = model["export_csv"]
        self.organism = (
            None if model["organism"] in ("*", None) else model["organism"]
        )
        self.intact_fields = model["intact_fields"]
        self.biogrid_fields = model["biogrid_fields"]
        self.string_fields = model["string_fields"]
        self.add_prefix = model["add_prefix"]
        self.test_mode = model["test_mode"]

        self.swissprots = set(uniprot._all_uniprots("*", True))

        self.check_status_and_properties: dict[str, dict] = {
            "intact": {
                "downloaded": False,
                "processed": False,
                "properties_dict": None,
                "dataframe": None,
            },
            "biogrid": {
                "downloaded": False,
                "processed": False,
                "properties_dict": None,
                "dataframe": None,
            },
            "string": {
                "downloaded": False,
                "processed": False,
                "properties_dict": None,
                "dataframe": None,
            },
        }

        if model["export_csv"]:
            self.output_dir = model["output_dir"]

    @validate_call
    def download_ppi_data(self, cache: bool = False, debug: bool = False,
                          retries: int = 6) -> None:
        """
        Wrapper function to download PPI data using pypath; used to access
        settings.
        Args:
            cache: if True, it uses the cached version of the data, otherwise
            forces download.
            debug: if True, turns on debug mode in pypath.
            retries: number of retries in case of download error.
        """
        
        with ExitStack() as stack:
            stack.enter_context(settings.settings.context(retries=retries))

            if debug:
                stack.enter_context(curl.debug_on())
            if not cache:
                stack.enter_context(curl.cache_off())

            self.download_intact_data()
            self.download_biogrid_data()
            self.download_string_data()

    def process_ppi_data(self) -> None:
        self.intact_process()
        self.biogrid_process()
        self.string_process()
        
    def download_intact_data(self) -> None:
        """
        Wrapper function to download IntAct data using pypath; used to access
        settings.

        To do: Make arguments of intact.intact_interactions selectable for user.
        """

        logger.debug("Started downloading IntAct data")
        logger.info(
            f"This is the link of IntAct data we downloaded:{urls.urls['intact']['mitab']}. Please check if it is up to date"
        )
        t0 = time()

        self.intact_ints = intact.intact_interactions(
            miscore=0,
            organism=self.organism,
            complex_expansion=True,
            only_proteins=True,
        )

        t1 = time()
        logger.info(
            f"IntAct data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

        if self.test_mode:
            self.intact_ints = self.intact_ints[:100]

        self.check_status_and_properties["intact"]["downloaded"] = True

    @validate_call
    def intact_process(
        self, rename_selected_fields: dict[str, str] = None
    ) -> None:
        """
        Processor function for IntAct data. It drops duplicate and reciprocal duplicate protein pairs and collects pubmed ids of duplicated pairs. Also, it filters
        protein pairs found in swissprot.

        Args:
            rename_selected_fields : List of new field names for selected fields. If not defined, default field names will be used.
        """
        if self.intact_fields is None:
            selected_fields = [field.value for field in IntactEdgeField]
        else:
            selected_fields = [field.value for field in self.intact_fields]

        # if rename_selected_fields is not defined create column names from this dictionary
        default_field_names = {
            "source": "source",
            "pubmeds": "pubmed_ids",
            "mi_score": "intact_score",
            "methods": "method",
            "interaction_types": "interaction_type",
        }

        self.intact_field_new_names = {}

        if rename_selected_fields:
            if len(selected_fields) != len(rename_selected_fields):
                raise Exception(
                    "Length of selected_fields variable should be equal to length of rename_selected_fields variable"
                )

            for field_old_name, field_new_name in list(
                zip(selected_fields, rename_selected_fields)
            ):
                self.intact_field_new_names[field_old_name] = field_new_name

            self.intact_field_new_names["id_a"] = "uniprot_a"
            self.intact_field_new_names["id_b"] = "uniprot_b"

        else:
            for field_old_name in selected_fields:
                self.intact_field_new_names[field_old_name] = (
                    default_field_names[field_old_name]
                )

            self.intact_field_new_names["id_a"] = "uniprot_a"
            self.intact_field_new_names["id_b"] = "uniprot_b"

        logger.debug("Started processing IntAct data")
        t1 = time()

        # create dataframe
        intact_df = pd.DataFrame.from_records(
            self.intact_ints, columns=self.intact_ints[0]._fields
        )

        # turn list columns to string
        for list_column in ["pubmeds", "methods", "interaction_types"]:
            intact_df[list_column] = [
                ";".join(map(str, l)) for l in intact_df[list_column]
            ]

        intact_df.fillna(value=np.nan, inplace=True)

        # add source database info
        intact_df["source"] = "IntAct"

        # filter selected fields
        intact_df = intact_df[list(self.intact_field_new_names.keys())]

        # rename columns
        intact_df.rename(columns = self.intact_field_new_names, inplace=True)

        # drop rows if uniprot_a or uniprot_b is not a swiss-prot protein
        intact_df = intact_df[
            (intact_df["uniprot_a"].isin(self.swissprots))
            & (intact_df["uniprot_b"].isin(self.swissprots))
        ]
        intact_df.reset_index(drop=True, inplace=True)

        if "pubmeds" in self.intact_field_new_names.keys():
            # assing pubmed ids that contain unassigned to NaN value
            intact_df[self.intact_field_new_names["pubmeds"]].loc[
                intact_df[self.intact_field_new_names["pubmeds"]]
                .astype(str)
                .str.contains("unassigned", na=False)
            ] = np.nan

        # drop duplicates if same a x b pair exists multiple times
        # keep the pair with the highest score and collect pubmed ids of duplicated a x b pairs in that pair's pubmed id column
        # if a x b pair has same interaction type with b x a pair, drop b x a pair
        if "mi_score" in self.intact_field_new_names.keys():
            intact_df.sort_values(
                by=self.intact_field_new_names["mi_score"],
                ascending=False,
                inplace=True,
            )

        intact_df_unique = intact_df.dropna(
            subset=["uniprot_a", "uniprot_b"]
        ).reset_index(drop=True)

        def aggregate_pubmeds(element):
            element = "|".join([str(e) for e in set(element.dropna())])
            return np.nan if not element else element

        agg_dict = {}
        for e in self.intact_field_new_names.values():
            if e == self.intact_field_new_names["pubmeds"]:
                agg_dict[e] = aggregate_pubmeds
            else:
                agg_dict[e] = "first"

        intact_df_unique = intact_df_unique.groupby(
            ["uniprot_a", "uniprot_b"], sort=False, as_index=False
        ).aggregate(agg_dict)

        # intact_df_unique["pubmed_id"].replace("", np.nan, inplace=True) # replace empty string with NaN

        if "interaction_types" in self.intact_field_new_names.keys():
            intact_df_unique = intact_df_unique[
                ~intact_df_unique[
                    [
                        "uniprot_a",
                        "uniprot_b",
                        self.intact_field_new_names["interaction_types"],
                    ]
                ]
                .apply(frozenset, axis=1)
                .duplicated()
            ].reset_index(drop=True)
        else:
            intact_df_unique = intact_df_unique[
                ~intact_df_unique[["uniprot_a", "uniprot_b"]]
                .apply(frozenset, axis=1)
                .duplicated()
            ].reset_index(drop=True)

        t2 = time()
        logger.info(
            f"IntAct data is processed in {round((t2-t1) / 60, 2)} mins"
        )
        logger.debug(
            f"Total number of interactions for IntAct is {intact_df_unique.shape[0]}"
        )

        self.check_status_and_properties["intact"]["processed"] = True
        self.check_status_and_properties["intact"][
            "dataframe"
        ] = intact_df_unique
        self.check_status_and_properties["intact"][
            "properties_dict"
        ] = self.intact_field_new_names

    def download_biogrid_data(self) -> None:
        """
        Wrapper function to download BioGRID data using pypath; used to access
        settings.

        To do: Make arguments of biogrid.biogrid_all_interactions selectable for user.
        """

        logger.info(
            f"This is the link of BioGRID data we downloaded:{urls.urls['biogrid']['all']}. Please check if it is up to date"
        )
        logger.debug("Started downloading BioGRID data")
        t0 = time()

        # download biogrid data
        self.biogrid_ints = biogrid.biogrid_all_interactions(
            self.organism, 9999999999, False
        )

        # download these fields for mapping from gene symbol to uniprot id
        try:
            self.uniprot_to_gene = uniprot.uniprot_data("gene_names", "*", True)
            self.uniprot_to_tax = uniprot.uniprot_data("organism_id", "*", True)
            # Check if results are dicts, if not, use empty dicts
            if not isinstance(self.uniprot_to_gene, dict):
                logger.warning(f"uniprot_to_gene returned {type(self.uniprot_to_gene)} instead of dict. Using empty mapping.")
                self.uniprot_to_gene = {}
            if not isinstance(self.uniprot_to_tax, dict):
                logger.warning(f"uniprot_to_tax returned {type(self.uniprot_to_tax)} instead of dict. Using empty mapping.")
                self.uniprot_to_tax = {}
        except (ValueError, AttributeError, TypeError) as e:
            logger.warning(f"Failed to download uniprot gene/tax mapping data: {e}. Using empty mappings.")
            self.uniprot_to_gene = {}
            self.uniprot_to_tax = {}

        t1 = time()
        logger.info(
            f"BioGRID data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

        if self.test_mode:
            self.biogrid_ints = self.biogrid_ints[:100]

        self.check_status_and_properties["biogrid"]["downloaded"] = True

    @validate_call
    def biogrid_process(
        self, rename_selected_fields: dict[str, str] = None
    ) -> None:
        """
        Processor function for BioGRID data. It drops duplicate and reciprocal duplicate protein pairs and collects pubmed ids of duplicated pairs. In addition, it
        maps entries to uniprot ids using gene name and tax id information in the BioGRID data. Also, it filters protein pairs found in swissprot.

         Args:
            rename_selected_fields : List of new field names for selected fields. If not defined, default field names will be used.
        """

        if self.biogrid_fields is None:
            selected_fields = [field.value for field in BiogridEdgeField]
        else:
            selected_fields = [field.value for field in self.biogrid_fields]

        default_field_names = {
            "source": "source",
            "pmid": "pubmed_ids",
            "experimental_system": "method",
        }

        self.biogrid_field_new_names = {}

        if rename_selected_fields:
            if len(selected_fields) != len(rename_selected_fields):
                raise Exception(
                    "Length of selected_fields variable should be equal to length of rename_selected_fields variable"
                )

            for field_old_name, field_new_name in list(
                zip(selected_fields, rename_selected_fields)
            ):
                self.biogrid_field_new_names[field_old_name] = field_new_name

            self.biogrid_field_new_names["uniprot_a"] = "uniprot_a"
            self.biogrid_field_new_names["uniprot_b"] = "uniprot_b"
        else:
            for field_old_name in selected_fields:
                self.biogrid_field_new_names[field_old_name] = (
                    default_field_names[field_old_name]
                )

            self.biogrid_field_new_names["uniprot_a"] = "uniprot_a"
            self.biogrid_field_new_names["uniprot_b"] = "uniprot_b"

        logger.debug("Started processing BioGRID data")
        t1 = time()

        # create dataframe
        biogrid_df = pd.DataFrame.from_records(
            self.biogrid_ints, columns=self.biogrid_ints[0]._fields
        )

        # biogrid id (gene symbols) to uniprot id mapping
        biogrid_df["partner_a"] = biogrid_df["partner_a"].str.upper()
        biogrid_df["partner_b"] = biogrid_df["partner_b"].str.upper()

        gene_to_uniprot = collections.defaultdict(list)
        for k, v in self.uniprot_to_gene.items():
            for gene in v.split():
                gene_to_uniprot[gene.upper()].append(k)

        prot_a_uniprots = []
        for prot, tax in zip(biogrid_df["partner_a"], biogrid_df["tax_a"]):
            uniprot_id_a = (
                ";".join(
                    [
                        _id
                        for _id in gene_to_uniprot[prot]
                        if tax == self.uniprot_to_tax[_id]
                    ]
                )
                if prot in gene_to_uniprot
                else None
            )
            prot_a_uniprots.append(uniprot_id_a)

        prot_b_uniprots = []
        for prot, tax in zip(biogrid_df["partner_b"], biogrid_df["tax_b"]):
            uniprot_id_b = (
                ";".join(
                    [
                        _id
                        for _id in gene_to_uniprot[prot]
                        if tax == self.uniprot_to_tax[_id]
                    ]
                )
                if prot in gene_to_uniprot
                else None
            )
            prot_b_uniprots.append(uniprot_id_b)

        biogrid_df["uniprot_a"] = prot_a_uniprots
        biogrid_df["uniprot_b"] = prot_b_uniprots

        biogrid_df.fillna(value=np.nan, inplace=True)

        # add source database info
        biogrid_df["source"] = "BioGRID"
        # filter selected fields
        biogrid_df = biogrid_df[list(self.biogrid_field_new_names.keys())]
        # rename columns
        biogrid_df.rename(columns=self.biogrid_field_new_names, inplace=True)

        # Convert to string and drop rows with None/NaN values first
        biogrid_df = biogrid_df.dropna(subset=["uniprot_a", "uniprot_b"]).reset_index(drop=True)
        
        # Ensure columns are string type
        biogrid_df["uniprot_a"] = biogrid_df["uniprot_a"].astype(str)
        biogrid_df["uniprot_b"] = biogrid_df["uniprot_b"].astype(str)
        
        # drop rows that have semicolon (";")
        biogrid_df.drop(
            biogrid_df[
                (biogrid_df["uniprot_a"].str.contains(";", na=False))
                | (biogrid_df["uniprot_b"].str.contains(";", na=False))
            ].index,
            axis=0,
            inplace=True,
        )
        biogrid_df.reset_index(drop=True, inplace=True)

        # drop rows if uniprot_a or uniprot_b is not a swiss-prot protein
        biogrid_df = biogrid_df[
            (biogrid_df["uniprot_a"].isin(self.swissprots))
            & (biogrid_df["uniprot_b"].isin(self.swissprots))
        ]
        biogrid_df.reset_index(drop=True, inplace=True)

        # drop duplicates if same a x b pair exists multiple times
        # keep the first pair and collect pubmed ids of duplicated a x b pairs in that pair's pubmed id column
        # if a x b pair has same experimental system type with b x a pair, drop b x a pair
        biogrid_df_unique = biogrid_df.dropna(
            subset=["uniprot_a", "uniprot_b"]
        ).reset_index(drop=True)

        def aggregate_pubmeds(element):
            element = "|".join([str(e) for e in set(element.dropna())])
            return np.nan if not element else element

        agg_dict = {}
        for e in self.biogrid_field_new_names.values():
            if e == self.biogrid_field_new_names["pmid"]:
                agg_dict[e] = aggregate_pubmeds
            else:
                agg_dict[e] = "first"

        biogrid_df_unique = biogrid_df_unique.groupby(
            ["uniprot_a", "uniprot_b"], sort=False, as_index=False
        ).aggregate(agg_dict)
        # biogrid_df_unique["pubmed_id"].replace("", np.nan, inplace=True)

        if "experimental_system" in self.biogrid_field_new_names.keys():
            biogrid_df_unique = biogrid_df_unique[
                ~biogrid_df_unique[
                    [
                        "uniprot_a",
                        "uniprot_b",
                        self.biogrid_field_new_names["experimental_system"],
                    ]
                ]
                .apply(frozenset, axis=1)
                .duplicated()
            ].reset_index(drop=True)
        else:
            biogrid_df_unique = biogrid_df_unique[
                ~biogrid_df_unique[["uniprot_a", "uniprot_b"]]
                .apply(frozenset, axis=1)
                .duplicated()
            ].reset_index(drop=True)

        t2 = time()
        logger.info(
            f"BioGRID data is processed in {round((t2-t1) / 60, 2)} mins"
        )
        logger.debug(
            f"Total number of interactions for BioGRID is {biogrid_df_unique.shape[0]}"
        )

        self.check_status_and_properties["biogrid"]["processed"] = True
        self.check_status_and_properties["biogrid"][
            "dataframe"
        ] = biogrid_df_unique
        self.check_status_and_properties["biogrid"][
            "properties_dict"
        ] = self.biogrid_field_new_names

    def download_string_data(self) -> None:
        """
        Wrapper function to download STRING data using pypath; used to access
        settings.

        To do: Make arguments of string.string_links_interactions selectable for user.
        """

        t0 = time()

        if self.organism is None:
            string_species = string.string_species()
            self.tax_ids = list(string_species.keys())
        else:
            self.tax_ids = [self.organism]

        # map string ids to swissprot ids
        try:
            uniprot_to_string = uniprot.uniprot_data("xref_string", "*", True)
            # Check if result is a dict, if not (e.g., empty list), use empty dict
            if not isinstance(uniprot_to_string, dict):
                logger.warning(f"uniprot_to_string returned {type(uniprot_to_string)} instead of dict. Using empty mapping.")
                uniprot_to_string = {}
        except (ValueError, AttributeError, TypeError) as e:
            logger.warning(f"Failed to download uniprot-string mapping: {e}. Using empty mapping.")
            uniprot_to_string = {}

        self.string_to_uniprot = collections.defaultdict(list)
        for k, v in uniprot_to_string.items():
            if v:  # Check if v is not None or empty
                for string_id in list(filter(None, str(v).split(";"))):
                    if "." in string_id:
                        self.string_to_uniprot[string_id.split(".")[1]].append(k)

        self.string_ints = []

        logger.debug("Started downloading STRING data")
        logger.info(
            f"This is the link of STRING data we downloaded:{urls.urls['string']['links']}. Please check if it is up to date"
        )

        # this tax id give an error
        tax_ids_to_be_skipped = [
            "4565",
            "8032",
        ]

        # Filter out problematic tax IDs
        valid_tax_ids = [tax for tax in self.tax_ids if tax not in tax_ids_to_be_skipped]
        
        # In test mode, only download a subset
        if self.test_mode:
            valid_tax_ids = valid_tax_ids[:10]  # Only process first 10 species in test mode
            logger.info(f"Test mode: limiting to first 10 species out of {len(self.tax_ids)} total")

        # Parallel download function
        from concurrent.futures import ThreadPoolExecutor, as_completed
        from threading import Lock
        
        results_lock = Lock()
        self.string_ints = []
        
        def download_single_species(tax):
            """Download STRING data for a single species"""
            try:
                organism_string_ints = [
                    i
                    for i in string.string_links_interactions(
                        ncbi_tax_id=int(tax),
                        score_threshold="high_confidence",
                    )
                    if i.protein_a in self.string_to_uniprot
                    and i.protein_b in self.string_to_uniprot
                ]
                
                return (tax, organism_string_ints)
            except Exception as e:
                logger.warning(f"Failed to download STRING data for tax_id {tax}: {e}")
                return (tax, [])
        
        # Use parallel downloads with up to 8 workers
        max_workers = min(8, len(valid_tax_ids))
        logger.info(f"Downloading STRING data for {len(valid_tax_ids)} species using {max_workers} parallel workers...")
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all download tasks
            future_to_tax = {executor.submit(download_single_species, tax): tax for tax in valid_tax_ids}
            
            # Process completed downloads with progress bar
            for future in tqdm(as_completed(future_to_tax), total=len(valid_tax_ids), desc="Downloading STRING"):
                tax, organism_string_ints = future.result()
                
                if organism_string_ints:
                    with results_lock:
                        self.string_ints.extend(organism_string_ints)
                    
                    logger.debug(
                        f"Downloaded STRING data for tax_id {tax}: {len(organism_string_ints)} interactions (total: {len(self.string_ints)})"
                    )

        t1 = time()
        logger.info(
            f"STRING data is downloaded in {round((t1-t0) / 60, 2)} mins"
        )

        if self.test_mode:
            self.string_ints = self.string_ints[:100]

        self.check_status_and_properties["string"]["downloaded"] = True

    @validate_call
    def string_process(
        self, rename_selected_fields: dict[str, str] = None
    ) -> None:
        """
        Processor function for STRING data. It drops duplicate and reciprocal duplicate protein pairs. In addition, it maps entries to uniprot ids
        using crossreferences to STRING in the Uniprot data. Also, it filters protein pairs found in swissprot.

         Args:
            rename_selected_fields : List of new field names for selected fields. If not defined, default field names will be used.
        """

        if self.string_fields is None:
            selected_fields = [field.value for field in StringEdgeField]
        else:
            selected_fields = [field.value for field in self.string_fields]

        default_field_names = {
            "source": "source",
            "combined_score": "string_combined_score",
            "physical_combined_score": "string_physical_combined_score",
        }

        self.string_field_new_names = {}

        if rename_selected_fields:
            if len(selected_fields) != len(rename_selected_fields):
                raise Exception(
                    "Length of selected_fields variable should be equal to length of rename_selected_fields variable"
                )

            for field_old_name, field_new_name in list(
                zip(selected_fields, rename_selected_fields)
            ):
                self.string_field_new_names[field_old_name] = field_new_name

            self.string_field_new_names["uniprot_a"] = "uniprot_a"
            self.string_field_new_names["uniprot_b"] = "uniprot_b"
        else:
            for field_old_name in selected_fields:
                self.string_field_new_names[field_old_name] = (
                    default_field_names[field_old_name]
                )

            self.string_field_new_names["uniprot_a"] = "uniprot_a"
            self.string_field_new_names["uniprot_b"] = "uniprot_b"

        logger.debug("Started processing STRING data")
        t1 = time()

        # Check if string_ints is empty
        if not self.string_ints:
            logger.warning("No STRING interactions found. This might be due to failed uniprot-string mapping.")
            logger.warning("Creating empty STRING dataframe.")
            
            # Create empty dataframe with expected structure
            string_df_unique = pd.DataFrame(columns=list(self.string_field_new_names.values()))
            
            t2 = time()
            logger.info(f"STRING data is processed in {round((t2-t1) / 60, 2)} mins")
            logger.debug(f"Total number of interactions for STRING is {string_df_unique.shape[0]}")
            
            self.check_status_and_properties["string"]["processed"] = True
            self.check_status_and_properties["string"]["dataframe"] = string_df_unique
            self.check_status_and_properties["string"]["properties_dict"] = self.string_field_new_names
            return

        # create dataframe
        string_df = pd.DataFrame.from_records(
            self.string_ints, columns=self.string_ints[0]._fields
        )

        prot_a_uniprots = []
        for protein in string_df["protein_a"]:
            id_a = ";".join(self.string_to_uniprot[protein])
            # if protein in self.string_to_uniprot else None)
            # now that we filtered interactions in line 307, we should not get KeyError here
            prot_a_uniprots.append(id_a)

        prot_b_uniprots = []
        for protein in string_df["protein_b"]:
            id_b = ";".join(self.string_to_uniprot[protein])
            # if protein in self.string_to_uniprot else None)
            # now that we filtered interactions in line 307, we should not get KeyError here
            prot_b_uniprots.append(id_b)

        string_df["uniprot_a"] = prot_a_uniprots
        string_df["uniprot_b"] = prot_b_uniprots

        string_df.fillna(value=np.nan, inplace=True)

        # add source database info
        string_df["source"] = "STRING"
        # filter selected fields
        string_df = string_df[list(self.string_field_new_names.keys())]
        # rename columns
        string_df.rename(columns=self.string_field_new_names, inplace=True)

        # filter with swissprot ids
        # we already filtered interactions in line 307, we can remove this part or keep it for a double check
        # string_df = string_df[(string_df["uniprot_a"].isin(self.swissprots)) & (string_df["uniprot_b"].isin(self.swissprots))]
        # string_df.reset_index(drop=True, inplace=True)

        # drop duplicates if same a x b pair exists in b x a format
        # keep the one with the highest combined score
        if "combined_score" in self.string_field_new_names.keys():
            string_df.sort_values(
                by=self.string_field_new_names["combined_score"],
                ascending=False,
                inplace=True,
            )
            string_df_unique = (
                string_df.dropna(subset=["uniprot_a", "uniprot_b"])
                .drop_duplicates(
                    subset=["uniprot_a", "uniprot_b"], keep="first"
                )
                .reset_index(drop=True)
            )
            string_df_unique = string_df_unique[
                ~string_df_unique[
                    [
                        "uniprot_a",
                        "uniprot_b",
                        self.string_field_new_names["combined_score"],
                    ]
                ]
                .apply(frozenset, axis=1)
                .duplicated()
            ].reset_index(drop=True)
        else:
            string_df_unique = string_df_unique[
                ~string_df_unique[["uniprot_a", "uniprot_b"]]
                .apply(frozenset, axis=1)
                .duplicated()
            ].reset_index(drop=True)

        t2 = time()
        logger.info(
            f"STRING data is processed in {round((t2-t1) / 60, 2)} mins"
        )
        logger.debug(
            f"Total number of interactions for STRING is {string_df_unique.shape[0]}"
        )

        self.check_status_and_properties["string"]["processed"] = True
        self.check_status_and_properties["string"][
            "dataframe"
        ] = string_df_unique
        self.check_status_and_properties["string"][
            "properties_dict"
        ] = self.string_field_new_names

    def merge_all(self) -> pd.DataFrame:
        """
        Merge function for all 3 databases. Merge dataframes according to uniprot_a and uniprot_b (i.e., protein pairs) columns.

        """

        t1 = time()
        logger.debug(
            "started merging interactions from all 3 databases (IntAct, BioGRID, STRING)"
        )

        def merge_pubmed_ids(elem):
            """
            Merges pubmed id columns
            """
            if len(elem.dropna().tolist()) > 0:
                new_list = []
                for e in elem.dropna().tolist():
                    if "|" in e:
                        new_list.extend(e.split("|"))
                    else:
                        new_list.append(e)

                return "|".join(list(set(new_list)))
            else:
                return np.nan

        # during the merging, it changes datatypes of some columns from int to float. So it needs to be reverted
        def float_to_int(element):
            """
            Forces to change data type from float to int in dataframe
            """
            if "." in str(element):
                dot_index = str(element).index(".")
                element = str(element)[:dot_index]
                return element
            else:
                return element

        # check which databases will be merged
        dbs_will_be_merged = []
        for db in self.check_status_and_properties.keys():
            if (
                self.check_status_and_properties[db]["downloaded"]
                and self.check_status_and_properties[db]["processed"]
            ):
                dbs_will_be_merged.append(db)

        seen_dbs = set()
        for db in dbs_will_be_merged:

            if db in seen_dbs:
                continue

            if dbs_will_be_merged.index(db) == 0:
                if db == "intact" and dbs_will_be_merged[1] == "biogrid":
                    seen_dbs.add(db)
                    seen_dbs.add(dbs_will_be_merged[1])

                    df1 = self.check_status_and_properties[db]["dataframe"]
                    df2 = self.check_status_and_properties[
                        dbs_will_be_merged[1]
                    ]["dataframe"]

                    merged_df = pd.merge(
                        df1, df2, on=["uniprot_a", "uniprot_b"], how="outer"
                    )

                    # if source column exists in both intact and biogrid merge them
                    if self.intact_field_new_names.get(
                        "source", None
                    ) and self.biogrid_field_new_names.get("source", None):
                        # if they have the same name
                        if (
                            self.intact_field_new_names["source"]
                            == self.biogrid_field_new_names["source"]
                        ):
                            merged_df["source"] = merged_df[
                                [
                                    self.intact_field_new_names["source"]
                                    + "_x",
                                    self.biogrid_field_new_names["source"]
                                    + "_y",
                                ]
                            ].apply(lambda x: "|".join(x.dropna()), axis=1)

                            merged_df.drop(
                                columns=[
                                    self.intact_field_new_names["source"]
                                    + "_x",
                                    self.biogrid_field_new_names["source"]
                                    + "_y",
                                ],
                                inplace=True,
                            )

                        # if they dont have the same name
                        else:
                            merged_df["source"] = merged_df[
                                [
                                    self.intact_field_new_names["source"],
                                    self.biogrid_field_new_names["source"],
                                ]
                            ].apply(lambda x: "|".join(x.dropna()), axis=1)

                            merged_df.drop(
                                columns=[
                                    self.intact_field_new_names["source"],
                                    self.biogrid_field_new_names["source"],
                                ],
                                inplace=True,
                            )

                    # if pubmeds and pmid column exist in intact and biogrid merge them
                    if self.intact_field_new_names.get(
                        "pubmeds", None
                    ) and self.biogrid_field_new_names.get("pmid", None):
                        # if they have the same name
                        if (
                            self.intact_field_new_names["pubmeds"]
                            == self.biogrid_field_new_names["pmid"]
                        ):
                            merged_df["pubmed_ids"] = merged_df[
                                [
                                    self.intact_field_new_names["pubmeds"]
                                    + "_x",
                                    self.biogrid_field_new_names["pmid"] + "_y",
                                ]
                            ].apply(merge_pubmed_ids, axis=1)

                            merged_df.drop(
                                columns=[
                                    self.intact_field_new_names["pubmeds"]
                                    + "_x",
                                    self.biogrid_field_new_names["pmid"] + "_y",
                                ],
                                inplace=True,
                            )

                        # if they dont have the same name
                        else:
                            merged_df["pubmed_ids"] = merged_df[
                                [
                                    self.intact_field_new_names["pubmeds"],
                                    self.biogrid_field_new_names["pmid"],
                                ]
                            ].apply(merge_pubmed_ids, axis=1)

                            merged_df.drop(
                                columns=[
                                    self.intact_field_new_names["pubmeds"],
                                    self.biogrid_field_new_names["pmid"],
                                ],
                                inplace=True,
                            )

                    # if methods and experimental_system column exist in intact and biogrid merge them
                    if self.intact_field_new_names.get(
                        "methods", None
                    ) and self.biogrid_field_new_names.get(
                        "experimental_system", None
                    ):
                        # if they have the same name
                        if (
                            self.intact_field_new_names["methods"]
                            == self.biogrid_field_new_names[
                                "experimental_system"
                            ]
                        ):
                            merged_df["method"] = merged_df[
                                [
                                    self.intact_field_new_names["methods"]
                                    + "_x",
                                    self.biogrid_field_new_names[
                                        "experimental_system"
                                    ]
                                    + "_y",
                                ]
                            ].apply(
                                lambda x: (
                                    x.dropna().tolist()[0]
                                    if len(x.dropna().tolist()) > 0
                                    else np.nan
                                ),
                                axis=1,
                            )

                            merged_df.drop(
                                columns=[
                                    self.intact_field_new_names["methods"]
                                    + "_x",
                                    self.biogrid_field_new_names[
                                        "experimental_system"
                                    ]
                                    + "_y",
                                ],
                                inplace=True,
                            )
                        # if they dont have the same name
                        else:
                            merged_df["method"] = merged_df[
                                [
                                    self.intact_field_new_names["methods"],
                                    self.biogrid_field_new_names[
                                        "experimental_system"
                                    ],
                                ]
                            ].apply(
                                lambda x: (
                                    x.dropna().tolist()[0]
                                    if len(x.dropna().tolist()) > 0
                                    else np.nan
                                ),
                                axis=1,
                            )

                            merged_df.drop(
                                columns=[
                                    self.intact_field_new_names["methods"],
                                    self.biogrid_field_new_names[
                                        "experimental_system"
                                    ],
                                ],
                                inplace=True,
                            )

                else:
                    seen_dbs.add(db)
                    seen_dbs.add(dbs_will_be_merged[1])

                    df1 = self.check_status_and_properties[db]["dataframe"]
                    df2 = self.check_status_and_properties[
                        dbs_will_be_merged[1]
                    ]["dataframe"]

                    merged_df = pd.merge(
                        df1, df2, on=["uniprot_a", "uniprot_b"], how="outer"
                    )

                    # if source column exists in both intact and string merge them
                    if self.check_status_and_properties[db][
                        "properties_dict"
                    ].get("source", None) and self.check_status_and_properties[
                        dbs_will_be_merged[1]
                    ][
                        "properties_dict"
                    ].get(
                        "source", None
                    ):
                        # if they have the same name
                        if (
                            self.check_status_and_properties[db][
                                "properties_dict"
                            ]["source"]
                            == self.check_status_and_properties[
                                dbs_will_be_merged[1]
                            ]["properties_dict"]["source"]
                        ):
                            merged_df["source"] = merged_df[
                                [
                                    self.check_status_and_properties[db][
                                        "properties_dict"
                                    ]["source"]
                                    + "_x",
                                    self.check_status_and_properties[
                                        dbs_will_be_merged[1]
                                    ]["properties_dict"]["source"]
                                    + "_y",
                                ]
                            ].apply(lambda x: "|".join(x.dropna()), axis=1)

                            merged_df.drop(
                                columns=[
                                    self.check_status_and_properties[db][
                                        "properties_dict"
                                    ]["source"]
                                    + "_x",
                                    self.check_status_and_properties[
                                        dbs_will_be_merged[1]
                                    ]["properties_dict"]["source"]
                                    + "_y",
                                ],
                                inplace=True,
                            )

                        # if they dont have the same name
                        else:
                            merged_df["source"] = merged_df[
                                [
                                    self.check_status_and_properties[db][
                                        "properties_dict"
                                    ]["source"],
                                    self.check_status_and_properties[
                                        dbs_will_be_merged[1]
                                    ]["properties_dict"]["source"],
                                ]
                            ].apply(lambda x: "|".join(x.dropna()), axis=1)

                            merged_df.drop(
                                columns=[
                                    self.check_status_and_properties[db][
                                        "properties_dict"
                                    ]["source"],
                                    self.check_status_and_properties[
                                        dbs_will_be_merged[1]
                                    ]["properties_dict"]["source"],
                                ],
                                inplace=True,
                            )

                    # if combined_score field exists in dataframe force its data data type become int
                    if self.check_status_and_properties[dbs_will_be_merged[1]][
                        "properties_dict"
                    ].get("combined_score", None):
                        merged_df[
                            self.string_field_new_names["combined_score"]
                        ] = merged_df[
                            self.string_field_new_names["combined_score"]
                        ].astype(
                            str, errors="ignore"
                        )
                        merged_df[
                            self.string_field_new_names["combined_score"]
                        ] = merged_df[
                            self.string_field_new_names["combined_score"]
                        ].apply(
                            float_to_int
                        )

                    # if physical_combined_score field exists in dataframe force its data data type become int
                    if self.check_status_and_properties[dbs_will_be_merged[1]][
                        "properties_dict"
                    ].get("physical_combined_score", None):
                        merged_df[
                            self.string_field_new_names[
                                "physical_combined_score"
                            ]
                        ] = merged_df[
                            self.string_field_new_names[
                                "physical_combined_score"
                            ]
                        ].astype(
                            str, errors="ignore"
                        )
                        merged_df[
                            self.string_field_new_names[
                                "physical_combined_score"
                            ]
                        ] = merged_df[
                            self.string_field_new_names[
                                "physical_combined_score"
                            ]
                        ].apply(
                            float_to_int
                        )

            else:
                seen_dbs.add(db)

                df2 = self.check_status_and_properties[db]["dataframe"]
                source_flag = "source" in list(merged_df.columns)

                merged_df = pd.merge(
                    merged_df, df2, on=["uniprot_a", "uniprot_b"], how="outer"
                )

                # if source column exists in both merged_df and string merge them
                if source_flag and self.string_field_new_names.get(
                    "source", None
                ):

                    # if they have the same name
                    if "source" == self.string_field_new_names["source"]:
                        merged_df["source"] = merged_df[
                            [
                                "source_x",
                                self.string_field_new_names["source"] + "_y",
                            ]
                        ].apply(lambda x: "|".join(x.dropna()), axis=1)

                        merged_df.drop(
                            columns=[
                                "source_x",
                                self.string_field_new_names["source"] + "_y",
                            ],
                            inplace=True,
                        )

                    # if they dont have the same name
                    else:
                        merged_df["source"] = merged_df[
                            ["source", self.string_field_new_names["source"]]
                        ].apply(lambda x: "|".join(x.dropna()), axis=1)

                        merged_df.drop(
                            columns=[
                                "source",
                                self.string_field_new_names["source"],
                            ],
                            inplace=True,
                        )

                # if combined_score field exists in dataframe force its data data type become int
                if self.string_field_new_names.get("combined_score", None):
                    merged_df[self.string_field_new_names["combined_score"]] = (
                        merged_df[
                            self.string_field_new_names["combined_score"]
                        ].astype(str, errors="ignore")
                    )
                    merged_df[self.string_field_new_names["combined_score"]] = (
                        merged_df[
                            self.string_field_new_names["combined_score"]
                        ].apply(float_to_int)
                    )

                # if physical_combined_score field exists in dataframe force its data data type become int
                if self.string_field_new_names.get(
                    "physical_combined_score", None
                ):
                    merged_df[
                        self.string_field_new_names["physical_combined_score"]
                    ] = merged_df[
                        self.string_field_new_names["physical_combined_score"]
                    ].astype(
                        str, errors="ignore"
                    )
                    merged_df[
                        self.string_field_new_names["physical_combined_score"]
                    ] = merged_df[
                        self.string_field_new_names["physical_combined_score"]
                    ].apply(
                        float_to_int
                    )

        logger.debug("Merged all interactions")
        t2 = time()
        logger.info(
            f"All data is merged and processed in {round((t2-t1) / 60, 2)} mins"
        )
        logger.debug(
            f"Total number of interactions for PPI data is {merged_df.shape[0]}"
        )

        if self.export_csv:
            if self.output_dir:
                full_path = os.path.join(self.output_dir, "PPI.csv")
            else:
                full_path = os.path.join(os.getcwd(), "PPI.csv")

            merged_df.to_csv(full_path, index=False)
            logger.info(f"PPI data is written: {full_path}")

        return merged_df

    @validate_call
    def add_prefix_to_id(
        self, prefix: str = "uniprot", identifier: str = None, sep: str = ":"
    ) -> str:
        """
        Adds prefix to uniprot id
        """
        if self.add_prefix and identifier:
            return normalize_curie(prefix + sep + identifier)

        return identifier

    @validate_call
    def get_ppi_edges(
        self, label: str = "Protein_interacts_with_protein"
    ) -> list[tuple]:
        """
        Get PPI edges from merged data
        Args:
            label: label of protein-protein interaction edges
        """
        merged_df = self.merge_all()

        # create edge list
        edge_list = []
        for _, row in tqdm(merged_df.iterrows()):
            _dict = row.to_dict()

            _source = self.add_prefix_to_id(identifier=str(row["uniprot_a"]))
            _target = self.add_prefix_to_id(identifier=str(row["uniprot_b"]))

            del _dict["uniprot_a"], _dict["uniprot_b"]

            _props = {}
            for k, v in _dict.items():
                if str(v) != "nan":
                    if isinstance(v, str) and "|" in v:
                        _props[str(k).replace(" ", "_").lower()] = v.replace(
                            "'", "^"
                        ).split("|")
                    else:
                        _props[str(k).replace(" ", "_").lower()] = str(
                            v
                        ).replace("'", "^")

            edge_list.append((None, _source, _target, label, _props))

        return edge_list
