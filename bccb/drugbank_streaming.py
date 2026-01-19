"""
Memory-efficient DrugBank XML parser using iterative parsing.
Replaces pypath's DrugbankFull which loads entire 1.5GB XML into memory.
"""
from __future__ import annotations

import os
import collections
from zipfile import ZipFile
from lxml import etree
from typing import Optional
from pypath.share import curl
from pypath.resources import urls

from biocypher._logger import logger


def _drugbank_download(
    url: str,
    user: str,
    passwd: str,
    credentials_fname: Optional[str] = None,
):
    """Download DrugBank file with authentication."""
    from pypath.inputs.drugbank import _drugbank_download as pypath_download
    return pypath_download(
        url=url,
        user=user,
        passwd=passwd,
        credentials_fname=credentials_fname,
    )


class DrugbankStreaming:
    """
    Memory-efficient DrugBank parser using iterparse.
    Processes drugs one at a time instead of loading entire XML.
    """
    
    NS = {'db': 'http://www.drugbank.ca'}
    
    def __init__(
        self,
        user: Optional[str] = None,
        passwd: Optional[str] = None,
        credentials_fname: Optional[str] = None,
    ):
        self.user = user
        self.passwd = passwd
        self.credentials_fname = credentials_fname
        self._xml_path = None
        self._ensure_xml_extracted()
    
    def _ensure_xml_extracted(self):
        """Download and extract XML if needed."""
        # First check if XML already exists in cache
        import pypath.share.settings as settings_module
        cache_dir = settings_module.get('cachedir')
        xml_path = os.path.join(cache_dir, 'full database.xml')
        
        if os.path.exists(xml_path):
            self._xml_path = xml_path
            logger.debug(f"DrugBank XML already exists: {self._xml_path}")
            return
        
        # Otherwise download
        result = _drugbank_download(
            url=urls.urls['drugbank']['full_database'],
            user=self.user,
            passwd=self.passwd,
            credentials_fname=self.credentials_fname,
        )
        
        if result is None or result.fileobj is None:
            # Check if file exists anyway (from previous download)
            if os.path.exists(xml_path):
                self._xml_path = xml_path
                logger.debug(f"DrugBank XML found in cache: {self._xml_path}")
                return
            raise RuntimeError("Failed to download DrugBank data")
        
        path = result.fileobj.name
        xml_path = os.path.join(os.path.dirname(path), 'full database.xml')
        
        if not os.path.exists(xml_path):
            with ZipFile(path, 'r') as zip_ref:
                zip_ref.extractall(os.path.dirname(path))
        
        self._xml_path = xml_path
        logger.debug(f"DrugBank XML path: {self._xml_path}")
    
    def iter_drugs(self):
        """
        Iterate over drugs in the XML file without loading entire file.
        Yields drug elements one at a time.
        """
        context = etree.iterparse(
            self._xml_path,
            events=('end',),
            tag='{http://www.drugbank.ca}drug',
        )
        
        for event, elem in context:
            # Only process top-level drug elements
            if elem.getparent() is not None and elem.getparent().tag == '{http://www.drugbank.ca}drugbank':
                yield elem
                # Clear element to free memory
                elem.clear()
                # Also clear parent's references to processed children
                while elem.getprevious() is not None:
                    del elem.getparent()[0]
    
    def drugbank_drugs_full(
        self,
        fields: list[str] | None = None,
    ) -> list[tuple]:
        """
        Returns list of namedtuples with drug information.
        Memory-efficient: processes one drug at a time.
        """
        basic_fields = [
            'drugbank_id', 'type', 'name', 'description', 'cas_number', 'unii',
            'average_mass', 'monoisotopic_mass', 'state', 'groups',
            'general_references', 'atc_codes',
        ]
        
        if fields is None:
            fields = basic_fields
        
        # Ensure drugbank_id is always included
        if 'drugbank_id' not in fields:
            fields = ['drugbank_id'] + list(fields)
        
        DrugbankDrug = collections.namedtuple('DrugbankDrug', fields)
        result = []
        
        for drug in self.iter_drugs():
            drug_data = {}
            
            for field in fields:
                if field == 'drugbank_id':
                    id_elem = drug.find('db:drugbank-id[@primary="true"]', self.NS)
                    drug_data[field] = id_elem.text if id_elem is not None else None
                elif field == 'type':
                    drug_data[field] = drug.get('type')
                elif field == 'groups':
                    groups = drug.findall('db:groups/db:group', self.NS)
                    drug_data[field] = [g.text for g in groups]
                elif field == 'general_references':
                    refs = drug.findall('db:general-references/db:articles/db:article/db:pubmed-id', self.NS)
                    drug_data[field] = [r.text for r in refs if r.text]
                elif field == 'atc_codes':
                    codes = drug.findall('db:atc-codes/db:atc-code', self.NS)
                    drug_data[field] = [c.get('code') for c in codes if c.get('code')]
                else:
                    elem = drug.find(f'db:{field.replace("_", "-")}', self.NS)
                    drug_data[field] = elem.text if elem is not None else None
            
            result.append(DrugbankDrug(**drug_data))
        
        return result
    
    def drugbank_external_ids_full(self) -> dict:
        """
        Returns dict mapping drugbank_id to external identifiers.
        """
        result = {}
        
        for drug in self.iter_drugs():
            drug_id_elem = drug.find('db:drugbank-id[@primary="true"]', self.NS)
            if drug_id_elem is None:
                continue
            drug_id = drug_id_elem.text
            
            external_ids = {}
            for ext_id in drug.findall('db:external-identifiers/db:external-identifier', self.NS):
                resource = ext_id.find('db:resource', self.NS)
                identifier = ext_id.find('db:identifier', self.NS)
                if resource is not None and identifier is not None:
                    external_ids[resource.text] = identifier.text
            
            result[drug_id] = external_ids
        
        return result
    
    def drugbank_properties_full(self) -> dict:
        """
        Returns dict mapping drugbank_id to calculated properties (SMILES, InChI, etc.)
        """
        result = {}
        
        for drug in self.iter_drugs():
            drug_id_elem = drug.find('db:drugbank-id[@primary="true"]', self.NS)
            if drug_id_elem is None:
                continue
            drug_id = drug_id_elem.text
            
            properties = {}
            for prop in drug.findall('db:calculated-properties/db:property', self.NS):
                kind = prop.find('db:kind', self.NS)
                value = prop.find('db:value', self.NS)
                if kind is not None and value is not None:
                    properties[kind.text] = value.text
            
            result[drug_id] = properties
        
        return result
    
    def drugbank_targets_full(self, fields: list = None) -> list[tuple]:
        """
        Returns list of drug-target interactions.
        The fields parameter is accepted for compatibility but not used.
        """
        field_names = ['drugbank_id', 'target_id', 'target_name', 'target_uniprot', 
                  'actions', 'known_action', 'organism', 'references', 'polypeptide']
        DrugbankTarget = collections.namedtuple('DrugbankTarget', field_names)
        result = []
        
        for drug in self.iter_drugs():
            drug_id_elem = drug.find('db:drugbank-id[@primary="true"]', self.NS)
            if drug_id_elem is None:
                continue
            drug_id = drug_id_elem.text
            
            for target in drug.findall('db:targets/db:target', self.NS):
                target_id_elem = target.find('db:id', self.NS)
                target_name_elem = target.find('db:name', self.NS)
                organism_elem = target.find('db:organism', self.NS)
                known_action_elem = target.find('db:known-action', self.NS)
                
                # Get UniProt IDs and polypeptide info
                uniprot_ids = []
                polypeptides = []
                for polypep in target.findall('db:polypeptide', self.NS):
                    uniprot_id = polypep.get('id')
                    if uniprot_id:
                        uniprot_ids.append(uniprot_id)
                        polypeptides.append(polypep)
                
                # Get actions
                actions = [a.text for a in target.findall('db:actions/db:action', self.NS)]
                
                # Get references (pubmed IDs)
                references = []
                for ref in target.findall('db:references/db:articles/db:article', self.NS):
                    pubmed = ref.find('db:pubmed-id', self.NS)
                    if pubmed is not None and pubmed.text:
                        references.append(pubmed.text)
                
                for idx, uniprot_id in enumerate(uniprot_ids):
                    result.append(DrugbankTarget(
                        drugbank_id=drug_id,
                        target_id=target_id_elem.text if target_id_elem is not None else None,
                        target_name=target_name_elem.text if target_name_elem is not None else None,
                        target_uniprot=uniprot_id,
                        actions=actions,
                        known_action=known_action_elem.text if known_action_elem is not None else None,
                        organism=organism_elem.text if organism_elem is not None else None,
                        references=references,
                        polypeptide=polypeptides[idx] if idx < len(polypeptides) else None,
                    ))
        
        return result
