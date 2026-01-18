"""
PyPath compatibility layer for handling API changes between versions.

This module provides wrapper functions that work with both old and new pypath APIs.
"""
from collections import namedtuple
import logging

logger = logging.getLogger(__name__)

# Try to import pypath modules
try:
    from pypath.inputs import chembl
except ImportError:
    chembl = None

try:
    from pypath.inputs import unichem as unichem_module
except ImportError:
    unichem_module = None


# UniChem ID type mapping (name -> numeric ID)
UNICHEM_ID_TYPES = {
    'chembl': 1,
    'drugbank': 2,
    'pdb': 3,
    'gtopdb': 4,
    'pubchem_dotf': 5,
    'kegg': 6,
    'chebi': 7,
    'nih_ncc': 8,
    'zinc': 9,
    'emolecules': 10,
    'atlas': 12,
    'fda_srs': 14,
    'surechembl': 15,
    'pharmgkb': 17,
    'hmdb': 18,
    'selleck': 20,
    'pubchem_tp': 21,
    'pubchem': 22,
    'mcule': 23,
    'nmrshiftdb': 24,
    'lincs': 25,
    'actor': 26,
    'recon': 27,
    'molport': 28,
    'nikkaji': 29,
    'bindingdb': 31,
    'epa_comptox': 32,
    'lipidmaps': 33,
    'drugcentral': 34,
    'carotenoiddb': 35,
    'metabolights': 36,
    'brenda': 37,
    'rhea': 38,
    'chemicalbook': 39,
    'swisslipids': 41,
    'dailymed': 45,
    'clinicaltrials': 46,
    'rxnorm': 47,
    'medchemexpress': 48,
    'probesdrugs': 49,
    'ccdc': 50,
}


def unichem_mapping(id_type_a, id_type_b):
    """
    Get UniChem ID mapping between two databases.
    
    Accepts both string names (old API) and numeric IDs (new API).
    """
    if unichem_module is None:
        raise ImportError("pypath.inputs.unichem not available")
    
    # Convert string names to numeric IDs if needed
    if isinstance(id_type_a, str):
        id_type_a = UNICHEM_ID_TYPES.get(id_type_a.lower())
        if id_type_a is None:
            raise ValueError(f"Unknown UniChem ID type: {id_type_a}")
    
    if isinstance(id_type_b, str):
        id_type_b = UNICHEM_ID_TYPES.get(id_type_b.lower())
        if id_type_b is None:
            raise ValueError(f"Unknown UniChem ID type: {id_type_b}")
    
    return unichem_module.unichem_mapping(id_type_a, id_type_b)


# Define namedtuples that match the OLD API structure for backward compatibility
# These will be populated from new API data

OldChemblActivity = namedtuple('ChemblActivity', [
    'chembl',           # molecule_chembl_id
    'target_chembl',    # target_chembl_id
    'assay_chembl',     # assay_id (was assay_chembl_id)
    'document',         # document_chembl_id
    'standard_type',    # standard_type
    'standard_relation', # standard_relation
    'standard_value',   # standard_value
    'standard_units',   # standard_units
    'pchembl',          # pchembl_value
    'activity_id',      # activity_id
])

OldChemblTarget = namedtuple('ChemblTarget', [
    'target_chembl_id',  # chembl_id
    'target_type',       # target_type
    'preferred_name',    # preferred_name
    'accession',         # from components
    'organism',          # organism
    'tax_id',            # ncbi_taxa_id
])

OldChemblAssay = namedtuple('ChemblAssay', [
    'assay_chembl_id',   # assay_chembl_id
    'assay_type',        # assay_type
    'confidence_score',  # confidence_score
    'target_chembl_id',  # target_chembl_id
])

OldChemblMechanism = namedtuple('ChemblMechanism', [
    'chembl',            # molecule_chembl_id
    'action_type',       # action_type
    'direct_interaction', # direct_interaction
    'mechanism_of_action', # NOT in new API - use None
    'disease_efficacy',  # NOT in new API - use None
    'target_chembl_id',  # target_chembl_id
])


def _has_new_api():
    """Check if we have the new pypath API (functions without 'chembl_' prefix)."""
    if chembl is None:
        return False
    return hasattr(chembl, 'molecule') and not hasattr(chembl, 'chembl_molecules')


def chembl_molecules(max_pages=None):
    """
    Get ChEMBL molecules.
    
    Returns a generator/list of molecule records.
    """
    if chembl is None:
        raise ImportError("pypath.inputs.chembl not available")
    
    if _has_new_api():
        # New API: chembl.molecule()
        logger.debug("Using new pypath API: chembl.molecule()")
        return chembl.molecule(max_pages=max_pages)
    else:
        # Old API: chembl.chembl_molecules()
        logger.debug("Using old pypath API: chembl.chembl_molecules()")
        return chembl.chembl_molecules()


def chembl_activities(standard_relation=None, max_pages=None):
    """
    Get ChEMBL activities.
    
    Returns a generator/list of activity records adapted to old API format.
    """
    if chembl is None:
        raise ImportError("pypath.inputs.chembl not available")
    
    if _has_new_api():
        # New API: chembl.activity()
        logger.debug("Using new pypath API: chembl.activity()")
        try:
            for act in chembl.activity(max_pages=max_pages):
                # Filter by standard_relation if specified
                if standard_relation and act.standard_relation != standard_relation:
                    continue
                
                # Convert to old format
                yield OldChemblActivity(
                    chembl=act.molecule_chembl_id,
                    target_chembl=act.target_chembl_id,
                    assay_chembl=act.assay_id,
                    document=act.document_chembl_id,
                    standard_type=act.standard_type,
                    standard_relation=act.standard_relation,
                    standard_value=act.standard_value,
                    standard_units=act.standard_units,
                    pchembl=act.pchembl_value,
                    activity_id=act.activity_id,
                )
        except (AttributeError, TypeError) as e:
            logger.warning(f"ChemBL activity download failed: {e}. Returning empty results.")
            return
    else:
        # Old API: chembl.chembl_activities()
        logger.debug("Using old pypath API: chembl.chembl_activities()")
        if standard_relation:
            return chembl.chembl_activities(standard_relation=standard_relation)
        else:
            return chembl.chembl_activities()


def chembl_targets(max_pages=None):
    """
    Get ChEMBL targets.
    
    Returns a generator/list of target records adapted to old API format.
    """
    if chembl is None:
        raise ImportError("pypath.inputs.chembl not available")
    
    if _has_new_api():
        # New API: chembl.target()
        logger.debug("Using new pypath API: chembl.target()")
        for target in chembl.target(max_pages=max_pages):
            # Extract UniProt accession from components
            accession = None
            if target.components:
                for comp in target.components:
                    # New API uses ChemblComponent namedtuple with uniprot_accession
                    if hasattr(comp, 'uniprot_accession') and comp.uniprot_accession:
                        accession = comp.uniprot_accession
                        break
                    elif isinstance(comp, dict):
                        accession = comp.get('uniprot_accession') or comp.get('accession')
                        if accession:
                            break
            
            yield OldChemblTarget(
                target_chembl_id=target.chembl_id,
                target_type=target.target_type,
                preferred_name=target.preferred_name,
                accession=accession,
                organism=target.organism,
                tax_id=target.ncbi_taxa_id,
            )
    else:
        # Old API: chembl.chembl_targets()
        logger.debug("Using old pypath API: chembl.chembl_targets()")
        return chembl.chembl_targets()


def chembl_assays(max_pages=None):
    """
    Get ChEMBL assays.
    
    Returns a generator/list of assay records adapted to old API format.
    """
    if chembl is None:
        raise ImportError("pypath.inputs.chembl not available")
    
    if _has_new_api():
        # New API: chembl.assay()
        logger.debug("Using new pypath API: chembl.assay()")
        for assay in chembl.assay(max_pages=max_pages):
            yield OldChemblAssay(
                assay_chembl_id=assay.assay_chembl_id,
                assay_type=assay.assay_type,
                confidence_score=assay.confidence_score,
                target_chembl_id=assay.target_chembl_id,
            )
    else:
        # Old API: chembl.chembl_assays()
        logger.debug("Using old pypath API: chembl.chembl_assays()")
        return chembl.chembl_assays()


def chembl_documents(max_pages=None):
    """
    Get ChEMBL documents.
    
    Returns a dict mapping document_chembl_id -> pubmed_id.
    """
    if chembl is None:
        raise ImportError("pypath.inputs.chembl not available")
    
    if _has_new_api():
        # New API: chembl.document()
        logger.debug("Using new pypath API: chembl.document()")
        result = {}
        for doc in chembl.document(max_pages=max_pages):
            if doc.document_chembl_id and doc.pubmed_id:
                result[doc.document_chembl_id] = doc.pubmed_id
        return result
    else:
        # Old API: chembl.chembl_documents()
        logger.debug("Using old pypath API: chembl.chembl_documents()")
        return chembl.chembl_documents()


def chembl_mechanisms(max_pages=None):
    """
    Get ChEMBL mechanisms.
    
    Returns a generator/list of mechanism records adapted to old API format.
    """
    if chembl is None:
        raise ImportError("pypath.inputs.chembl not available")
    
    if _has_new_api():
        # New API: chembl.mechanism()
        logger.debug("Using new pypath API: chembl.mechanism()")
        for mech in chembl.mechanism(max_pages=max_pages):
            yield OldChemblMechanism(
                chembl=mech.molecule_chembl_id,
                action_type=mech.action_type,
                direct_interaction=mech.direct_interaction,
                mechanism_of_action=None,  # Not available in new API
                disease_efficacy=None,      # Not available in new API
                target_chembl_id=mech.target_chembl_id,
            )
    else:
        # Old API: chembl.chembl_mechanisms()
        logger.debug("Using old pypath API: chembl.chembl_mechanisms()")
        return chembl.chembl_mechanisms()


def chembl_drug_indications(max_pages=None):
    """
    Get ChEMBL drug indications.
    
    Returns a generator/list of indication records.
    """
    if chembl is None:
        raise ImportError("pypath.inputs.chembl not available")
    
    if _has_new_api():
        # New API: chembl.indication()
        logger.debug("Using new pypath API: chembl.indication()")
        return chembl.indication(max_pages=max_pages)
    else:
        # Old API: chembl.chembl_drug_indications()
        logger.debug("Using old pypath API: chembl.chembl_drug_indications()")
        return chembl.chembl_drug_indications()
