"""
PyPath Cache Configuration

Configures pypath to use categorized cache directories for better organization.
Each adapter gets its own subdirectory for clean data separation.
"""
import os
import tempfile
from pypath.share import settings

# Monkey-patch pypath_common to use logs subdirectory
try:
    import pypath_common._logger as _logger_module
    _original_new_logger = _logger_module.new_logger

    def _patched_new_logger(name, settings_obj, logdir=None, **kwargs):
        # If logdir not specified, use logs/{name}_log instead of {name}_log
        if logdir is None:
            logdir = f'logs/{name}_log'
        return _original_new_logger(name, settings_obj, logdir=logdir, **kwargs)

    _logger_module.new_logger = _patched_new_logger
except Exception:
    pass  # If patching fails, continue with default behavior

# Base cache directory
PYPATH_CACHE_BASE = "/GenSIvePFS/users/data/pypath_cache"
PYPATH_PICKLE_DIR = "/GenSIvePFS/users/data/pypath_cache/pickles"
PYPATH_TEMP_DIR = "/GenSIvePFS/users/data/pypath_temp"

# Adapter-specific cache directories
CACHE_DIRS = {
    # PPI sources
    'ppi': f'{PYPATH_CACHE_BASE}/ppi',
    'string': f'{PYPATH_CACHE_BASE}/string',
    'intact': f'{PYPATH_CACHE_BASE}/intact',
    'biogrid': f'{PYPATH_CACHE_BASE}/biogrid',
    
    # Protein/Gene databases  
    'uniprot': f'{PYPATH_CACHE_BASE}/uniprot',
    'interpro': f'{PYPATH_CACHE_BASE}/interpro',
    
    # Ontologies
    'go': f'{PYPATH_CACHE_BASE}/go',
    
    # Drug/Compound databases
    'drug': f'{PYPATH_CACHE_BASE}/drug',
    'chembl': f'{PYPATH_CACHE_BASE}/chembl',
    'drugbank': f'{PYPATH_CACHE_BASE}/drugbank',
    'compound': f'{PYPATH_CACHE_BASE}/compound',
    'pharos': f'{PYPATH_CACHE_BASE}/pharos',
    
    # Disease/Phenotype
    'disease': f'{PYPATH_CACHE_BASE}/disease',
    'phenotype': f'{PYPATH_CACHE_BASE}/phenotype',
    'hpo': f'{PYPATH_CACHE_BASE}/hpo',
    
    # Pathway databases
    'pathway': f'{PYPATH_CACHE_BASE}/pathway',
    'reactome': f'{PYPATH_CACHE_BASE}/reactome',
    'kegg': f'{PYPATH_CACHE_BASE}/kegg',
    
    # Other databases
    'side_effect': f'{PYPATH_CACHE_BASE}/side_effect',
    'ec': f'{PYPATH_CACHE_BASE}/ec',
    'orthology': f'{PYPATH_CACHE_BASE}/orthology',
    'tfgene': f'{PYPATH_CACHE_BASE}/tfgene',
    
    # Fallback
    'other': f'{PYPATH_CACHE_BASE}/other',
}

# Legacy alias for backward compatibility
PYPATH_CACHE_DIR = PYPATH_CACHE_BASE


def get_cache_dir(category: str) -> str:
    """
    Get the cache directory path for a specific data category.
    
    Args:
        category: Data category name (e.g., 'string', 'uniprot', 'chembl')
        
    Returns:
        Path to the category-specific cache directory
    """
    category = category.lower()
    if category in CACHE_DIRS:
        path = CACHE_DIRS[category]
    else:
        path = CACHE_DIRS['other']
    
    # Ensure directory exists
    os.makedirs(path, exist_ok=True)
    return path


def set_adapter_cache(adapter_name: str) -> str:
    """
    Switch pypath's cache directory to the adapter-specific location.
    
    Call this at the beginning of each adapter's download method to ensure
    data is saved to the correct subdirectory.
    
    Args:
        adapter_name: Name of the adapter (e.g., 'ppi', 'drug', 'go')
        
    Returns:
        Path to the cache directory that was set
    """
    cache_dir = get_cache_dir(adapter_name)
    settings.setup(cachedir=cache_dir)
    return cache_dir


def setup_pypath_cache():
    """Configure pypath to use custom cache directory on the large filesystem."""
    # Create base and pickle directories
    os.makedirs(PYPATH_CACHE_BASE, exist_ok=True)
    os.makedirs(PYPATH_PICKLE_DIR, exist_ok=True)
    os.makedirs(PYPATH_TEMP_DIR, exist_ok=True)
    os.makedirs('logs', exist_ok=True)

    # Create all category directories
    for category_dir in CACHE_DIRS.values():
        os.makedirs(category_dir, exist_ok=True)

    # Set pypath settings (default to 'other' dir for uncategorized files, adapters will override)
    settings.setup(cachedir=CACHE_DIRS['other'])
    settings.setup(pickle_dir=PYPATH_PICKLE_DIR)

    # Set environment variables for temp files and log directories
    os.environ['TMPDIR'] = PYPATH_TEMP_DIR
    os.environ['TEMP'] = PYPATH_TEMP_DIR
    os.environ['TMP'] = PYPATH_TEMP_DIR
    os.environ['PYPATH_LOG'] = 'logs'

    # Reinitialize tempfile to use new directory
    tempfile.tempdir = PYPATH_TEMP_DIR
