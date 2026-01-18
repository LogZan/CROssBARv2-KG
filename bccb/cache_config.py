import os
import tempfile
from pypath.share import settings

# Use workspace directory for all pypath data
PYPATH_CACHE_DIR = "/GenSIvePFS/users/data/pypath_cache"
PYPATH_PICKLE_DIR = "/GenSIvePFS/users/data/pypath_cache/pickles"
PYPATH_TEMP_DIR = "/GenSIvePFS/users/data/pypath_temp"

def setup_pypath_cache():
    """Configure pypath to use custom cache directory on the large filesystem."""
    # Create directories
    os.makedirs(PYPATH_CACHE_DIR, exist_ok=True)
    os.makedirs(PYPATH_PICKLE_DIR, exist_ok=True)
    os.makedirs(PYPATH_TEMP_DIR, exist_ok=True)
    
    # Set pypath settings
    settings.setup(cachedir=PYPATH_CACHE_DIR)
    settings.setup(pickle_dir=PYPATH_PICKLE_DIR)
    
    # Also set TMPDIR environment variable for temp files
    os.environ['TMPDIR'] = PYPATH_TEMP_DIR
    os.environ['TEMP'] = PYPATH_TEMP_DIR
    os.environ['TMP'] = PYPATH_TEMP_DIR
    
    # Reinitialize tempfile to use new directory
    tempfile.tempdir = PYPATH_TEMP_DIR
