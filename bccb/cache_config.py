import os
from pypath.share import settings

PYPATH_CACHE_DIR = "/GenSIvePFS/users/data/pypath_cache"

def setup_pypath_cache():
    """Configure pypath to use custom cache directory."""
    os.makedirs(PYPATH_CACHE_DIR, exist_ok=True)
    settings.setup(cachedir=PYPATH_CACHE_DIR)
