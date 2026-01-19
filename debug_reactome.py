
import sys
import os

# Set pypath offline to minimalize network calls during init
os.environ['PYPATH_OFFLINE'] = '1'

# Add project root to path
sys.path.append(os.getcwd())

# Patch RaMP again 
try:
    from pypath.inputs.ramp import _rest
    import json
    def patched_ramp_id_types(entity_type=None):
        return set()
    _rest.ramp_id_types = patched_ramp_id_types
except ImportError:
    pass

from pypath.inputs import unichem
from pypath.share import settings

# Enable online for this test
settings.settings.offline = False

print("Testing unichem.unichem_mapping('DrugBank', 'ChEBI')...")
try:
    mapping = unichem.unichem_mapping("DrugBank", "ChEBI")
    print(f"Success with Case Sensitive. First 5: {list(mapping.items())[:5]}")
except Exception as e:
    print(f"Error: {e}")

print("Testing unichem.unichem_mapping(2, 7)...")
try:
    mapping = unichem.unichem_mapping(2, 7)
    print(f"Success with IDs. First 5: {list(mapping.items())[:5]}")
except Exception as e:
    print(f"Error: {e}")

print("\nListing available UniChem sources if possible...")
try:
    sources = unichem.unichem_sources()
    # Pypath unichem_sources returns {id: name}
    print("UniChem Sources:")
    for k, v in sources.items():
        print(f"{k}: {v}")
except Exception as e:
    print(f"Could not list sources: {e}")
