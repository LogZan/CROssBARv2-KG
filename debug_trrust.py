
import sys
import os

# Set pypath offline to minimalize network calls during init
os.environ['PYPATH_OFFLINE'] = '1'

# Add project root to path
sys.path.append(os.getcwd())

# Patch RaMP again just in case valid pypath is loaded from env
try:
    from pypath.inputs.ramp import _rest
    import json
    def patched_ramp_id_types(entity_type=None):
        return set()
    _rest.ramp_id_types = patched_ramp_id_types
except ImportError:
    pass

from pypath.inputs import trrust
from pypath.share import settings

# Enable online for this test
settings.settings.offline = False

print("Testing trrust.trrust_human()...")
try:
    human_data = trrust.trrust_human()
    print(f"Type: {type(human_data)}")
    if human_data is None:
        print("human_data is None")
    else:
        print(f"Length: {len(list(human_data))}")
except Exception as e:
    print(f"Error: {e}")
