
import sys
import os

# Set pypath offline to minimalize network calls during init
os.environ['PYPATH_OFFLINE'] = '1'

import pypath.share.settings as settings

# Disable RaMP or other problematic modules if possible via settings
# (pypath doesn't have a simple disable flag for ramp, but we can try to mock it if needed)
# However, the error is at import time of pypath.inputs.uniprot -> pypath.utils.orthology -> pypath.utils.mapping

# Let's try to patch the ramp request behavior before importing pypath.utils.mapping
# But mapping is imported deep in the stack.

# Alternative: Mock the response for RaMP.
import json
from unittest.mock import MagicMock

# Attempt to patch pypath.inputs.ramp._rest
# We need to do this carefully.

print("Starting debug script...")

try:
    from pypath.inputs import uniprot
    print("Successfully imported pypath.inputs.uniprot")
    
    # Enable online again for the specific call we want to test
    settings.settings.offline = False

    print("Testing uniprot.uniprot_data('gene_names', '*', True)...")
    # We might need to handle the network error for this call if servers are flaky
    gene_names = uniprot.uniprot_data("gene_names", "*", True)
    print(f"Type: {type(gene_names)}")
    if isinstance(gene_names, list):
        print(f"First 5 elements: {gene_names[:5]}")
    elif isinstance(gene_names, dict):
        print(f"First 5 keys: {list(gene_names.keys())[:5]}")

    print("\nTesting uniprot.uniprot_data('xref_string', '*', True)...")
    xref_string = uniprot.uniprot_data("xref_string", "*", True)
    print(f"Type: {type(xref_string)}")
    if isinstance(xref_string, list):
        print(f"First 5 elements: {xref_string[:5]}")

except Exception as e:
    print(f"Caught exception: {e}")
    import traceback
    traceback.print_exc()
