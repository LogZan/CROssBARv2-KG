
import sys
import os

# Add project root to path
sys.path.append(os.getcwd())

from pypath.inputs import uniprot
from pypath.share import settings

# Configure settings if needed (cache, etc.)
settings.setup(progress_bars=True)

print("Testing uniprot.uniprot_data('gene_names', '*', True)...")
try:
    gene_names = uniprot.uniprot_data("gene_names", "*", True)
    print(f"Type: {type(gene_names)}")
    if isinstance(gene_names, list):
        print(f"First 5 elements: {gene_names[:5]}")
    elif isinstance(gene_names, dict):
        print(f"First 5 keys: {list(gene_names.keys())[:5]}")
except Exception as e:
    print(f"Error: {e}")

print("\nTesting uniprot.uniprot_data('xref_string', '*', True)...")
try:
    xref_string = uniprot.uniprot_data("xref_string", "*", True)
    print(f"Type: {type(xref_string)}")
    if isinstance(xref_string, list):
        print(f"First 5 elements: {xref_string[:5]}")
except Exception as e:
    print(f"Error: {e}")
