#!/usr/bin/env python3
"""
Test organism parameter behavior across different pypath functions.
"""

import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from bccb import cache_config
cache_config.setup_pypath_cache()

from pypath.inputs import intact, biogrid, uniprot

print("=" * 60)
print("Testing organism parameter across pypath functions")
print("=" * 60)

# Test 1: IntAct with different organism values
print("\n1. IntAct - organism parameter")
for org_value, org_label in [(None, "None"), (9606, "9606"), ("*", '"*"')]:
    print(f"\n  Testing organism={org_label}:")
    try:
        result = intact.intact_interactions(
            miscore=0,
            organism=org_value,
            complex_expansion=True,
            only_proteins=True
        )
        print(f"    Result: {len(result)} interactions")
        if len(result) > 0:
            print(f"    Sample: {result[0]}")
    except Exception as e:
        print(f"    ERROR: {e}")

# Test 2: BioGRID with different organism values
print("\n2. BioGRID - organism parameter")
for org_value, org_label in [(None, "None"), (9606, "9606"), ("*", '"*"')]:
    print(f"\n  Testing organism={org_label}:")
    try:
        result = biogrid.biogrid_all_interactions(
            org_value, 9999999999, False
        )
        print(f"    Result: {len(result)} interactions")
    except Exception as e:
        print(f"    ERROR: {e}")

# Test 3: UniProt uniprot_data with different organism values
print("\n3. UniProt uniprot_data - organism parameter")
for org_value, org_label in [(None, "None"), (9606, "9606"), ("*", '"*"')]:
    print(f"\n  Testing organism={org_label}:")
    try:
        result = uniprot.uniprot_data(
            fields="xref_string",
            organism=org_value,
            reviewed=True
        )
        if isinstance(result, dict):
            print(f"    Result: dict with {len(result)} entries")
        elif isinstance(result, list):
            print(f"    Result: list with {len(result)} items")
        else:
            print(f"    Result: {type(result)}")
    except Exception as e:
        print(f"    ERROR: {e}")

print("\n" + "=" * 60)
print("Summary:")
print("  - Check which functions accept '*' vs None for all organisms")
print("  - Check if behavior is consistent across pypath functions")
print("=" * 60)
