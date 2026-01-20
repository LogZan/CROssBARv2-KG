#!/usr/bin/env python3
"""
Debug script to test UniProt-STRING cross-reference retrieval.
"""

import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

# Setup pypath cache
from bccb import cache_config
cache_config.setup_pypath_cache()

from pypath.inputs import uniprot

print("=" * 60)
print("Testing UniProt-STRING cross-reference retrieval")
print("=" * 60)

# Test 1: Try with keyword arguments
print("\nTest 1: fields='xref_string', organism='*', reviewed=True")
try:
    result = uniprot.uniprot_data(fields="xref_string", organism="*", reviewed=True)
    print(f"Result type: {type(result)}")
    if isinstance(result, dict):
        print(f"Result length: {len(result)}")
        if result:
            # Show first 3 entries
            for i, (k, v) in enumerate(list(result.items())[:3]):
                print(f"  {k}: {v}")
    elif isinstance(result, list):
        print(f"Result is a list with {len(result)} items")
        if result:
            print(f"  First items: {result[:3]}")
    else:
        print(f"Result: {result}")
except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()

# Test 2: Try with organism=9606 (human only)
print("\nTest 2: fields='xref_string', organism=9606, reviewed=True")
try:
    result = uniprot.uniprot_data(fields="xref_string", organism=9606, reviewed=True)
    print(f"Result type: {type(result)}")
    if isinstance(result, dict):
        print(f"Result length: {len(result)}")
        if result:
            for i, (k, v) in enumerate(list(result.items())[:3]):
                print(f"  {k}: {v}")
    elif isinstance(result, list):
        print(f"Result is a list with {len(result)} items")
    else:
        print(f"Result: {result}")
except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()

# Test 2.5: Try with organism=None (all organisms?)
print("\nTest 2.5: fields='xref_string', organism=None, reviewed=True")
try:
    result = uniprot.uniprot_data(fields="xref_string", organism=None, reviewed=True)
    print(f"Result type: {type(result)}")
    if isinstance(result, dict):
        print(f"Result length: {len(result)}")
        if result:
            for i, (k, v) in enumerate(list(result.items())[:3]):
                print(f"  {k}: {v}")
    elif isinstance(result, list):
        print(f"Result is a list with {len(result)} items")
    else:
        print(f"Result: {result}")
except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()

# Test 3: Check available fields
print("\nTest 3: Checking if 'xref_string' is a valid field")
try:
    from pypath.inputs.uniprot import UniprotQuery
    query = UniprotQuery()

    # Check if xref_string is in field synonyms
    if hasattr(query, '_FIELD_SYNONYMS'):
        if 'xref_string' in query._FIELD_SYNONYMS:
            print(f"  'xref_string' is a valid field")
            print(f"  Maps to: {query._FIELD_SYNONYMS['xref_string']}")
        else:
            print(f"  'xref_string' is NOT in _FIELD_SYNONYMS")
            # Search for string-related fields
            string_fields = [k for k in query._FIELD_SYNONYMS.keys() if 'string' in k.lower()]
            print(f"  Available STRING-related fields: {string_fields}")
except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()

# Test 4: Try alternative field names
print("\nTest 4: Trying alternative field names")
for field_name in ['xref_string', 'database(STRING)', 'xref', 'cross-references']:
    print(f"\n  Trying field: '{field_name}'")
    try:
        result = uniprot.uniprot_data(fields=field_name, organism=9606, reviewed=True)
        print(f"    Result type: {type(result)}")
        if isinstance(result, dict) and result:
            print(f"    Result length: {len(result)}")
            first_key = list(result.keys())[0]
            print(f"    Sample: {first_key}: {result[first_key]}")
            break  # Found working field
    except Exception as e:
        print(f"    ERROR: {e}")

print("\n" + "=" * 60)
print("Debug complete")
print("=" * 60)
