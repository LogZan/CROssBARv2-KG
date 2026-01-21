#!/usr/bin/env python3
"""
Test UniProt organism parameter behavior for uniprot_data and _all_uniprots.
"""

import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from bccb import cache_config
cache_config.setup_pypath_cache()

from pypath.inputs import uniprot


def _len_or_unknown(value):
    if isinstance(value, (list, dict, tuple, set)):
        return len(value)
    return None


def _sample_items(value, limit=10):
    if isinstance(value, dict):
        return list(value.keys())[:limit]
    if isinstance(value, (list, tuple)):
        return list(value)[:limit]
    if isinstance(value, set):
        return list(sorted(value))[:limit]
    return []


def test_uniprot_data_organism_param():
    print("=" * 60)
    print("Test uniprot.uniprot_data organism parameter")
    print("=" * 60)

    results = {}
    for org_value, label in [("*", '"*"'), (None, "None")]:
        print(f"\nCalling uniprot_data(fields='gene_names', organism={label}, reviewed=True)")
        try:
            result = uniprot.uniprot_data(
                fields="gene_names",
                organism=org_value,
                reviewed=True,
            )
            results[label] = result
            length = _len_or_unknown(result)
            print(f"  Result type: {type(result).__name__}, size: {length}")
            sample = _sample_items(result, limit=10)
            if sample:
                print(f"  Sample (up to 10): {sample}")
        except Exception as exc:
            results[label] = exc
            print(f"  ERROR: {exc}")

    star = results.get('"*"')
    none = results.get("None")

    failed = False
    if isinstance(star, Exception):
        print("  FAIL: organism='*' raised an error for uniprot_data.")
        failed = True
    else:
        star_len = _len_or_unknown(star)
        if star_len is None or star_len > 0:
            print("  FAIL: organism='*' did not return empty for uniprot_data.")
            failed = True
        else:
            print("  PASS: organism='*' returned empty for uniprot_data.")

    if isinstance(none, Exception):
        print("  FAIL: organism=None raised an error for uniprot_data.")
        failed = True
    else:
        none_len = _len_or_unknown(none)
        if none_len is None or none_len == 0:
            print("  FAIL: organism=None did not return data for uniprot_data.")
            failed = True
        else:
            print("  PASS: organism=None returned data for uniprot_data.")

    return not failed


def test_all_uniprots_organism_param():
    print("=" * 60)
    print("Test uniprot._all_uniprots organism parameter")
    print("=" * 60)

    failed = False

    print("\nCalling _all_uniprots(organism='*', swissprot=True)")
    try:
        star = uniprot._all_uniprots("*", True)
        star_len = _len_or_unknown(star)
        print(f"  Result type: {type(star).__name__}, size: {star_len}")
        sample = _sample_items(star, limit=10)
        if sample:
            print(f"  Sample (up to 10): {sample}")
        if star_len is None or star_len == 0:
            print("  FAIL: organism='*' did not return data for _all_uniprots.")
            failed = True
        else:
            print("  PASS: organism='*' returned data for _all_uniprots.")
    except Exception as exc:
        print(f"  FAIL: organism='*' raised an error: {exc}")
        failed = True

    print("\nCalling _all_uniprots(organism=None, swissprot=True)")
    try:
        _ = uniprot._all_uniprots(None, True)
        print("  FAIL: organism=None did not raise an error for _all_uniprots.")
        failed = True
    except Exception as exc:
        print(f"  PASS: organism=None raised an error: {exc}")

    return not failed


if __name__ == "__main__":
    ok_data = test_uniprot_data_organism_param()
    ok_all = test_all_uniprots_organism_param()
    if not (ok_data and ok_all):
        sys.exit(1)
    print("\nAll checks passed.")
