#!/usr/bin/env python3
"""
Performance test for InterPro adapter to identify bottlenecks.
"""

import sys
import time
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from bccb import cache_config
cache_config.setup_pypath_cache()

from bccb.interpro_adapter import InterPro

DATA_FILE = '/GenSIvePFS/users/data/InterPro/107.0/protein2ipr.dat'
TEST_LINES = 100000  # Test with 100k lines


def test_1_raw_file_read():
    """Test 1: Pure file reading speed (baseline)"""
    print("\n" + "="*60)
    print("Test 1: Raw file reading speed")
    print("="*60)

    start = time.time()
    count = 0

    with open(DATA_FILE, 'r', encoding='utf-8') as f:
        for line in f:
            count += 1
            if count >= TEST_LINES:
                break

    elapsed = time.time() - start
    speed = count / elapsed

    print(f"Lines read: {count:,}")
    print(f"Time: {elapsed:.2f}s")
    print(f"Speed: {speed:,.0f} lines/s")
    return speed


def test_2_file_read_and_parse():
    """Test 2: File reading + basic parsing"""
    print("\n" + "="*60)
    print("Test 2: File reading + parsing")
    print("="*60)

    start = time.time()
    count = 0

    with open(DATA_FILE, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue

            protein_id, entry_id, entry_name, _, start_pos, end_pos = parts[:6]
            protein_id = protein_id.upper()

            count += 1
            if count >= TEST_LINES:
                break

    elapsed = time.time() - start
    speed = count / elapsed

    print(f"Lines processed: {count:,}")
    print(f"Time: {elapsed:.2f}s")
    print(f"Speed: {speed:,.0f} lines/s")
    return speed


def test_3_with_prefix_and_props():
    """Test 3: File reading + parsing + prefix + props"""
    print("\n" + "="*60)
    print("Test 3: File reading + parsing + prefix + props")
    print("="*60)

    from bioregistry import normalize_curie

    start = time.time()
    count = 0

    with open(DATA_FILE, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue

            protein_id, entry_id, entry_name, _, start_pos, end_pos = parts[:6]
            protein_id = protein_id.upper()

            # Add prefix
            interpro_id = normalize_curie("interpro:" + entry_id)
            uniprot_id = normalize_curie("uniprot:" + protein_id)

            # Create props
            props = {}
            props["start"] = int(start_pos) if start_pos.isdigit() else start_pos
            props["end"] = int(end_pos) if end_pos.isdigit() else end_pos

            count += 1
            if count >= TEST_LINES:
                break

    elapsed = time.time() - start
    speed = count / elapsed

    print(f"Lines processed: {count:,}")
    print(f"Time: {elapsed:.2f}s")
    print(f"Speed: {speed:,.0f} lines/s")
    return speed


def test_4_with_generator():
    """Test 4: Full processing with generator (no BioCypher)"""
    print("\n" + "="*60)
    print("Test 4: Full processing with generator")
    print("="*60)

    interpro_adapter = InterPro(
        organism="*",
        test_mode=False,
        data_dir='/GenSIvePFS/users/data/InterPro/107.0'
    )

    interpro_adapter.download_interpro_data(cache=True)

    start = time.time()
    count = 0

    for edge in interpro_adapter.get_interpro_edges():
        count += 1
        if count >= TEST_LINES:
            break

    elapsed = time.time() - start
    speed = count / elapsed

    print(f"Edges generated: {count:,}")
    print(f"Time: {elapsed:.2f}s")
    print(f"Speed: {speed:,.0f} lines/s")
    return speed


if __name__ == "__main__":
    print("\n" + "="*60)
    print("InterPro Performance Test")
    print(f"Testing with {TEST_LINES:,} lines")
    print("="*60)

    results = {}

    results['raw_read'] = test_1_raw_file_read()
    results['parse'] = test_2_file_read_and_parse()
    results['prefix_props'] = test_3_with_prefix_and_props()
    results['generator'] = test_4_with_generator()

    print("\n" + "="*60)
    print("PERFORMANCE SUMMARY")
    print("="*60)
    print(f"1. Raw file read:        {results['raw_read']:>10,.0f} lines/s (baseline)")
    print(f"2. + Parsing:            {results['parse']:>10,.0f} lines/s ({results['parse']/results['raw_read']*100:.1f}% of baseline)")
    print(f"3. + Prefix & props:     {results['prefix_props']:>10,.0f} lines/s ({results['prefix_props']/results['raw_read']*100:.1f}% of baseline)")
    print(f"4. + Generator:          {results['generator']:>10,.0f} lines/s ({results['generator']/results['raw_read']*100:.1f}% of baseline)")

    print("\n" + "="*60)
    print("BOTTLENECK ANALYSIS")
    print("="*60)

    overhead_parse = (1 - results['parse']/results['raw_read']) * 100
    overhead_prefix = (1 - results['prefix_props']/results['parse']) * 100
    overhead_generator = (1 - results['generator']/results['prefix_props']) * 100

    print(f"Parsing overhead:        {overhead_parse:>6.1f}%")
    print(f"Prefix/props overhead:   {overhead_prefix:>6.1f}%")
    print(f"Generator overhead:      {overhead_generator:>6.1f}%")

    print("\nNote: BioCypher write overhead not tested here.")
    print("If production speed is much slower, BioCypher is the bottleneck.")
