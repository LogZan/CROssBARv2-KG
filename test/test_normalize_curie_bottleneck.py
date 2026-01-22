#!/usr/bin/env python3
"""
Test to verify normalize_curie is the bottleneck.
"""

import time
from bioregistry import normalize_curie

DATA_FILE = '/GenSIvePFS/users/data/InterPro/107.0/protein2ipr.dat'
TEST_LINES = 100000


def test_with_normalize_curie():
    """Test with normalize_curie (current implementation)"""
    print("\n" + "="*60)
    print("Test WITH normalize_curie")
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

            # WITH normalize_curie
            interpro_id = normalize_curie("interpro:" + entry_id)
            uniprot_id = normalize_curie("uniprot:" + protein_id)

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


def test_without_normalize_curie():
    """Test without normalize_curie (simple string concatenation)"""
    print("\n" + "="*60)
    print("Test WITHOUT normalize_curie")
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

            # WITHOUT normalize_curie - simple string concatenation
            interpro_id = "interpro:" + entry_id
            uniprot_id = "uniprot:" + protein_id

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


if __name__ == "__main__":
    print("\n" + "="*60)
    print("normalize_curie Bottleneck Test")
    print(f"Testing with {TEST_LINES:,} lines")
    print("="*60)

    speed_with = test_with_normalize_curie()
    speed_without = test_without_normalize_curie()

    print("\n" + "="*60)
    print("COMPARISON")
    print("="*60)
    print(f"WITH normalize_curie:    {speed_with:>10,.0f} lines/s")
    print(f"WITHOUT normalize_curie: {speed_without:>10,.0f} lines/s")
    print(f"\nSpeedup: {speed_without/speed_with:.1f}x faster")
    print(f"Time reduction: {(1 - speed_with/speed_without)*100:.1f}%")

    print("\n" + "="*60)
    print("ESTIMATED TOTAL TIME FOR 1.15 BILLION LINES")
    print("="*60)
    total_lines = 1_152_995_715
    time_with = total_lines / speed_with / 3600
    time_without = total_lines / speed_without / 3600
    print(f"WITH normalize_curie:    {time_with:.1f} hours")
    print(f"WITHOUT normalize_curie: {time_without:.1f} hours")
    print(f"Time saved: {time_with - time_without:.1f} hours")
