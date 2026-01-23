#!/usr/bin/env python3
"""
Profile PPI adapter to identify the exact bottleneck causing 15.82s/it.
This test will measure each operation separately to find the slowest part.
"""

import sys
import time
import pandas as pd
from pathlib import Path
from bioregistry import normalize_curie

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

# Setup pypath cache before any pypath imports
from bccb import cache_config
cache_config.setup_pypath_cache()

from bccb.ppi_adapter import PPI

print("\n" + "="*80)
print("PPI BOTTLENECK PROFILING")
print("="*80)

# Configuration
CACHE = True
TEST_MODE = True  # Use test mode for faster profiling
ORGANISM = 9606  # Human

print("\nStep 1: Initialize PPI adapter...")
start = time.time()
ppi_adapter = PPI(
    organism=ORGANISM,
    output_dir=None,
    export_csv=False,
    test_mode=TEST_MODE
)
print(f"✓ Initialization took {time.time() - start:.2f}s")

print("\nStep 2: Download PPI data (using cache)...")
start = time.time()
ppi_adapter.download_ppi_data(cache=CACHE)
print(f"✓ Download took {time.time() - start:.2f}s")

print("\nStep 3: Process PPI data...")
start = time.time()
ppi_adapter.process_ppi_data()
print(f"✓ Processing took {time.time() - start:.2f}s")

print("\nStep 4: Merge all dataframes...")
start = time.time()
merged_df = ppi_adapter.merge_all()
print(f"✓ Merge took {time.time() - start:.2f}s")
print(f"  Total edges to process: {len(merged_df):,}")

# ========================================================================
# Now profile the get_ppi_edges function step by step
# ========================================================================

print("\n" + "="*80)
print("PROFILING get_ppi_edges FUNCTION")
print("="*80)

# Test with small sample first
SAMPLE_SIZE = 100
sample_df = merged_df.head(SAMPLE_SIZE)

print(f"\nTesting with {SAMPLE_SIZE} edges...")
print("-" * 80)

# Test 1: Just iterating
print("\nTest 1: Just iterating through dataframe...")
start = time.time()
count = 0
for _, row in sample_df.iterrows():
    count += 1
elapsed = time.time() - start
per_item = (elapsed / count) * 1000
print(f"  Time: {elapsed:.4f}s for {count} rows")
print(f"  Per item: {per_item:.2f}ms")
print(f"  Estimated for {len(merged_df):,} rows: {elapsed * len(merged_df) / count / 60:.2f} minutes")

# Test 2: Adding row.to_dict()
print("\nTest 2: Iterating + row.to_dict()...")
start = time.time()
count = 0
for _, row in sample_df.iterrows():
    _dict = row.to_dict()
    count += 1
elapsed = time.time() - start
per_item = (elapsed / count) * 1000
print(f"  Time: {elapsed:.4f}s for {count} rows")
print(f"  Per item: {per_item:.2f}ms")
print(f"  Estimated for {len(merged_df):,} rows: {elapsed * len(merged_df) / count / 60:.2f} minutes")

# Test 3: Adding normalize_curie (THE SUSPECTED BOTTLENECK)
print("\nTest 3: Iterating + row.to_dict() + normalize_curie (2x per row)...")
start = time.time()
count = 0
for _, row in sample_df.iterrows():
    _dict = row.to_dict()
    _source = normalize_curie("uniprot:" + str(row["uniprot_a"]))
    _target = normalize_curie("uniprot:" + str(row["uniprot_b"]))
    count += 1
elapsed = time.time() - start
per_item = (elapsed / count) * 1000
print(f"  Time: {elapsed:.4f}s for {count} rows")
print(f"  Per item: {per_item:.2f}ms")
print(f"  Estimated for {len(merged_df):,} rows: {elapsed * len(merged_df) / count / 60:.2f} minutes")

# Test 4: Without normalize_curie (baseline)
print("\nTest 4: Iterating + row.to_dict() + simple string concat (no normalize_curie)...")
start = time.time()
count = 0
for _, row in sample_df.iterrows():
    _dict = row.to_dict()
    _source = "uniprot:" + str(row["uniprot_a"])
    _target = "uniprot:" + str(row["uniprot_b"])
    count += 1
elapsed = time.time() - start
per_item = (elapsed / count) * 1000
print(f"  Time: {elapsed:.4f}s for {count} rows")
print(f"  Per item: {per_item:.2f}ms")
print(f"  Estimated for {len(merged_df):,} rows: {elapsed * len(merged_df) / count / 60:.2f} minutes")

# Test 5: Full edge creation (current implementation)
print("\nTest 5: Full edge creation (current implementation)...")
start = time.time()
edge_list = []
count = 0
label = "Protein_interacts_with_protein"
for _, row in sample_df.iterrows():
    _dict = row.to_dict()
    
    _source = normalize_curie("uniprot:" + str(row["uniprot_a"]))
    _target = normalize_curie("uniprot:" + str(row["uniprot_b"]))
    
    del _dict["uniprot_a"], _dict["uniprot_b"]
    
    _props = {}
    for k, v in _dict.items():
        if str(v) != "nan":
            if isinstance(v, str) and "|" in v:
                _props[str(k).replace(" ", "_").lower()] = v.replace(
                    "'", "^"
                ).split("|")
            else:
                _props[str(k).replace(" ", "_").lower()] = str(
                    v
                ).replace("'", "^")
    
    edge_list.append((None, _source, _target, label, _props))
    count += 1
elapsed = time.time() - start
per_item = (elapsed / count) * 1000
print(f"  Time: {elapsed:.4f}s for {count} rows")
print(f"  Per item: {per_item:.2f}ms")
print(f"  Estimated for {len(merged_df):,} rows: {elapsed * len(merged_df) / count / 60:.2f} minutes")

# Test 6: JUST normalize_curie calls in isolation
print("\nTest 6: JUST normalize_curie performance (isolated)...")
test_ids = sample_df["uniprot_a"].head(SAMPLE_SIZE).tolist()
start = time.time()
for protein_id in test_ids:
    _ = normalize_curie("uniprot:" + str(protein_id))
elapsed = time.time() - start
per_item = (elapsed / len(test_ids)) * 1000
print(f"  Time: {elapsed:.4f}s for {len(test_ids)} calls")
print(f"  Per call: {per_item:.2f}ms")
print(f"  For 2x{len(merged_df):,} calls: {elapsed * 2 * len(merged_df) / len(test_ids) / 60:.2f} minutes")

# Test 7: Caching normalize_curie results
print("\nTest 7: With normalize_curie caching...")
normalize_cache = {}
def cached_normalize_curie(curie_str):
    if curie_str not in normalize_cache:
        normalize_cache[curie_str] = normalize_curie(curie_str)
    return normalize_cache[curie_str]

start = time.time()
edge_list = []
count = 0
for _, row in sample_df.iterrows():
    _dict = row.to_dict()
    
    _source = cached_normalize_curie("uniprot:" + str(row["uniprot_a"]))
    _target = cached_normalize_curie("uniprot:" + str(row["uniprot_b"]))
    
    del _dict["uniprot_a"], _dict["uniprot_b"]
    
    _props = {}
    for k, v in _dict.items():
        if str(v) != "nan":
            if isinstance(v, str) and "|" in v:
                _props[str(k).replace(" ", "_").lower()] = v.replace(
                    "'", "^"
                ).split("|")
            else:
                _props[str(k).replace(" ", "_").lower()] = str(
                    v
                ).replace("'", "^")
    
    edge_list.append((None, _source, _target, label, _props))
    count += 1
elapsed = time.time() - start
per_item = (elapsed / count) * 1000
print(f"  Time: {elapsed:.4f}s for {count} rows")
print(f"  Per item: {per_item:.2f}ms")
print(f"  Cache size: {len(normalize_cache)}")
print(f"  Estimated for {len(merged_df):,} rows: {elapsed * len(merged_df) / count / 60:.2f} minutes")

print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print(f"Total edges in merged dataframe: {len(merged_df):,}")
print("\nBased on the profiling, the slowest operation is likely:")
print("  - If Test 3 is much slower than Test 4: normalize_curie is the bottleneck")
print("  - If Test 1-2 are slow: iterrows() is inefficient (should use apply/vectorize)")
print("  - Check the per-item times above to identify the culprit")
print("\nRecommendations:")
print("  1. If normalize_curie is slow: Add caching (Test 7 shows improvement)")
print("  2. If iterrows is slow: Use df.apply() or vectorized operations")
print("  3. Consider batch processing instead of row-by-row")
print("="*80)
