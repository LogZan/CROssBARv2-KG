#!/usr/bin/env python3
"""
Simple profiling test to identify PPI bottleneck.
Focus on the core operations without full data processing.
"""

import time
import pandas as pd
from bioregistry import normalize_curie

print("\n" + "="*80)
print("SIMPLE PPI BOTTLENECK TEST")
print("="*80)

# Create a small test dataframe simulating PPI edges
print("\nCreating test dataframe with 10,000 edges...")
test_data = {
    'uniprot_a': ['P12345'] * 10000,
    'uniprot_b': ['Q67890'] * 10000,
    'source': ['IntAct'] * 10000,
    'score': [0.9] * 10000,
    'pmid': ['123456'] * 10000,
}
df = pd.DataFrame(test_data)
print(f"✓ Test dataframe created: {len(df)} rows")

# Simulate what happens in get_ppi_edges
print("\n" + "="*80)
print("TEST 1: Current implementation (iterrows + normalize_curie)")
print("="*80)

start = time.time()
edge_list = []
label = "Protein_interacts_with_protein"

for _, row in df.iterrows():
    _dict = row.to_dict()
    
    # This is the suspected bottleneck
    _source = normalize_curie("uniprot:" + str(row["uniprot_a"]))
    _target = normalize_curie("uniprot:" + str(row["uniprot_b"]))
    
    del _dict["uniprot_a"], _dict["uniprot_b"]
    
    _props = {}
    for k, v in _dict.items():
        if str(v) != "nan":
            _props[k] = str(v)
    
    edge_list.append((None, _source, _target, label, _props))

elapsed = time.time() - start
per_item_ms = (elapsed / len(df)) * 1000

print(f"Time: {elapsed:.2f}s for {len(df)} edges")
print(f"Per edge: {per_item_ms:.2f}ms")
print(f"Speed: {len(df)/elapsed:.0f} edges/s")

# Estimate for full dataset
print(f"\n📊 Extrapolation for different dataset sizes:")
for size in [100_000, 500_000, 1_000_000, 5_000_000]:
    estimated_time = (elapsed / len(df)) * size
    print(f"  {size:>10,} edges: {estimated_time/60:>8.1f} minutes ({estimated_time/3600:>6.2f} hours)")

# ========================================================================
print("\n" + "="*80)
print("TEST 2: With caching (cache normalize_curie results)")
print("="*80)

normalize_cache = {}
def cached_normalize(curie_str):
    if curie_str not in normalize_cache:
        normalize_cache[curie_str] = normalize_curie(curie_str)
    return normalize_cache[curie_str]

start = time.time()
edge_list = []

for _, row in df.iterrows():
    _dict = row.to_dict()
    
    _source = cached_normalize("uniprot:" + str(row["uniprot_a"]))
    _target = cached_normalize("uniprot:" + str(row["uniprot_b"]))
    
    del _dict["uniprot_a"], _dict["uniprot_b"]
    
    _props = {}
    for k, v in _dict.items():
        if str(v) != "nan":
            _props[k] = str(v)
    
    edge_list.append((None, _source, _target, label, _props))

elapsed_cached = time.time() - start
per_item_ms_cached = (elapsed_cached / len(df)) * 1000

print(f"Time: {elapsed_cached:.2f}s for {len(df)} edges")
print(f"Per edge: {per_item_ms_cached:.2f}ms")
print(f"Speed: {len(df)/elapsed_cached:.0f} edges/s")
print(f"Cache size: {len(normalize_cache)} unique CURIEs")
print(f"Cache hit rate: {(1 - len(normalize_cache)/(2*len(df)))*100:.1f}%")

print(f"\n📊 Extrapolation for different dataset sizes:")
for size in [100_000, 500_000, 1_000_000, 5_000_000]:
    estimated_time = (elapsed_cached / len(df)) * size
    print(f"  {size:>10,} edges: {estimated_time/60:>8.1f} minutes ({estimated_time/3600:>6.2f} hours)")

# ========================================================================
print("\n" + "="*80)
print("TEST 3: Vectorized approach (no iterrows)")
print("="*80)

start = time.time()

# Pre-compute all normalized IDs using caching
df['norm_a'] = df['uniprot_a'].apply(lambda x: cached_normalize("uniprot:" + str(x)))
df['norm_b'] = df['uniprot_b'].apply(lambda x: cached_normalize("uniprot:" + str(x)))

# Build edges
edge_list = []
for _, row in df.iterrows():
    _source = row['norm_a']
    _target = row['norm_b']
    _props = {
        'source': str(row['source']),
        'score': str(row['score']),
        'pmid': str(row['pmid']),
    }
    edge_list.append((None, _source, _target, label, _props))

elapsed_vectorized = time.time() - start
per_item_ms_vectorized = (elapsed_vectorized / len(df)) * 1000

print(f"Time: {elapsed_vectorized:.2f}s for {len(df)} edges")
print(f"Per edge: {per_item_ms_vectorized:.2f}ms")
print(f"Speed: {len(df)/elapsed_vectorized:.0f} edges/s")

print(f"\n📊 Extrapolation for different dataset sizes:")
for size in [100_000, 500_000, 1_000_000, 5_000_000]:
    estimated_time = (elapsed_vectorized / len(df)) * size
    print(f"  {size:>10,} edges: {estimated_time/60:>8.1f} minutes ({estimated_time/3600:>6.2f} hours)")

# ========================================================================
print("\n" + "="*80)
print("TEST 4: Without normalize_curie (simple string concatenation)")
print("="*80)

start = time.time()
edge_list = []

for _, row in df.iterrows():
    _dict = row.to_dict()
    
    # No normalize_curie - just simple string concatenation
    _source = "uniprot:" + str(row["uniprot_a"])
    _target = "uniprot:" + str(row["uniprot_b"])
    
    del _dict["uniprot_a"], _dict["uniprot_b"]
    
    _props = {}
    for k, v in _dict.items():
        if str(v) != "nan":
            _props[k] = str(v)
    
    edge_list.append((None, _source, _target, label, _props))

elapsed_no_normalize = time.time() - start
per_item_ms_no_normalize = (elapsed_no_normalize / len(df)) * 1000

print(f"Time: {elapsed_no_normalize:.2f}s for {len(df)} edges")
print(f"Per edge: {per_item_ms_no_normalize:.2f}ms")
print(f"Speed: {len(df)/elapsed_no_normalize:.0f} edges/s")

print(f"\n📊 Extrapolation for different dataset sizes:")
for size in [100_000, 500_000, 1_000_000, 5_000_000]:
    estimated_time = (elapsed_no_normalize / len(df)) * size
    print(f"  {size:>10,} edges: {estimated_time/60:>8.1f} minutes ({estimated_time/3600:>6.2f} hours)")

# ========================================================================
print("\n" + "="*80)
print("SUMMARY COMPARISON")
print("="*80)

speedup_cached = elapsed / elapsed_cached
speedup_vectorized = elapsed / elapsed_vectorized
speedup_no_normalize = elapsed / elapsed_no_normalize

print(f"\nMethod                          Time      Per Edge    Speed        Speedup")
print(f"-" * 80)
print(f"1. Current (iterrows + norm)   {elapsed:6.2f}s   {per_item_ms:7.2f}ms   {len(df)/elapsed:6.0f}/s      1.00x")
print(f"2. With caching                {elapsed_cached:6.2f}s   {per_item_ms_cached:7.2f}ms   {len(df)/elapsed_cached:6.0f}/s    {speedup_cached:6.2f}x")
print(f"3. Vectorized + cache          {elapsed_vectorized:6.2f}s   {per_item_ms_vectorized:7.2f}ms   {len(df)/elapsed_vectorized:6.0f}/s    {speedup_vectorized:6.2f}x")
print(f"4. No normalize_curie          {elapsed_no_normalize:6.2f}s   {per_item_ms_no_normalize:7.2f}ms   {len(df)/elapsed_no_normalize:6.0f}/s    {speedup_no_normalize:6.2f}x")

print("\n" + "="*80)
print("FINDINGS")
print("="*80)
print("\n✓ Identified bottleneck:")
if speedup_cached > 2:
    print(f"  - normalize_curie is a MAJOR bottleneck ({speedup_cached:.1f}x speedup with caching)")
    print(f"  - Adding cache reduces time from {elapsed:.1f}s to {elapsed_cached:.1f}s")
elif speedup_no_normalize > 2:
    print(f"  - normalize_curie is the bottleneck ({speedup_no_normalize:.1f}x speedup without it)")
else:
    print(f"  - iterrows() itself is slow, consider full vectorization")

print("\n✓ Recommended solution:")
if speedup_cached > 1.5:
    print("  - Add caching to normalize_curie calls")
    print("  - This will handle duplicate protein IDs efficiently")
if speedup_vectorized > speedup_cached:
    print("  - Consider vectorized approach for even better performance")
else:
    print("  - Current row-by-row approach is acceptable with caching")

print("\n" + "="*80)
