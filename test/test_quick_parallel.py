#!/usr/bin/env python3
"""
Quick test to verify parallel processing works with 3 organisms.
"""

import sys
import os
from time import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from bccb.ppi_adapter import PPI

print("="*70)
print("QUICK PARALLEL TEST - Single Organism (Human)")
print("="*70)
print()

# Test with a single organism that has lots of data: human
HUMAN = 9606

# Test 1: Serial
print(">>> Test 1: Serial (n_workers=1)")
ppi1 = PPI(organism=HUMAN, test_mode=False, export_csv=False)

t0 = time()
ppi1.download_string_data(cache=True, n_workers=1)
t_serial = time() - t0
n_serial = len(ppi1.string_ints)

print(f"   Time: {t_serial:.2f}s")
print(f"   Interactions: {n_serial:,}\n")

# Test 2: Parallel with 4 workers
print(">>> Test 2: Parallel (n_workers=4)")
ppi2 = PPI(organism=HUMAN, test_mode=False, export_csv=False)

t0 = time()
ppi2.download_string_data(cache=True, n_workers=4)
t_parallel = time() - t0
n_parallel = len(ppi2.string_ints)

print(f"   Time: {t_parallel:.2f}s")
print(f"   Interactions: {n_parallel:,}\n")

# Results
print("="*70)
print("RESULTS")
print("="*70)

speedup = t_serial / t_parallel if t_parallel > 0 else 0

print(f"\nSerial:   {t_serial:>8.2f}s  ({n_serial:,} interactions)")
print(f"Parallel: {t_parallel:>8.2f}s  ({n_parallel:,} interactions)")
print(f"\nSpeedup:    {speedup:.2f}x")

if n_serial == n_parallel:
    print(f"\n✅ PASS: Results match ({n_serial:,} interactions)")
else:
    print(f"\n❌ FAIL: Results differ (serial: {n_serial:,}, parallel: {n_parallel:,})")

# Note about single organism
print(f"\n💡 NOTE: With a single organism (n_workers={4}), there's no parallelization benefit")
print("   since each organism is processed sequentially.")
print("   Parallel processing benefits come when processing MULTIPLE organisms.")

print("\n" + "="*70)

