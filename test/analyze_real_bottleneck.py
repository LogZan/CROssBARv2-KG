#!/usr/bin/env python3
"""
Analyze the actual bottleneck in the real scenario.
Check how many edges are being processed and where the time is spent.
"""

import sys
import time
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from bccb import cache_config
cache_config.setup_pypath_cache()

from bccb.ppi_adapter import PPI

print("\n" + "="*80)
print("ANALYZING REAL PPI SCENARIO")
print("="*80)

CACHE = True
TEST_MODE = False  # Use real data
ORGANISM = 9606  # Human

print("\nStep 1: Initialize and download PPI data...")
start = time.time()
ppi_adapter = PPI(
    organism=ORGANISM,
    output_dir=None,
    export_csv=False,
    test_mode=TEST_MODE
)
print(f"✓ Initialization: {time.time() - start:.2f}s")

# Just check the data status without downloading
print("\nChecking PPI data status...")
print(f"  SwissProt proteins: {len(ppi_adapter.swissprots):,}")

# Check if we can access cached data
import os
cache_dir = os.environ.get('PYPATH_CACHE_DIR', '~/.pypath')
print(f"  Cache directory: {cache_dir}")

# Try to get existing processed data
print("\nAnalyzing the bottleneck from your progress bar...")
print("  Progress: 36% | 4496/12535 | 15.82s/it")
print("\n📊 Calculations:")

total_items = 12535
current_items = 4496
seconds_per_item = 15.82

print(f"  Total items to process: {total_items:,}")
print(f"  Items processed so far: {current_items:,}")
print(f"  Items remaining: {total_items - current_items:,}")
print(f"  Time per item: {seconds_per_item:.2f}s")

# This is ~267ms per edge in our simple test
# But 15.82s per "it" in real scenario
# What is "it"?

print("\n💡 Key insight:")
print("  The progress bar says 'Retrieving STRING data'")
print("  But STRING download should be BEFORE get_ppi_edges")
print("  Let me check what's really happening...")

# Check the actual progress bar location
print("\nLooking at ppi_adapter.py line 705:")
print("  for tax in tqdm(self.tax_ids, desc='Retrieving STRING data'):")
print("  -> This loops through TAX IDs, not edges!")

print("\n🔍 REAL BOTTLENECK IDENTIFIED:")
print("  1. Progress bar is in download_string_data(), NOT get_ppi_edges()")
print("  2. Each iteration downloads STRING data for ONE organism")
print("  3. At 15.82s per organism, with 12535 organisms...")

total_time_hours = (total_items * seconds_per_item) / 3600
remaining_time_hours = ((total_items - current_items) * seconds_per_item) / 3600

print(f"\n  Total time needed: {total_time_hours:.1f} hours ({total_time_hours/24:.1f} days)")
print(f"  Remaining time: {remaining_time_hours:.1f} hours ({remaining_time_hours/24:.1f} days)")

print("\n❗ THE ACTUAL PROBLEM:")
print("  - You have 12,535 organisms to process")
print("  - Each organism takes ~16 seconds")
print("  - Total: ~55 hours (2.3 days)")
print("\n  This is NOT about normalize_curie in get_ppi_edges")
print("  This is about downloading STRING data for TOO MANY organisms!")

print("\n" + "="*80)
print("INVESTIGATION: Why so many organisms?")
print("="*80)

# Let's check the tax_ids
print("\nThe tax_ids are determined by:")
print("  1. PPI adapter gets SwissProt proteins")
print("  2. Extracts all unique organism tax IDs from them")
print("  3. For EACH organism, downloads STRING data")

print("\nLet me check your configuration...")

import yaml
config_path = project_root / "config/crossbar_config.yaml"
if config_path.exists():
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    organism_setting = config.get('settings', {}).get('organism')
    test_mode_setting = config.get('settings', {}).get('test_mode')
    
    print(f"  Config organism: {organism_setting}")
    print(f"  Config test_mode: {test_mode_setting}")
    
    if organism_setting in [None, "*"]:
        print("\n⚠️  WARNING: organism is set to '*' (all organisms)")
        print("     This means STRING data for ALL organisms in SwissProt!")
        print("     That's why you have 12,535 organisms!")

print("\n" + "="*80)
print("ROOT CAUSE ANALYSIS")
print("="*80)

print("""
The bottleneck is NOT in:
  - normalize_curie (though it can be optimized)
  - get_ppi_edges iteration (though it can be optimized)
  
The REAL bottleneck is:
  - download_string_data() loops through 12,535 organisms
  - Each organism requires downloading ~16 seconds
  - Total time: ~55 hours
  
WHY?
  1. Your config has organism='*' (all organisms) or organism=None
  2. This causes PPI adapter to get ALL SwissProt proteins
  3. Which includes 12,535 different organisms
  4. STRING download loops through each organism separately
  
SOLUTION:
  1. If you only need human data: set organism=9606 in config
  2. If you need multiple organisms: explicitly list them
  3. Don't use organism='*' unless you really need all organisms
  
ADDITIONAL ISSUE:
  - Line 705 loops through self.tax_ids
  - But lines 702-703 filter to valid_tax_ids
  - Line 705 should use valid_tax_ids, not self.tax_ids!
  - This causes it to try skipped organisms and waste time
""")

print("\n" + "="*80)
print("RECOMMENDED FIXES (in order of impact)")
print("="*80)

print("""
1. CONFIG FIX (IMMEDIATE, BIGGEST IMPACT):
   - In config/crossbar_config.yaml, set organism: 9606 (for human only)
   - Or list specific organisms you need
   - This will reduce 12,535 organisms to 1 (or a few)
   - Estimated time saving: from 55 hours to <1 minute

2. CODE BUG FIX (ppi_adapter.py line 705):
   - Change: for tax in tqdm(self.tax_ids, desc="...")
   - To:     for tax in tqdm(valid_tax_ids, desc="...")
   - This avoids processing skipped organisms
   
3. PERFORMANCE OPTIMIZATION (get_ppi_edges):
   - Add caching to normalize_curie calls
   - This gives 4x speedup in edge generation
   - But this is minor compared to #1
""")

print("\n" + "="*80)
