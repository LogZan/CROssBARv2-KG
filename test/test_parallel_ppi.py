#!/usr/bin/env python3
"""
Test parallel processing of STRING data in PPI adapter.

This script tests:
1. Serial processing (n_workers=1)
2. Parallel processing with auto-detect (n_workers=None)
3. Parallel processing with fixed workers (n_workers=4, 8, etc.)
4. Speed comparison and validation

Usage:
    python test/test_parallel_ppi.py
"""

import sys
import os
from time import time
import logging

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from bccb.ppi_adapter import PPI

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def test_ppi_parallel(n_workers, test_name, num_organisms=10):
    """
    Test PPI adapter with specified number of workers.
    
    Args:
        n_workers: Number of workers (1=serial, None=auto, >1=parallel)
        test_name: Name for this test
        num_organisms: Number of organisms to process (default: 10)
    
    Returns:
        tuple: (elapsed_time, num_interactions)
    """
    logger.info(f"\n{'='*60}")
    logger.info(f"Test: {test_name} (n_workers={n_workers})")
    logger.info(f"{'='*60}")
    
    # Create PPI adapter - use specific organisms for consistent testing
    # Use test mode to limit to a few organisms for quick testing
    ppi = PPI(
        organism=None,  # Process multiple organisms
        test_mode=False,  # Don't use test mode so we can control organism list
        export_csv=False
    )
    
    t0 = time()
    
    # Download STRING data with specified workers
    ppi.download_string_data(cache=True, n_workers=n_workers)
    
    t1 = time()
    elapsed = t1 - t0
    
    num_interactions = len(ppi.string_ints)
    
    logger.info(f"\n{test_name} Results:")
    logger.info(f"  Time: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
    logger.info(f"  Interactions: {num_interactions:,}")
    logger.info(f"  Speed: {num_interactions/elapsed:.1f} interactions/sec")
    
    return elapsed, num_interactions


def run_comparison_tests():
    """
    Run comprehensive comparison tests.
    """
    print("\n" + "="*80)
    print("PPI ADAPTER PARALLEL PROCESSING PERFORMANCE TEST")
    print("="*80)
    
    results = {}
    
    # Test 1: Serial processing (baseline)
    print("\n>>> Running Test 1: Serial Processing (n_workers=1)")
    t_serial, n_serial = test_ppi_parallel(
        n_workers=1,
        test_name="Serial Processing"
    )
    results['serial'] = {'time': t_serial, 'interactions': n_serial}
    
    # Test 2: Auto-detect parallel
    print("\n>>> Running Test 2: Auto-detect Parallel (n_workers=None)")
    t_auto, n_auto = test_ppi_parallel(
        n_workers=None,
        test_name="Auto-detect Parallel"
    )
    results['auto'] = {'time': t_auto, 'interactions': n_auto}
    
    # Test 3: Fixed 4 workers
    print("\n>>> Running Test 3: Fixed 4 Workers (n_workers=4)")
    t_4w, n_4w = test_ppi_parallel(
        n_workers=4,
        test_name="4 Workers"
    )
    results['4_workers'] = {'time': t_4w, 'interactions': n_4w}
    
    # Test 4: Fixed 8 workers
    print("\n>>> Running Test 4: Fixed 8 Workers (n_workers=8)")
    t_8w, n_8w = test_ppi_parallel(
        n_workers=8,
        test_name="8 Workers"
    )
    results['8_workers'] = {'time': t_8w, 'interactions': n_8w}
    
    # Print summary
    print("\n" + "="*80)
    print("PERFORMANCE SUMMARY")
    print("="*80)
    
    print(f"\n{'Configuration':<25} {'Time (s)':<15} {'Speedup':<15} {'Interactions':<15}")
    print("-" * 80)
    
    baseline_time = results['serial']['time']
    
    for name, data in results.items():
        t = data['time']
        n = data['interactions']
        speedup = baseline_time / t if t > 0 else 0
        
        display_name = {
            'serial': 'Serial (baseline)',
            'auto': 'Auto-detect',
            '4_workers': '4 Workers',
            '8_workers': '8 Workers'
        }.get(name, name)
        
        print(f"{display_name:<25} {t:>10.2f} s    {speedup:>8.2f}x      {n:>12,}")
    
    # Validation
    print("\n" + "="*80)
    print("VALIDATION")
    print("="*80)
    
    base_interactions = results['serial']['interactions']
    all_same = all(data['interactions'] == base_interactions for data in results.values())
    
    if all_same:
        print(f"✅ PASS: All configurations produced same number of interactions ({base_interactions:,})")
    else:
        print("❌ FAIL: Different configurations produced different interaction counts!")
        for name, data in results.items():
            print(f"  {name}: {data['interactions']:,}")
    
    # Calculate efficiency
    print("\n" + "="*80)
    print("PARALLEL EFFICIENCY")
    print("="*80)
    
    for name, data in results.items():
        if name == 'serial':
            continue
        
        workers = {
            'auto': 'auto',
            '4_workers': 4,
            '8_workers': 8
        }.get(name, 1)
        
        speedup = baseline_time / data['time']
        
        if isinstance(workers, int):
            efficiency = (speedup / workers) * 100
            print(f"{name}: {speedup:.2f}x speedup / {workers} workers = {efficiency:.1f}% efficiency")
        else:
            print(f"{name}: {speedup:.2f}x speedup")
    
    print("\n" + "="*80)
    print("TEST COMPLETE")
    print("="*80)
    
    return results


def quick_test():
    """
    Quick test with just 5 organisms to verify functionality.
    """
    print("\n" + "="*80)
    print("QUICK FUNCTIONALITY TEST (5 organisms)")
    print("="*80)
    
    from bccb.ppi_adapter import PPI
    
    # Test with very limited organisms
    ppi = PPI(organism=None, test_mode=True, export_csv=False)
    
    print("\n>>> Testing Serial (n_workers=1)")
    t0 = time()
    ppi.download_string_data(cache=True, n_workers=1)
    t_serial = time() - t0
    n_serial = len(ppi.string_ints)
    print(f"Serial: {t_serial:.2f}s, {n_serial} interactions")
    
    # Reset
    ppi2 = PPI(organism=None, test_mode=True, export_csv=False)
    
    print("\n>>> Testing Parallel (n_workers=4)")
    t0 = time()
    ppi2.download_string_data(cache=True, n_workers=4)
    t_parallel = time() - t0
    n_parallel = len(ppi2.string_ints)
    print(f"Parallel: {t_parallel:.2f}s, {n_parallel} interactions")
    
    speedup = t_serial / t_parallel if t_parallel > 0 else 0
    print(f"\nSpeedup: {speedup:.2f}x")
    
    if n_serial == n_parallel:
        print(f"✅ PASS: Same results ({n_serial} interactions)")
    else:
        print(f"❌ FAIL: Different results (serial: {n_serial}, parallel: {n_parallel})")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Test PPI parallel processing')
    parser.add_argument('--quick', action='store_true', help='Run quick test with 1 organism')
    parser.add_argument('--full', action='store_true', help='Run full comparison test')
    
    args = parser.parse_args()
    
    if args.quick:
        quick_test()
    elif args.full:
        run_comparison_tests()
    else:
        # Default: run quick test
        print("Running quick test (use --full for comprehensive comparison)")
        quick_test()
