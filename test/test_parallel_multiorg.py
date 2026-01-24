#!/usr/bin/env python3
"""
Test parallel processing of STRING data with multiple organisms.

This script tests parallel vs serial processing with a controlled set of organisms.

Usage:
    python test/test_parallel_multiorg.py
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
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def test_with_organisms(organism_list, n_workers, test_name):
    """
    Test PPI processing with specific organisms.
    
    Args:
        organism_list: List of NCBI taxonomy IDs to process
        n_workers: Number of workers
        test_name: Name for this test
    
    Returns:
        tuple: (elapsed_time, num_interactions, num_organisms_processed)
    """
    logger.info(f"\n{'='*70}")
    logger.info(f"Test: {test_name}")
    logger.info(f"Organisms: {len(organism_list)}, Workers: {n_workers}")
    logger.info(f"{'='*70}")
    
    # Create PPI adapter
    ppi = PPI(
        organism=None,
        test_mode=False,
        export_csv=False
    )
    
    # Manually set the organism list
    ppi.tax_ids = organism_list
    ppi.organism = None  # Ensure we use the tax_ids list
    
    t0 = time()
    
    # Download STRING data
    ppi.download_string_data(cache=True, n_workers=n_workers)
    
    t1 = time()
    elapsed = t1 - t0
    
    num_interactions = len(ppi.string_ints)
    
    logger.info(f"Results:")
    logger.info(f"  Time: {elapsed:.2f}s ({elapsed/60:.2f} min)")
    logger.info(f"  Interactions: {num_interactions:,}")
    logger.info(f"  Speed: {len(organism_list)/elapsed:.2f} organisms/sec")
    if num_interactions > 0:
        logger.info(f"  Throughput: {num_interactions/elapsed:.1f} interactions/sec")
    
    return elapsed, num_interactions, len(organism_list)


def main():
    """
    Run comparison test with multiple organisms.
    """
    print("\n" + "="*80)
    print("PARALLEL PPI PROCESSING TEST - MULTIPLE ORGANISMS")
    print("="*80)
    
    # Select a diverse set of organisms with known STRING data
    # These are common model organisms that should have cached data
    test_organisms = [
        "9606",   # Homo sapiens (human)
        "10090",  # Mus musculus (mouse)
        "10116",  # Rattus norvegicus (rat)
        "7227",   # Drosophila melanogaster (fruit fly)
        "6239",   # Caenorhabditis elegans (worm)
        "7955",   # Danio rerio (zebrafish)
        "4932",   # Saccharomyces cerevisiae (yeast)
        "9031",   # Gallus gallus (chicken)
        "3702",   # Arabidopsis thaliana (plant)
        "511145", # Escherichia coli (bacteria)
    ]
    
    print(f"\nTesting with {len(test_organisms)} organisms:")
    for org in test_organisms:
        print(f"  - {org}")
    
    results = {}
    
    # Test 1: Serial processing
    print("\n>>> Test 1: Serial Processing (n_workers=1)")
    t_serial, n_serial, _ = test_with_organisms(
        test_organisms,
        n_workers=1,
        test_name="Serial (baseline)"
    )
    results['serial'] = {'time': t_serial, 'interactions': n_serial}
    
    # Test 2: Parallel with 4 workers
    print("\n>>> Test 2: Parallel Processing (n_workers=4)")
    t_4w, n_4w, _ = test_with_organisms(
        test_organisms,
        n_workers=4,
        test_name="4 Workers"
    )
    results['4_workers'] = {'time': t_4w, 'interactions': n_4w}
    
    # Test 3: Parallel with 8 workers
    print("\n>>> Test 3: Parallel Processing (n_workers=8)")
    t_8w, n_8w, _ = test_with_organisms(
        test_organisms,
        n_workers=8,
        test_name="8 Workers"
    )
    results['8_workers'] = {'time': t_8w, 'interactions': n_8w}
    
    # Test 4: Auto-detect
    print("\n>>> Test 4: Auto-detect (n_workers=None)")
    t_auto, n_auto, _ = test_with_organisms(
        test_organisms,
        n_workers=None,
        test_name="Auto-detect"
    )
    results['auto'] = {'time': t_auto, 'interactions': n_auto}
    
    # Print summary
    print("\n" + "="*80)
    print("PERFORMANCE SUMMARY")
    print("="*80)
    
    print(f"\n{'Configuration':<20} {'Time':<15} {'Speedup':<12} {'Efficiency':<12} {'Interactions'}")
    print("-" * 80)
    
    baseline_time = results['serial']['time']
    
    configs = [
        ('serial', 'Serial (baseline)', 1),
        ('4_workers', '4 Workers', 4),
        ('8_workers', '8 Workers', 8),
        ('auto', 'Auto-detect', None),
    ]
    
    for key, name, workers in configs:
        data = results[key]
        t = data['time']
        n = data['interactions']
        speedup = baseline_time / t if t > 0 else 0
        
        if workers and workers > 1:
            efficiency = (speedup / workers) * 100
            eff_str = f"{efficiency:.1f}%"
        else:
            eff_str = "-"
        
        print(f"{name:<20} {t:>8.2f}s      {speedup:>6.2f}x      {eff_str:>8}     {n:>12,}")
    
    # Validation
    print("\n" + "="*80)
    print("VALIDATION")
    print("="*80)
    
    base_interactions = results['serial']['interactions']
    
    print(f"\nInteraction counts:")
    for key, name, _ in configs:
        count = results[key]['interactions']
        match = "✅" if count == base_interactions else "❌"
        print(f"  {match} {name:<20}: {count:>12,}")
    
    all_match = all(results[k]['interactions'] == base_interactions for k, _, _ in configs)
    
    if all_match:
        print(f"\n✅ PASS: All configurations produced identical results ({base_interactions:,} interactions)")
    else:
        print(f"\n❌ WARNING: Some configurations produced different results!")
    
    # Recommendations
    print("\n" + "="*80)
    print("RECOMMENDATIONS")
    print("="*80)
    
    best_speedup = max((baseline_time / results[k]['time'], name) 
                       for k, name, w in configs if k != 'serial')
    
    speedup_val, best_config = best_speedup
    
    print(f"\nBest configuration: {best_config}")
    print(f"Speedup achieved: {speedup_val:.2f}x")
    print(f"Time saved: {baseline_time - results['4_workers']['time']:.1f}s "
          f"({(1 - results['4_workers']['time']/baseline_time)*100:.1f}%)")
    
    # Estimate full processing time
    total_organisms = 12535
    avg_time_per_org_serial = baseline_time / len(test_organisms)
    avg_time_per_org_parallel = results['4_workers']['time'] / len(test_organisms)
    
    est_serial_hours = (total_organisms * avg_time_per_org_serial) / 3600
    est_parallel_hours = (total_organisms * avg_time_per_org_parallel) / 3600
    
    print(f"\n📊 Estimated time for all {total_organisms} organisms:")
    print(f"  Serial:   {est_serial_hours:.1f} hours")
    print(f"  4-workers: {est_parallel_hours:.1f} hours")
    print(f"  Time saved: {est_serial_hours - est_parallel_hours:.1f} hours")
    
    print("\n" + "="*80)
    print("TEST COMPLETE")
    print("="*80)


if __name__ == "__main__":
    main()
