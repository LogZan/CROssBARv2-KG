#!/usr/bin/env python3
"""
Debug script to analyze memory usage in Drug and Compound adapters.
This helps identify which specific operations cause OOM.
"""

import sys
import gc
import os
import psutil
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

def get_memory_mb():
    """Get current memory usage in MB."""
    process = psutil.Process()
    return process.memory_info().rss / 1024 / 1024

def print_memory(label):
    """Print memory usage with label."""
    mem = get_memory_mb()
    print(f"[MEMORY] {label}: {mem:.1f} MB")
    return mem

def test_drug_adapter_step_by_step():
    """Test Drug adapter step by step to find memory bottleneck."""
    print("=" * 80)
    print("Testing Drug Adapter Memory Usage Step by Step")
    print("=" * 80)
    
    initial_mem = print_memory("Initial")
    
    # Step 1: Import
    print("\n--- Step 1: Import modules ---")
    from bccb.drug_adapter import Drug, DrugNodeField, DrugEdgeType
    print_memory("After import")
    
    # Step 2: Initialize (without loading SwissProt)
    print("\n--- Step 2: Initialize Drug adapter ---")
    drug_adapter = Drug(
        drugbank_user="zengchuanlong23@mails.ucas.ac.cn",
        drugbank_passwd="iHDTbZ3vpFaWzC",
        test_mode=True,
        low_memory_mode=True
    )
    print_memory("After Drug.__init__")
    print(f"  _swissprots is None (lazy): {drug_adapter._swissprots is None}")
    
    # Step 3: Test SwissProt loading
    print("\n--- Step 3: Load SwissProt (via property) ---")
    _ = drug_adapter.swissprots
    mem_after_swissprot = print_memory("After loading SwissProt")
    print(f"  SwissProt count: {len(drug_adapter.swissprots)}")
    
    # Clean up
    drug_adapter._swissprots = None
    gc.collect()
    print_memory("After clearing SwissProt")
    
    # Step 4: Test DrugBank download
    print("\n--- Step 4: Download DrugBank data ---")
    from pypath.share import curl, settings
    from contextlib import ExitStack
    
    with ExitStack() as stack:
        stack.enter_context(settings.context(retries=3))
        
        print("  Downloading DrugBank node data...")
        drug_adapter.download_drugbank_node_data()
        print_memory("After download_drugbank_node_data")
        
        print(f"  drugbank_drugs_detailed count: {len(drug_adapter.drugbank_drugs_detailed) if hasattr(drug_adapter, 'drugbank_drugs_detailed') else 'N/A'}")
        print(f"  drugbank_drugs_external_ids count: {len(drug_adapter.drugbank_drugs_external_ids) if hasattr(drug_adapter, 'drugbank_drugs_external_ids') else 'N/A'}")
    
    # Step 5: Test UniChem mappings one by one
    print("\n--- Step 5: Test UniChem mappings individually ---")
    from pypath.inputs import unichem
    
    unichem_mappings = [
        ("drugbank", "zinc"),
        ("chembl", "drugbank"),
        ("drugbank", "bindingdb"),
        ("drugbank", "clinicaltrials"),
        ("drugbank", "chebi"),
        ("drugbank", "pubchem"),
    ]
    
    for src, tgt in unichem_mappings:
        gc.collect()
        mem_before = get_memory_mb()
        print(f"\n  Loading unichem.unichem_mapping('{src}', '{tgt}')...")
        try:
            mapping = unichem.unichem_mapping(src, tgt)
            mem_after = get_memory_mb()
            print(f"    Entries: {len(mapping)}, Memory: {mem_before:.1f} -> {mem_after:.1f} MB (delta: {mem_after - mem_before:.1f} MB)")
            del mapping
            gc.collect()
        except Exception as e:
            print(f"    ERROR: {e}")
    
    print_memory("After all UniChem tests")
    
    # Cleanup
    del drug_adapter
    gc.collect()
    print_memory("Final (after cleanup)")


def test_compound_adapter_step_by_step():
    """Test Compound adapter step by step to find memory bottleneck."""
    print("\n" + "=" * 80)
    print("Testing Compound Adapter Memory Usage Step by Step")
    print("=" * 80)
    
    initial_mem = print_memory("Initial")
    
    # Step 1: Import
    print("\n--- Step 1: Import modules ---")
    from bccb.compound_adapter import Compound
    print_memory("After import")
    
    # Step 2: Initialize
    print("\n--- Step 2: Initialize Compound adapter ---")
    compound_adapter = Compound(
        stitch_organism=9606,
        test_mode=True,
        low_memory_mode=True
    )
    print_memory("After Compound.__init__")
    
    # Step 3: Test ChEMBL downloads
    print("\n--- Step 3: Test ChEMBL data downloads ---")
    from pypath.inputs import chembl
    
    gc.collect()
    mem_before = get_memory_mb()
    print("  Loading chembl.chembl_molecules()...")
    try:
        molecules = chembl.chembl_molecules()
        molecules_list = list(molecules)[:1000]  # Only first 1000 for testing
        mem_after = get_memory_mb()
        print(f"    Sample count: {len(molecules_list)}, Memory delta: {mem_after - mem_before:.1f} MB")
        del molecules, molecules_list
        gc.collect()
    except Exception as e:
        print(f"    ERROR: {e}")
    
    gc.collect()
    mem_before = get_memory_mb()
    print("  Loading chembl.chembl_activities() (first 1000)...")
    try:
        activities = chembl.chembl_activities(standard_relation="=")
        count = 0
        for act in activities:
            count += 1
            if count >= 1000:
                break
        mem_after = get_memory_mb()
        print(f"    Sample count: {count}, Memory delta: {mem_after - mem_before:.1f} MB")
        del activities
        gc.collect()
    except Exception as e:
        print(f"    ERROR: {e}")
    
    # Cleanup
    del compound_adapter
    gc.collect()
    print_memory("Final (after cleanup)")


def test_go_adapter_memory():
    """Test GO adapter memory usage - this was where previous run failed."""
    print("\n" + "=" * 80)
    print("Testing GO Adapter Memory Usage")
    print("=" * 80)
    
    initial_mem = print_memory("Initial")
    
    from bccb.go_adapter import GO
    print_memory("After import")
    
    go_adapter = GO(organism=9606, test_mode=True)
    print_memory("After GO.__init__")
    
    from pypath.share import curl, settings
    from contextlib import ExitStack
    
    with ExitStack() as stack:
        stack.enter_context(settings.context(retries=3))
        go_adapter.download_go_data(cache=True)
    
    print_memory("After download_go_data")
    
    del go_adapter
    gc.collect()
    print_memory("Final (after cleanup)")


if __name__ == "__main__":
    print(f"System Memory: {psutil.virtual_memory().total / 1024 / 1024 / 1024:.1f} GB")
    print(f"Available Memory: {psutil.virtual_memory().available / 1024 / 1024 / 1024:.1f} GB")
    print()
    
    # Run tests
    try:
        test_drug_adapter_step_by_step()
    except Exception as e:
        print(f"Drug adapter test failed: {e}")
        import traceback
        traceback.print_exc()
    
    gc.collect()
    
    try:
        test_compound_adapter_step_by_step()
    except Exception as e:
        print(f"Compound adapter test failed: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "=" * 80)
    print("Memory debugging complete!")
    print("=" * 80)
