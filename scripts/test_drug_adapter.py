#!/usr/bin/env python3
"""
Test script for Drug adapter to verify Pharos streaming and all DTI sources work correctly.
Run with: conda activate crossbarv2 && python scripts/test_drug_adapter.py
"""

import sys
import gc
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

# Setup pypath cache before any pypath imports
from bccb import cache_config
cache_config.setup_pypath_cache()

# Patch pypath settings for compatibility
try:
    from pypath.share import settings
    if not hasattr(settings, 'context') and hasattr(settings, 'settings'):
        settings.context = settings.settings.context
except ImportError:
    pass

from bccb.drug_adapter import Drug, DrugNodeField

# Configuration
CACHE = True
TEST_MODE = True
output_dir_path = str(project_root / "biocypher-out")

# Credentials
drugbank_user = "zengchuanlong23@mails.ucas.ac.cn"
drugbank_passwd = "iHDTbZ3vpFaWzC"


def test_pharos_streaming():
    """Test Pharos DTI streaming specifically."""
    print("=" * 60)
    print("Testing Pharos DTI Streaming")
    print("=" * 60)
    
    drug_node_fields = [f for f in DrugNodeField if f != DrugNodeField.SELFORMER_EMBEDDING]
    
    drug_adapter = Drug(
        drugbank_user=drugbank_user,
        drugbank_passwd=drugbank_passwd,
        node_fields=drug_node_fields,
        export_csv=False,
        output_dir=output_dir_path,
        test_mode=TEST_MODE,
        low_memory_mode=True,
    )
    
    # Need to initialize drugcentral_to_drugbank mapping first
    print("\nInitializing DrugBank data for ID mappings...")
    drug_adapter.download_drugbank_node_data()
    drug_adapter.process_drugbank_node_data()
    
    # Test streaming method directly
    print("\nTesting _stream_pharos_dti_data generator...")
    record_count = 0
    try:
        for record in drug_adapter._stream_pharos_dti_data(chunk_size=100, max_chunks=3):
            record_count += 1
            if record_count <= 5:
                print(f"  Sample record {record_count}: uniprot={record[0]}, activity={record[1]}")
    except Exception as e:
        print(f"ERROR in streaming: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    print(f"\nStreaming test: collected {record_count} records from 3 chunks")
    
    # Cleanup
    del drug_adapter
    gc.collect()
    
    print("Pharos streaming test: PASSED\n")
    return True


def test_full_drug_adapter():
    """Test full Drug adapter workflow as used in create_crossbar.py."""
    print("=" * 60)
    print("Testing Full Drug Adapter Workflow")
    print("=" * 60)
    
    drug_node_fields = [f for f in DrugNodeField if f != DrugNodeField.SELFORMER_EMBEDDING]
    
    try:
        drug_adapter = Drug(
            drugbank_user=drugbank_user,
            drugbank_passwd=drugbank_passwd,
            node_fields=drug_node_fields,
            export_csv=False,
            output_dir=output_dir_path,
            test_mode=TEST_MODE,
            low_memory_mode=True,
        )
        print("Drug adapter initialized successfully")
        
        # Download data
        print("\nDownloading drug data...")
        drug_adapter.download_drug_data(cache=CACHE)
        print("Download completed")
        
        # Process data
        print("\nProcessing drug data...")
        drug_adapter.process_drug_data()
        print("Processing completed")
        
        # Get nodes
        print("\nGenerating drug nodes...")
        node_count = 0
        for node in drug_adapter.get_drug_nodes():
            node_count += 1
            if node_count <= 3:
                print(f"  Sample node: {node[0][:50]}...")
        print(f"Generated {node_count} drug nodes")
        
        # Get edges
        print("\nGenerating edges...")
        edge_count = 0
        for edge in drug_adapter.get_edges():
            edge_count += 1
            if edge_count <= 3:
                # Handle None values safely
                src = edge[0] if edge[0] else "None"
                tgt = edge[1] if edge[1] else "None"
                print(f"  Sample edge: {src[:30]}... -> {tgt[:30]}...")
        print(f"Generated {edge_count} edges")
        
        print("\nFull drug adapter test: PASSED")
        return True
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'drug_adapter' in dir():
            del drug_adapter
        gc.collect()
        print("Memory cleaned up")


def test_pharos_process():
    """Test process_pharos_dti_data specifically."""
    print("=" * 60)
    print("Testing Pharos DTI Processing")
    print("=" * 60)
    
    drug_node_fields = [f for f in DrugNodeField if f != DrugNodeField.SELFORMER_EMBEDDING]
    
    try:
        drug_adapter = Drug(
            drugbank_user=drugbank_user,
            drugbank_passwd=drugbank_passwd,
            node_fields=drug_node_fields,
            export_csv=False,
            output_dir=output_dir_path,
            test_mode=TEST_MODE,
            low_memory_mode=True,
        )
        
        # Need drugcentral_to_drugbank mapping for processing
        print("Downloading DrugBank node data (for ID mappings)...")
        drug_adapter.download_drugbank_node_data()
        drug_adapter.process_drugbank_node_data()
        
        print("\nProcessing Pharos DTI data...")
        pharos_df = drug_adapter.process_pharos_dti_data()
        
        print(f"\nPharos DTI DataFrame shape: {pharos_df.shape}")
        print(f"Columns: {list(pharos_df.columns)}")
        if len(pharos_df) > 0:
            print(f"\nSample rows:")
            print(pharos_df.head(3).to_string())
        
        print("\nPharos processing test: PASSED")
        return True
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'drug_adapter' in dir():
            del drug_adapter
        gc.collect()
        print("Memory cleaned up")


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("Drug Adapter Test Suite")
    print("=" * 60 + "\n")
    
    results = {}
    
    # Test 1: Pharos streaming
    results['pharos_streaming'] = test_pharos_streaming()
    
    # Test 2: Pharos processing
    results['pharos_process'] = test_pharos_process()
    
    # Test 3: Full workflow
    results['full_workflow'] = test_full_drug_adapter()
    
    # Summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    
    all_passed = True
    for test_name, passed in results.items():
        status = "PASSED" if passed else "FAILED"
        print(f"  {test_name}: {status}")
        if not passed:
            all_passed = False
    
    print("=" * 60)
    if all_passed:
        print("All tests PASSED! Drug adapter is ready for create_crossbar.py")
        sys.exit(0)
    else:
        print("Some tests FAILED! Please fix issues before running create_crossbar.py")
        sys.exit(1)
