#!/usr/bin/env python3
"""
Test script for TFGene (Transcription Factor - Gene) adapter to verify it works correctly.
Run with: conda activate crossbarv2 && python scripts/test_tfgene_adapter.py
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

from bccb.tfgen_adapter import TFGene

# Configuration
CACHE = True
TEST_MODE = True
output_dir_path = str(project_root / "biocypher-out")


def test_tfgene_download():
    """Test TFGene data download."""
    print("=" * 60)
    print("Testing TFGene Data Download")
    print("=" * 60)
    
    try:
        tfgene_adapter = TFGene(
            export_csv=False,
            output_dir=output_dir_path,
            test_mode=TEST_MODE
        )
        
        print("\nDownloading TFGene data...")
        tfgene_adapter.download_tfgen_data(cache=CACHE)
        print("Download completed")
        
        # Check what was downloaded
        if hasattr(tfgene_adapter, 'collectri_edges'):
            print(f"CollecTRI edges: {len(tfgene_adapter.collectri_edges) if tfgene_adapter.collectri_edges else 0}")
        
        if hasattr(tfgene_adapter, 'dorothea_edges'):
            print(f"DoRothEA edges: {len(tfgene_adapter.dorothea_edges) if tfgene_adapter.dorothea_edges else 0}")
        
        if hasattr(tfgene_adapter, 'trrust_edges'):
            print(f"TRRUST edges: {len(tfgene_adapter.trrust_edges) if tfgene_adapter.trrust_edges else 0}")
        
        print("\nTFGene download test: PASSED")
        return True
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'tfgene_adapter' in dir():
            del tfgene_adapter
        gc.collect()
        print("Memory cleaned up")


def test_tfgene_edges():
    """Test TFGene edge generation."""
    print("=" * 60)
    print("Testing TFGene Edge Generation")
    print("=" * 60)
    
    try:
        tfgene_adapter = TFGene(
            export_csv=False,
            output_dir=output_dir_path,
            test_mode=TEST_MODE
        )
        
        print("\nDownloading TFGene data...")
        tfgene_adapter.download_tfgen_data(cache=CACHE)
        
        print("\nGenerating TFGene edges...")
        edge_count = 0
        for edge in tfgene_adapter.get_edges():
            edge_count += 1
            if edge_count <= 5:
                src = str(edge[0])[:30] if edge[0] else "None"
                tgt = str(edge[1])[:30] if edge[1] else "None"
                print(f"  Edge {edge_count}: {src}... -> {tgt}...")
        
        print(f"Generated {edge_count} TFGene edges")
        
        print("\nTFGene edges test: PASSED")
        return True
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'tfgene_adapter' in dir():
            del tfgene_adapter
        gc.collect()
        print("Memory cleaned up")


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("TFGene Adapter Test Suite")
    print("=" * 60 + "\n")
    
    results = {}
    
    # Test 1: Download
    results['download'] = test_tfgene_download()
    
    # Test 2: Edges (TFGene only has edges, no nodes)
    results['edges'] = test_tfgene_edges()
    
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
        print("All tests PASSED! TFGene adapter is ready for create_crossbar.py")
        sys.exit(0)
    else:
        print("Some tests FAILED! Please fix issues before running create_crossbar.py")
        sys.exit(1)
