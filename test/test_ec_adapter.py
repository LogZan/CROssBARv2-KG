#!/usr/bin/env python3
"""
Test script for EC (Enzyme Commission) adapter to verify it works correctly.
Run with: conda activate crossbarv2 && python scripts/test_ec_adapter.py
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

from bccb.ec_adapter import EC

# Configuration
CACHE = True
TEST_MODE = True
output_dir_path = str(project_root / "biocypher-out")


def test_ec_download():
    """Test EC data download."""
    print("=" * 60)
    print("Testing EC Data Download")
    print("=" * 60)
    
    try:
        ec_adapter = EC(
            export_csv=False,
            output_dir=output_dir_path,
            test_mode=TEST_MODE
        )
        
        print("\nDownloading EC data...")
        ec_adapter.download_ec_data(cache=CACHE)
        print("Download completed")
        
        # Check what was downloaded
        if hasattr(ec_adapter, 'enzymes'):
            print(f"Enzymes: {len(ec_adapter.enzymes) if ec_adapter.enzymes else 0}")
        
        if hasattr(ec_adapter, 'enzyme_classes'):
            print(f"Enzyme classes: {len(ec_adapter.enzyme_classes) if ec_adapter.enzyme_classes else 0}")
        
        print("\nEC download test: PASSED")
        return True
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'ec_adapter' in dir():
            del ec_adapter
        gc.collect()
        print("Memory cleaned up")


def test_ec_nodes():
    """Test EC node generation."""
    print("=" * 60)
    print("Testing EC Node Generation")
    print("=" * 60)
    
    try:
        ec_adapter = EC(
            export_csv=False,
            output_dir=output_dir_path,
            test_mode=TEST_MODE
        )
        
        print("\nDownloading EC data...")
        ec_adapter.download_ec_data(cache=CACHE)
        
        print("\nGenerating EC nodes...")
        node_count = 0
        for node in ec_adapter.get_nodes():
            node_count += 1
            if node_count <= 5:
                node_id = node[0] if node[0] else "None"
                print(f"  Node {node_count}: {node_id[:50]}...")
        
        print(f"Generated {node_count} EC nodes")
        
        print("\nEC nodes test: PASSED")
        return True
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'ec_adapter' in dir():
            del ec_adapter
        gc.collect()
        print("Memory cleaned up")


def test_ec_edges():
    """Test EC edge generation."""
    print("=" * 60)
    print("Testing EC Edge Generation")
    print("=" * 60)
    
    try:
        ec_adapter = EC(
            export_csv=False,
            output_dir=output_dir_path,
            test_mode=TEST_MODE
        )
        
        print("\nDownloading EC data...")
        ec_adapter.download_ec_data(cache=CACHE)
        
        print("\nGenerating EC edges...")
        edge_count = 0
        for edge in ec_adapter.get_edges():
            edge_count += 1
            if edge_count <= 5:
                src = str(edge[0])[:30] if edge[0] else "None"
                tgt = str(edge[1])[:30] if edge[1] else "None"
                print(f"  Edge {edge_count}: {src}... -> {tgt}...")
        
        print(f"Generated {edge_count} EC edges")
        
        print("\nEC edges test: PASSED")
        return True
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'ec_adapter' in dir():
            del ec_adapter
        gc.collect()
        print("Memory cleaned up")


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("EC Adapter Test Suite")
    print("=" * 60 + "\n")
    
    results = {}
    
    # Test 1: Download
    results['download'] = test_ec_download()
    
    # Test 2: Nodes
    results['nodes'] = test_ec_nodes()
    
    # Test 3: Edges
    results['edges'] = test_ec_edges()
    
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
        print("All tests PASSED! EC adapter is ready for create_crossbar.py")
        sys.exit(0)
    else:
        print("Some tests FAILED! Please fix issues before running create_crossbar.py")
        sys.exit(1)
