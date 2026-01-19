#!/usr/bin/env python3
"""
Test script for Side Effect adapter to verify it works correctly.
Run with: conda activate crossbarv2 && python scripts/test_side_effect_adapter.py
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

from bccb.side_effect_adapter import SideEffect

# Configuration
CACHE = True
TEST_MODE = True
output_dir_path = str(project_root / "biocypher-out")

# Credentials
drugbank_user = "zengchuanlong23@mails.ucas.ac.cn"
drugbank_passwd = "iHDTbZ3vpFaWzC"


def test_side_effect_download():
    """Test Side Effect data download."""
    print("=" * 60)
    print("Testing Side Effect Data Download")
    print("=" * 60)
    
    try:
        side_effect_adapter = SideEffect(
            drugbank_user=drugbank_user,
            drugbank_passwd=drugbank_passwd,
            export_csv=False,
            output_dir=output_dir_path,
            test_mode=TEST_MODE,
        )
        
        print("\nDownloading side effect data...")
        side_effect_adapter.download_side_effect_data(cache=CACHE)
        print("Download completed")
        
        # Check what was downloaded
        if hasattr(side_effect_adapter, 'meddra_id_to_side_effect_name'):
            print(f"MedDRA side effects: {len(side_effect_adapter.meddra_id_to_side_effect_name)}")
        
        if hasattr(side_effect_adapter, 'sider_meddra_with_freq'):
            print(f"SIDER with frequency: {len(side_effect_adapter.sider_meddra_with_freq)}")
        
        if hasattr(side_effect_adapter, 'offsides'):
            print(f"OffSides data: {len(side_effect_adapter.offsides)}")
        
        if hasattr(side_effect_adapter, 'adrecs_data'):
            print(f"ADReCS data: {len(side_effect_adapter.adrecs_data)}")
        
        print("\nSide Effect download test: PASSED")
        return True
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'side_effect_adapter' in dir():
            del side_effect_adapter
        gc.collect()
        print("Memory cleaned up")


def test_side_effect_nodes():
    """Test Side Effect node generation."""
    print("=" * 60)
    print("Testing Side Effect Node Generation")
    print("=" * 60)
    
    try:
        side_effect_adapter = SideEffect(
            drugbank_user=drugbank_user,
            drugbank_passwd=drugbank_passwd,
            export_csv=False,
            output_dir=output_dir_path,
            test_mode=TEST_MODE,
        )
        
        print("\nDownloading side effect data...")
        side_effect_adapter.download_side_effect_data(cache=CACHE)
        
        print("\nGenerating side effect nodes...")
        node_count = 0
        for node in side_effect_adapter.get_nodes():
            node_count += 1
            if node_count <= 5:
                node_id = node[0] if node[0] else "None"
                print(f"  Node {node_count}: {node_id[:50]}...")
        
        print(f"Generated {node_count} side effect nodes")
        
        print("\nSide Effect nodes test: PASSED")
        return True
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'side_effect_adapter' in dir():
            del side_effect_adapter
        gc.collect()
        print("Memory cleaned up")


def test_side_effect_edges():
    """Test Side Effect edge generation."""
    print("=" * 60)
    print("Testing Side Effect Edge Generation")
    print("=" * 60)
    
    try:
        side_effect_adapter = SideEffect(
            drugbank_user=drugbank_user,
            drugbank_passwd=drugbank_passwd,
            export_csv=False,
            output_dir=output_dir_path,
            test_mode=TEST_MODE,
        )
        
        print("\nDownloading side effect data...")
        side_effect_adapter.download_side_effect_data(cache=CACHE)
        
        print("\nGenerating side effect edges...")
        edge_count = 0
        for edge in side_effect_adapter.get_edges():
            edge_count += 1
            if edge_count <= 5:
                src = edge[0] if edge[0] else "None"
                tgt = edge[1] if edge[1] else "None"
                print(f"  Edge {edge_count}: {src[:30]}... -> {tgt[:30]}...")
        
        print(f"Generated {edge_count} side effect edges")
        
        print("\nSide Effect edges test: PASSED")
        return True
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'side_effect_adapter' in dir():
            del side_effect_adapter
        gc.collect()
        print("Memory cleaned up")


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("Side Effect Adapter Test Suite")
    print("=" * 60 + "\n")
    
    results = {}
    
    # Test 1: Download
    results['download'] = test_side_effect_download()
    
    # Test 2: Nodes
    results['nodes'] = test_side_effect_nodes()
    
    # Test 3: Edges
    results['edges'] = test_side_effect_edges()
    
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
        print("All tests PASSED! Side Effect adapter is ready for create_crossbar.py")
        sys.exit(0)
    else:
        print("Some tests FAILED! Please fix issues before running create_crossbar.py")
        sys.exit(1)
