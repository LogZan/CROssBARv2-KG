#!/usr/bin/env python3
"""
Test script for GO adapter to verify streaming and memory optimization works correctly.
Run with: conda activate crossbarv2 && python scripts/test_go_adapter.py
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

from bccb.go_adapter import GO

# Configuration
CACHE = True
TEST_MODE = True


def test_go_streaming():
    """Test GO annotation streaming specifically."""
    print("=" * 60)
    print("Testing GO Annotation Streaming")
    print("=" * 60)
    
    try:
        go_adapter = GO(
            organism=9606,
            test_mode=TEST_MODE
        )
        
        # Test streaming method directly with small limit
        print("\nTesting _stream_go_annotations...")
        allowed_proteins = {"P04637", "P00533", "P38398", "O00206", "P01116"}
        
        result = go_adapter._stream_go_annotations(
            organism=9606,
            fields=["qualifier", "go_id", "reference", "evidence_code"],
            allowed_proteins=allowed_proteins,
            max_proteins=5,
        )
        
        print(f"Streamed annotations for {len(result)} proteins")
        for protein_id, annotations in list(result.items())[:3]:
            print(f"  {protein_id}: {len(annotations)} annotations")
            if annotations:
                sample = list(annotations)[0]
                print(f"    Sample: {sample}")
        
        print("\nGO streaming test: PASSED")
        return True
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'go_adapter' in dir():
            del go_adapter
        gc.collect()
        print("Memory cleaned up")


def test_go_download():
    """Test GO data download with streaming."""
    print("=" * 60)
    print("Testing GO Data Download")
    print("=" * 60)
    
    try:
        go_adapter = GO(
            organism=9606,
            test_mode=TEST_MODE
        )
        
        print("\nDownloading GO data...")
        go_adapter.download_go_data(cache=CACHE)
        print("Download completed")
        
        # Check what was downloaded
        if hasattr(go_adapter, 'go_ontology'):
            print(f"GO ontology loaded: {type(go_adapter.go_ontology)}")
        
        if hasattr(go_adapter, 'go_annots'):
            print(f"GO annotations: {len(go_adapter.go_annots)} proteins")
            for protein_id in list(go_adapter.go_annots.keys())[:3]:
                print(f"  {protein_id}: {len(go_adapter.go_annots[protein_id])} annotations")
        
        if hasattr(go_adapter, 'swissprots'):
            print(f"SwissProt proteins: {len(go_adapter.swissprots)}")
        
        print("\nGO download test: PASSED")
        return True
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'go_adapter' in dir():
            del go_adapter
        gc.collect()
        print("Memory cleaned up")


def test_go_nodes():
    """Test GO node generation."""
    print("=" * 60)
    print("Testing GO Node Generation")
    print("=" * 60)
    
    try:
        go_adapter = GO(
            organism=9606,
            test_mode=TEST_MODE
        )
        
        print("\nDownloading GO data...")
        go_adapter.download_go_data(cache=CACHE)
        
        print("\nGenerating GO nodes...")
        node_count = 0
        for node in go_adapter.get_go_nodes():
            node_count += 1
            if node_count <= 3:
                print(f"  Node {node_count}: {node[0][:50] if node[0] else 'None'}...")
            if node_count >= 100:
                break
        
        print(f"Generated {node_count} GO nodes (limited to 100)")
        
        print("\nGO nodes test: PASSED")
        return True
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'go_adapter' in dir():
            del go_adapter
        gc.collect()
        print("Memory cleaned up")


def test_go_edges():
    """Test GO edge generation."""
    print("=" * 60)
    print("Testing GO Edge Generation")
    print("=" * 60)
    
    try:
        go_adapter = GO(
            organism=9606,
            test_mode=TEST_MODE
        )
        
        print("\nDownloading GO data...")
        go_adapter.download_go_data(cache=CACHE)
        
        print("\nGenerating GO edges...")
        edge_count = 0
        for edge in go_adapter.get_go_edges():
            edge_count += 1
            if edge_count <= 3:
                src = edge[0] if edge[0] else "None"
                tgt = edge[1] if edge[1] else "None"
                print(f"  Edge {edge_count}: {src[:30]}... -> {tgt[:30]}...")
            if edge_count >= 100:
                break
        
        print(f"Generated {edge_count} GO edges (limited to 100)")
        
        print("\nGO edges test: PASSED")
        return True
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'go_adapter' in dir():
            del go_adapter
        gc.collect()
        print("Memory cleaned up")


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("GO Adapter Test Suite")
    print("=" * 60 + "\n")
    
    results = {}
    
    # Test 1: Streaming
    results['streaming'] = test_go_streaming()
    
    # Test 2: Download
    results['download'] = test_go_download()
    
    # Test 3: Nodes
    results['nodes'] = test_go_nodes()
    
    # Test 4: Edges
    results['edges'] = test_go_edges()
    
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
        print("All tests PASSED! GO adapter is ready for create_crossbar.py")
        sys.exit(0)
    else:
        print("Some tests FAILED! Please fix issues before running create_crossbar.py")
        sys.exit(1)
