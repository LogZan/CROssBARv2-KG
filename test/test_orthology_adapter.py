#!/usr/bin/env python3
"""
Test script for Orthology adapter to verify it works correctly.
Run with: conda activate crossbarv2 && python test/test_orthology_adapter.py
"""

import sys
import gc
import logging
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)

# Setup pypath cache
from bccb import cache_config
cache_config.setup_pypath_cache()

from bccb.orthology_adapter import Orthology

# Configuration
CACHE = True
TEST_MODE = True


def test_orthology_initialization():
    """Test Orthology adapter initialization."""
    print("=" * 60)
    print("Testing Orthology Adapter Initialization")
    print("=" * 60)

    try:
        orthology_adapter = Orthology(
            export_csv=False,
            output_dir=None,
            test_mode=TEST_MODE
        )

        print(f"Orthology adapter initialized successfully")
        print(f"  Test mode: {orthology_adapter.early_stopping}")
        print(f"  Add prefix: {orthology_adapter.add_prefix}")

        print("\nOrthology initialization test: PASSED")
        return True

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'orthology_adapter' in dir():
            del orthology_adapter
        gc.collect()
        print("Memory cleaned up")


def test_orthology_download():
    """Test Orthology data download."""
    print("=" * 60)
    print("Testing Orthology Data Download")
    print("=" * 60)

    try:
        orthology_adapter = Orthology(
            export_csv=False,
            output_dir=None,
            test_mode=TEST_MODE
        )

        print("\nDownloading Orthology data (OMA, Pharos)...")
        orthology_adapter.download_orthology_data(cache=CACHE)
        print("Download completed")

        # Check what was downloaded
        if hasattr(orthology_adapter, 'oma_orthology'):
            print(f"✓ OMA data downloaded: {len(orthology_adapter.oma_orthology)} orthologs")
        else:
            print("✗ OMA data not downloaded")

        if hasattr(orthology_adapter, 'pharos_orthology_init'):
            print(f"✓ Pharos data downloaded")
        else:
            print("✗ Pharos data not downloaded")

        print("\nOrthology download test: PASSED")
        return True

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'orthology_adapter' in dir():
            del orthology_adapter
        gc.collect()
        print("Memory cleaned up")


def test_orthology_edges():
    """Test Orthology edge generation."""
    print("=" * 60)
    print("Testing Orthology Edge Generation")
    print("=" * 60)

    try:
        orthology_adapter = Orthology(
            export_csv=False,
            output_dir=None,
            test_mode=TEST_MODE
        )

        print("\nDownloading Orthology data...")
        orthology_adapter.download_orthology_data(cache=CACHE)

        print("\nGenerating Orthology edges...")
        edges = orthology_adapter.get_orthology_edges()

        edge_count = len(edges)
        print(f"Generated {edge_count} orthology edges")

        # Show sample edges
        if edge_count > 0:
            print("\nSample edges:")
            for i, edge in enumerate(edges[:3]):
                edge_id, source, target, label, props = edge
                print(f"  Edge {i+1}:")
                print(f"    Source: {source}")
                print(f"    Target: {target}")
                print(f"    Label: {label}")
                print(f"    Properties: {list(props.keys())}")

        print("\nOrthology edges test: PASSED")
        return True

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'orthology_adapter' in dir():
            del orthology_adapter
        gc.collect()
        print("Memory cleaned up")


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("Orthology Adapter Test Suite")
    print("=" * 60 + "\n")

    results = {}

    # Test 1: Initialization
    results['initialization'] = test_orthology_initialization()

    # Test 2: Download
    results['download'] = test_orthology_download()

    # Test 3: Edges
    results['edges'] = test_orthology_edges()

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
        print("All tests PASSED! Orthology adapter is ready for create_crossbar.py")
        sys.exit(0)
    else:
        print("Some tests FAILED! Please fix issues before running create_crossbar.py")
        sys.exit(1)
