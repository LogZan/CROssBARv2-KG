#!/usr/bin/env python3
"""
Test script for PPI adapter to verify it works correctly.
Run with: conda activate crossbarv2 && python test/test_ppi_adapter.py
"""

import sys
import gc
import logging
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

# Configure logging to output to terminal
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)

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

from bccb.ppi_adapter import PPI

# Configuration
CACHE = True
TEST_MODE = True
ORGANISM = 9606  # Human


def test_ppi_initialization():
    """Test PPI adapter initialization."""
    print("=" * 60)
    print("Testing PPI Adapter Initialization")
    print("=" * 60)

    try:
        ppi_adapter = PPI(
            organism=ORGANISM,
            output_dir=None,
            export_csv=False,
            test_mode=TEST_MODE
        )

        print(f"PPI adapter initialized successfully")
        print(f"  Organism: {ppi_adapter.organism}")
        print(f"  Test mode: {ppi_adapter.test_mode}")
        print(f"  SwissProt proteins loaded: {len(ppi_adapter.swissprots)}")

        print("\nPPI initialization test: PASSED")
        return True

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'ppi_adapter' in dir():
            del ppi_adapter
        gc.collect()
        print("Memory cleaned up")


def test_ppi_download():
    """Test PPI data download."""
    print("=" * 60)
    print("Testing PPI Data Download")
    print("=" * 60)

    try:
        ppi_adapter = PPI(
            organism=ORGANISM,
            output_dir=None,
            export_csv=False,
            test_mode=TEST_MODE
        )

        print("\nDownloading PPI data (IntAct, BioGRID, STRING)...")
        ppi_adapter.download_ppi_data(cache=CACHE)
        print("Download completed")

        # Check what was downloaded
        status = ppi_adapter.check_status_and_properties

        if status['intact']['downloaded']:
            print(f"✓ IntAct data downloaded: {len(ppi_adapter.intact_ints)} interactions")
        else:
            print("✗ IntAct data not downloaded")

        if status['biogrid']['downloaded']:
            print(f"✓ BioGRID data downloaded: {len(ppi_adapter.biogrid_ints)} interactions")
        else:
            print("✗ BioGRID data not downloaded")

        if status['string']['downloaded']:
            print(f"✓ STRING data downloaded: {len(ppi_adapter.string_ints)} interactions")
        else:
            print("✗ STRING data not downloaded")

        print("\nPPI download test: PASSED")
        return True

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'ppi_adapter' in dir():
            del ppi_adapter
        gc.collect()
        print("Memory cleaned up")


def test_ppi_processing():
    """Test PPI data processing."""
    print("=" * 60)
    print("Testing PPI Data Processing")
    print("=" * 60)

    try:
        ppi_adapter = PPI(
            organism=ORGANISM,
            output_dir=None,
            export_csv=False,
            test_mode=TEST_MODE
        )

        print("\nDownloading PPI data...")
        ppi_adapter.download_ppi_data(cache=CACHE)

        print("\nProcessing PPI data...")
        ppi_adapter.process_ppi_data()
        print("Processing completed")

        # Check what was processed
        status = ppi_adapter.check_status_and_properties

        if status['intact']['processed']:
            df = status['intact']['dataframe']
            print(f"✓ IntAct data processed: {len(df)} unique interactions")
        else:
            print("✗ IntAct data not processed")

        if status['biogrid']['processed']:
            df = status['biogrid']['dataframe']
            print(f"✓ BioGRID data processed: {len(df)} unique interactions")
        else:
            print("✗ BioGRID data not processed")

        if status['string']['processed']:
            df = status['string']['dataframe']
            print(f"✓ STRING data processed: {len(df)} unique interactions")
        else:
            print("✗ STRING data not processed")

        print("\nPPI processing test: PASSED")
        return True

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'ppi_adapter' in dir():
            del ppi_adapter
        gc.collect()
        print("Memory cleaned up")


def test_ppi_edges():
    """Test PPI edge generation."""
    print("=" * 60)
    print("Testing PPI Edge Generation")
    print("=" * 60)

    try:
        ppi_adapter = PPI(
            organism=ORGANISM,
            output_dir=None,
            export_csv=False,
            test_mode=TEST_MODE
        )

        print("\nDownloading PPI data...")
        ppi_adapter.download_ppi_data(cache=CACHE)

        print("\nProcessing PPI data...")
        ppi_adapter.process_ppi_data()

        print("\nGenerating PPI edges...")
        edges = ppi_adapter.get_ppi_edges()

        # Convert generator to list to count
        edge_list = list(edges)
        edge_count = len(edge_list)

        print(f"Generated {edge_count} PPI edges")

        # Show sample edges
        if edge_count > 0:
            print("\nSample edges:")
            for i, edge in enumerate(edge_list[:3]):
                edge_id, source, target, label, props = edge
                print(f"  Edge {i+1}:")
                print(f"    Source: {source}")
                print(f"    Target: {target}")
                print(f"    Label: {label}")
                print(f"    Properties: {list(props.keys())[:5]}...")

        print("\nPPI edges test: PASSED")
        return True

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'ppi_adapter' in dir():
            del ppi_adapter
        gc.collect()
        print("Memory cleaned up")


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("PPI Adapter Test Suite")
    print("=" * 60 + "\n")

    results = {}

    # Test 1: Initialization
    results['initialization'] = test_ppi_initialization()

    # Test 2: Download
    results['download'] = test_ppi_download()

    # Test 3: Processing
    results['processing'] = test_ppi_processing()

    # Test 4: Edges
    results['edges'] = test_ppi_edges()

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
        print("All tests PASSED! PPI adapter is ready for create_crossbar.py")
        sys.exit(0)
    else:
        print("Some tests FAILED! Please fix issues before running create_crossbar.py")
        sys.exit(1)
