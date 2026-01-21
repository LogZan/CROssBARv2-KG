#!/usr/bin/env python3
"""
Test script for InterPro adapter.
Run with: conda activate crossbarv2 && python test/test_interpro_adapter.py
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

from bccb.interpro_adapter import InterPro

# Configuration
CACHE = True
TEST_MODE = True
USE_EMBEDDINGS = False  # Set to True to test with embeddings


def test_reviewed_api_10_proteins():
    """Test reviewed InterPro API in test mode (download + parse)."""
    print("=" * 60)
    print("Testing InterPro API (reviewed, test mode)")
    print("=" * 60)

    try:
        interpro_adapter = InterPro(
            organism="*",
            test_mode=TEST_MODE
        )

        print("\nDownloading reviewed proteins (test mode)...")
        interpro_adapter.download_domain_node_data(
            cache=CACHE
        )

        protein_count = len(getattr(interpro_adapter, "protein_domains", {}))
        entry_count = len(getattr(interpro_adapter, "interpro_entries", {}))
        print(f"Proteins with domains: {protein_count}")
        print(f"InterPro entries collected: {entry_count}")

        sample = None
        for protein_id, domains in getattr(interpro_adapter, "protein_domains", {}).items():
            if domains:
                sample = (protein_id, domains[0])
                break

        if not sample:
            raise RuntimeError("No domain mappings found in sampled proteins.")

        protein_id, domain = sample
        print(f"Sample protein: {protein_id}")
        print(f"Sample domain: {domain}")

        if domain.get("start") is None or domain.get("end") is None:
            raise RuntimeError("Parsed domain missing start/end.")

        print("\nInterPro API test: PASSED")
        return True

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'interpro_adapter' in dir():
            del interpro_adapter
        gc.collect()
        print("Memory cleaned up")


def test_interpro_download():
    """Test InterPro data download."""
    print("=" * 60)
    print("Testing InterPro Data Download")
    print("=" * 60)

    try:
        interpro_adapter = InterPro(
            organism="*",
            test_mode=TEST_MODE
        )

        print("\nDownloading InterPro data...")
        interpro_adapter.download_interpro_data(cache=CACHE)
        print("InterPro data download completed")

        # Check what was downloaded
        if hasattr(interpro_adapter, 'interpro_data'):
            print(f"InterPro data loaded: {len(interpro_adapter.interpro_data)} entries")

        print("\nInterPro download test: PASSED")
        return True

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'interpro_adapter' in dir():
            del interpro_adapter
        gc.collect()
        print("Memory cleaned up")


def test_domain_node_download():
    """Test domain node data download."""
    print("=" * 60)
    print("Testing Domain Node Data Download")
    print("=" * 60)

    try:
        interpro_adapter = InterPro(
            organism="*",
            test_mode=TEST_MODE
        )

        print("\nDownloading InterPro data...")
        interpro_adapter.download_interpro_data(cache=CACHE)

        print("\nDownloading domain node data...")
        embedding_path = None
        if USE_EMBEDDINGS:
            embedding_path = '/GenSIvePFS/users/data/embeddings/dom2vec_domain_embedding.h5'
            print(f"Using embedding: {embedding_path}")

        interpro_adapter.download_domain_node_data(
            dom2vec_embedding_path=embedding_path,
            cache=CACHE
        )
        print("Domain node data download completed")

        # Check what was downloaded
        if hasattr(interpro_adapter, 'domain_data'):
            print(f"Domain data loaded: {len(interpro_adapter.domain_data)} domains")

        print("\nDomain node download test: PASSED")
        return True

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'interpro_adapter' in dir():
            del interpro_adapter
        gc.collect()
        print("Memory cleaned up")


def test_interpro_nodes():
    """Test InterPro node generation."""
    print("=" * 60)
    print("Testing InterPro Node Generation")
    print("=" * 60)

    try:
        interpro_adapter = InterPro(
            organism="*",
            test_mode=TEST_MODE
        )

        print("\nDownloading data...")
        interpro_adapter.download_interpro_data(cache=CACHE)
        interpro_adapter.download_domain_node_data(cache=CACHE)

        print("\nGenerating InterPro nodes...")
        node_count = 0
        for node in interpro_adapter.get_interpro_nodes():
            node_count += 1
            if node_count <= 3:
                node_id = node[0] if node[0] else 'None'
                node_type = node[1] if len(node) > 1 else 'Unknown'
                print(f"  Node {node_count}: {node_id[:50]}... (type: {node_type})")
            if node_count >= 100:
                break

        print(f"Generated {node_count} InterPro nodes (limited to 100)")

        print("\nInterPro nodes test: PASSED")
        return True

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'interpro_adapter' in dir():
            del interpro_adapter
        gc.collect()
        print("Memory cleaned up")


def test_interpro_edges():
    """Test InterPro edge generation."""
    print("=" * 60)
    print("Testing InterPro Edge Generation")
    print("=" * 60)

    try:
        interpro_adapter = InterPro(
            organism="*",
            test_mode=TEST_MODE
        )

        print("\nDownloading data...")
        interpro_adapter.download_interpro_data(cache=CACHE)
        interpro_adapter.download_domain_node_data(cache=CACHE)

        print("\nGenerating InterPro edges...")
        edge_count = 0
        for edge in interpro_adapter.get_interpro_edges():
            edge_count += 1
            if edge_count <= 3:
                src = edge[1] if len(edge) > 1 and edge[1] else "None"
                tgt = edge[2] if len(edge) > 2 and edge[2] else "None"
                label = edge[3] if len(edge) > 3 and edge[3] else "Unknown"
                print(f"  Edge {edge_count}: {src[:30]}... --[{label}]--> {tgt[:30]}...")
            if edge_count >= 100:
                break

        print(f"Generated {edge_count} InterPro edges (limited to 100)")

        print("\nInterPro edges test: PASSED")
        return True

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'interpro_adapter' in dir():
            del interpro_adapter
        gc.collect()
        print("Memory cleaned up")


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("InterPro Adapter Test Suite")
    print("=" * 60 + "\n")

    results = {}

    # Test 0: InterPro API (reviewed, 10 proteins)
    results['api_reviewed_10'] = test_reviewed_api_10_proteins()

    # Test 1: InterPro Download
    results['interpro_download'] = test_interpro_download()

    # Test 2: Domain Node Download
    results['domain_download'] = test_domain_node_download()

    # Test 3: Nodes
    results['nodes'] = test_interpro_nodes()

    # Test 4: Edges
    results['edges'] = test_interpro_edges()

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
        print("All tests PASSED! InterPro adapter is ready for create_crossbar.py")
        sys.exit(0)
    else:
        print("Some tests FAILED! Please fix issues before running create_crossbar.py")
        sys.exit(1)
