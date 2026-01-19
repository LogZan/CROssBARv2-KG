#!/usr/bin/env python3
"""
Test script for Pathway adapter to verify it works correctly.
Run with: conda activate crossbarv2 && python scripts/test_pathway_adapter.py
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

from bccb.pathway_adapter import Pathway

# Configuration
CACHE = True
TEST_MODE = True
output_dir_path = str(project_root / "biocypher-out")

# Credentials
drugbank_user = "zengchuanlong23@mails.ucas.ac.cn"
drugbank_passwd = "iHDTbZ3vpFaWzC"


def test_pathway_download():
    """Test Pathway data download."""
    print("=" * 60)
    print("Testing Pathway Data Download")
    print("=" * 60)
    
    try:
        # Limit to human (hsa) only for faster testing
        pathway_adapter = Pathway(
            drugbank_user=drugbank_user,
            drugbank_passwd=drugbank_passwd,
            export_csv=False,
            output_dir=output_dir_path,
            test_mode=TEST_MODE,
            kegg_organism=["hsa"],  # Only human to speed up KEGG download
        )
        
        print("\nDownloading pathway data (limited to human organism)...")
        pathway_adapter.download_pathway_data(cache=CACHE)
        print("Download completed")
        
        # Check what was downloaded
        if hasattr(pathway_adapter, 'reactome_pathways'):
            print(f"Reactome pathways: {len(pathway_adapter.reactome_pathways)}")
        
        if hasattr(pathway_adapter, 'kegg_pathways'):
            print(f"KEGG pathways: {len(pathway_adapter.kegg_pathways)}")
        
        if hasattr(pathway_adapter, 'reactome_uniprot_pathway'):
            print(f"Reactome protein-pathway: {len(pathway_adapter.reactome_uniprot_pathway)}")
        
        if hasattr(pathway_adapter, 'kegg_gene_to_pathway'):
            print(f"KEGG gene-pathway: {len(pathway_adapter.kegg_gene_to_pathway)}")
        
        print("\nPathway download test: PASSED")
        return True
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'pathway_adapter' in dir():
            del pathway_adapter
        gc.collect()
        print("Memory cleaned up")


def test_pathway_nodes():
    """Test Pathway node generation."""
    print("=" * 60)
    print("Testing Pathway Node Generation")
    print("=" * 60)
    
    try:
        pathway_adapter = Pathway(
            drugbank_user=drugbank_user,
            drugbank_passwd=drugbank_passwd,
            export_csv=False,
            output_dir=output_dir_path,
            test_mode=TEST_MODE,
            kegg_organism=["hsa"],
        )
        
        print("\nDownloading pathway data...")
        pathway_adapter.download_pathway_data(cache=CACHE)
        
        print("\nGenerating pathway nodes...")
        node_count = 0
        for node in pathway_adapter.get_nodes():
            node_count += 1
            if node_count <= 5:
                node_id = node[0] if node[0] else "None"
                print(f"  Node {node_count}: {node_id[:50]}...")
        
        print(f"Generated {node_count} pathway nodes")
        
        print("\nPathway nodes test: PASSED")
        return True
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'pathway_adapter' in dir():
            del pathway_adapter
        gc.collect()
        print("Memory cleaned up")


def test_pathway_edges():
    """Test Pathway edge generation."""
    print("=" * 60)
    print("Testing Pathway Edge Generation")
    print("=" * 60)
    
    try:
        pathway_adapter = Pathway(
            drugbank_user=drugbank_user,
            drugbank_passwd=drugbank_passwd,
            export_csv=False,
            output_dir=output_dir_path,
            test_mode=TEST_MODE,
            kegg_organism=["hsa"],
        )
        
        print("\nDownloading pathway data...")
        pathway_adapter.download_pathway_data(cache=CACHE)
        
        print("\nGenerating pathway edges...")
        edge_count = 0
        for edge in pathway_adapter.get_edges():
            edge_count += 1
            if edge_count <= 5:
                src = edge[0] if edge[0] else "None"
                tgt = edge[1] if edge[1] else "None"
                print(f"  Edge {edge_count}: {src[:30]}... -> {tgt[:30]}...")
        
        print(f"Generated {edge_count} pathway edges")
        
        print("\nPathway edges test: PASSED")
        return True
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if 'pathway_adapter' in dir():
            del pathway_adapter
        gc.collect()
        print("Memory cleaned up")


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("Pathway Adapter Test Suite")
    print("=" * 60 + "\n")
    
    results = {}
    
    # Test 1: Download
    results['download'] = test_pathway_download()
    
    # Test 2: Nodes
    results['nodes'] = test_pathway_nodes()
    
    # Test 3: Edges
    results['edges'] = test_pathway_edges()
    
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
        print("All tests PASSED! Pathway adapter is ready for create_crossbar.py")
        sys.exit(0)
    else:
        print("Some tests FAILED! Please fix issues before running create_crossbar.py")
        sys.exit(1)
