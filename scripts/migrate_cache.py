#!/usr/bin/env python
"""
Migrate PyPath Cache Files to Categorized Directories

This script moves existing cache files to category-specific subdirectories.
"""
import os
import shutil
import logging
from pathlib import Path
from collections import defaultdict

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Cache directory
CACHE_BASE = "/GenSIvePFS/users/data/pypath_cache"

# Category directories
CATEGORIES = {
    'string': f'{CACHE_BASE}/string',
    'uniprot': f'{CACHE_BASE}/uniprot',
    'chembl': f'{CACHE_BASE}/chembl',
    'reactome': f'{CACHE_BASE}/reactome',
    'go': f'{CACHE_BASE}/go',
    'kegg': f'{CACHE_BASE}/kegg',
    'drugbank': f'{CACHE_BASE}/drugbank',
    'hpo': f'{CACHE_BASE}/hpo',
    'pharos': f'{CACHE_BASE}/pharos',
    'other': f'{CACHE_BASE}/other',
}


def categorize_file(filename: str) -> str:
    """Determine the category for a file based on its name."""
    name_lower = filename.lower()
    
    # STRING protein interaction files
    if '.protein.' in name_lower and ('.links.' in name_lower or '.physical.' in name_lower):
        return 'string'
    if 'species.v12' in name_lower:
        return 'string'
    
    # ChEMBL files
    if 'activity.json' in name_lower or 'mechanism.json' in name_lower:
        return 'chembl'
    if 'target.json' in name_lower or 'assay.json' in name_lower:
        return 'chembl'
    if 'document.json' in name_lower:
        return 'chembl'
    
    # Reactome files
    if 'reactome' in name_lower:
        return 'reactome'
    
    # GO files
    if name_lower.endswith('-terms') or '.obo' in name_lower:
        return 'go'
    if 'go.obo' in name_lower or 'gene_ontology' in name_lower:
        return 'go'
    
    # UniProt files
    if 'uniprot' in name_lower or 'swissprot' in name_lower:
        return 'uniprot'
    if 'sec_ac' in name_lower or 'speclist' in name_lower:
        return 'uniprot'
    
    # HPO files
    if 'hpo' in name_lower or 'phenotype' in name_lower:
        return 'hpo'
    
    # KEGG files
    if 'kegg' in name_lower:
        return 'kegg'
    
    # DrugBank files
    if 'drugbank' in name_lower or name_lower.startswith('d0'):
        return 'drugbank'
    
    # Pharos files
    if 'pharos' in name_lower or 'tcrd' in name_lower:
        return 'pharos'
    
    # Default to other
    return 'other'


def migrate_cache(dry_run: bool = False):
    """
    Migrate cache files to categorized directories.
    
    Args:
        dry_run: If True, only report what would be done without moving files
    """
    base_path = Path(CACHE_BASE)
    
    # Create category directories
    for category, path in CATEGORIES.items():
        os.makedirs(path, exist_ok=True)
        logger.info(f"Created directory: {path}")
    
    # Track statistics
    stats = defaultdict(int)
    moved_files = []
    
    # Get all files in base directory (not in subdirectories)
    files = [f for f in base_path.iterdir() if f.is_file()]
    logger.info(f"Found {len(files)} files to categorize")
    
    for filepath in files:
        filename = filepath.name
        category = categorize_file(filename)
        dest_dir = CATEGORIES[category]
        dest_path = Path(dest_dir) / filename
        
        stats[category] += 1
        
        if dry_run:
            logger.debug(f"Would move: {filename} -> {category}/")
        else:
            try:
                shutil.move(str(filepath), str(dest_path))
                moved_files.append((filename, category))
            except Exception as e:
                logger.error(f"Failed to move {filename}: {e}")
    
    # Print summary
    logger.info("=" * 60)
    logger.info("Migration Summary")
    logger.info("=" * 60)
    for category, count in sorted(stats.items(), key=lambda x: -x[1]):
        status = "(would move)" if dry_run else "(moved)"
        logger.info(f"  {category}: {count} files {status}")
    logger.info(f"  Total: {sum(stats.values())} files")
    logger.info("=" * 60)
    
    return stats


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Migrate cache files to categories')
    parser.add_argument('--dry-run', action='store_true', 
                        help='Only show what would be done')
    args = parser.parse_args()
    
    if args.dry_run:
        logger.info("DRY RUN - No files will be moved")
    
    migrate_cache(dry_run=args.dry_run)


if __name__ == "__main__":
    main()
