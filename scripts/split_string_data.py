#!/usr/bin/env python
"""
Split STRING Database File by Species

This script splits the large STRING protein.links.detailed file into
per-species files that can be used as pypath cache.

Usage:
    python split_string_data.py --input /path/to/protein.links.detailed.v12.0.txt.gz
"""

import os
import sys
import gzip
import hashlib
import argparse
import logging
from pathlib import Path
from collections import defaultdict
from time import time
from tqdm import tqdm

# Configuration
DEFAULT_INPUT = "/GenSIvePFS/users/data/STRING/protein.links.detailed.v12.0.txt.gz"
DEFAULT_OUTPUT = "/GenSIvePFS/users/data/pypath_cache"
STRING_VERSION = "v12.0"

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def generate_pypath_hash(tax_id: str, link_type: str = "links") -> str:
    """
    Generate a hash similar to what pypath uses for cache filenames.
    
    The hash is based on the URL pattern pypath would use to download the file.
    """
    # pypath uses URLs like:
    # https://stringdb-downloads.org/download/protein.links.detailed.v12.0/{taxid}.protein.links.detailed.v12.0.txt.gz
    url = f"https://stringdb-downloads.org/download/protein.{link_type}.detailed.{STRING_VERSION}/{tax_id}.protein.{link_type}.detailed.{STRING_VERSION}.txt.gz"
    return hashlib.md5(url.encode()).hexdigest()


def count_lines(filepath: str) -> int:
    """Count lines in a gzip file for progress bar."""
    logger.info("Counting lines in input file (this may take a while)...")
    count = 0
    try:
        with gzip.open(filepath, 'rt', encoding='utf-8') as f:
            for _ in f:
                count += 1
        return count
    except Exception as e:
        logger.warning(f"Could not count lines: {e}, will use estimate")
        # Estimate based on file size (roughly 100 bytes per line)
        file_size = os.path.getsize(filepath)
        # Compressed ratio is roughly 10:1
        return (file_size * 10) // 100


def get_tax_id(protein_id: str) -> str:
    """Extract tax_id from protein ID (format: taxid.protein_name)."""
    return protein_id.split('.')[0]


def split_string_file(
    input_file: str,
    output_dir: str,
    batch_size: int = 1000000,
    skip_existing: bool = True
):
    """
    Split the large STRING file into per-species files.
    
    Args:
        input_file: Path to protein.links.detailed.v12.0.txt.gz
        output_dir: Directory to write output files
        batch_size: Number of lines to process before writing to disk
        skip_existing: Skip species that already have cache files
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Find existing cache files to skip
    existing_tax_ids = set()
    if skip_existing:
        for f in output_path.glob("*-*.protein.links.detailed.*.txt.gz"):
            # Extract tax_id from filename like: hash-taxid.protein.links...
            parts = f.name.split('-', 1)
            if len(parts) == 2:
                tax_part = parts[1].split('.')[0]
                existing_tax_ids.add(tax_part)
        logger.info(f"Found {len(existing_tax_ids)} existing species in cache, will skip them")
    
    logger.info(f"Input file: {input_file}")
    logger.info(f"Output directory: {output_dir}")
    logger.info("=" * 80)
    logger.info("Phase 1: Reading and grouping data by species...")
    logger.info("=" * 80)
    
    # Count lines for progress
    total_lines = count_lines(input_file)
    logger.info(f"Total lines to process: {total_lines:,}")
    
    # Group data by tax_id
    species_data = defaultdict(list)
    header = None
    skipped_lines = 0
    processed_lines = 0
    
    start_time = time()
    
    with gzip.open(input_file, 'rt', encoding='utf-8') as f:
        # Read header
        header = f.readline().strip()
        
        with tqdm(total=total_lines - 1, desc="Reading", unit=" lines") as pbar:
            for line in f:
                processed_lines += 1
                pbar.update(1)
                
                line = line.strip()
                if not line:
                    continue
                
                # Extract tax_id from protein1
                parts = line.split('\t') if '\t' in line else line.split()
                if len(parts) < 2:
                    continue
                
                protein1 = parts[0]
                tax_id = get_tax_id(protein1)
                
                # Skip if already in cache
                if skip_existing and tax_id in existing_tax_ids:
                    skipped_lines += 1
                    continue
                
                species_data[tax_id].append(line)
                
                # Periodically log progress
                if processed_lines % 10000000 == 0:
                    elapsed = time() - start_time
                    rate = processed_lines / elapsed
                    logger.info(f"Progress: {processed_lines:,} lines, "
                               f"{len(species_data)} species, "
                               f"{rate:.0f} lines/sec")
    
    read_time = time() - start_time
    logger.info(f"\nPhase 1 complete in {read_time/60:.1f} minutes")
    logger.info(f"Processed {processed_lines:,} lines")
    logger.info(f"Skipped {skipped_lines:,} lines (already cached)")
    logger.info(f"Found {len(species_data)} new species to write")
    
    if not species_data:
        logger.info("No new species to write. All data already cached!")
        return
    
    logger.info("=" * 80)
    logger.info("Phase 2: Writing per-species files...")
    logger.info("=" * 80)
    
    write_start = time()
    written_species = 0
    total_interactions = 0
    
    # Sort species by number of interactions (largest first for better progress feedback)
    sorted_species = sorted(species_data.items(), key=lambda x: len(x[1]), reverse=True)
    
    with tqdm(total=len(sorted_species), desc="Writing", unit=" species") as pbar:
        for tax_id, lines in sorted_species:
            try:
                # Generate filename
                file_hash = generate_pypath_hash(tax_id)
                filename = f"{file_hash}-{tax_id}.protein.links.detailed.{STRING_VERSION}.txt.gz"
                filepath = output_path / filename
                
                # Write file
                with gzip.open(str(filepath), 'wt', encoding='utf-8') as out:
                    out.write(header + '\n')
                    for line in lines:
                        out.write(line + '\n')
                
                written_species += 1
                total_interactions += len(lines)
                
                pbar.update(1)
                pbar.set_postfix({
                    'species': tax_id,
                    'interactions': f"{len(lines):,}"
                })
                
            except Exception as e:
                logger.error(f"Failed to write species {tax_id}: {e}")
    
    write_time = time() - write_start
    total_time = time() - start_time
    
    logger.info("=" * 80)
    logger.info("Summary")
    logger.info("=" * 80)
    logger.info(f"Total time: {total_time/60:.1f} minutes")
    logger.info(f"  - Reading: {read_time/60:.1f} minutes")
    logger.info(f"  - Writing: {write_time/60:.1f} minutes")
    logger.info(f"Species written: {written_species}")
    logger.info(f"Total interactions: {total_interactions:,}")
    logger.info(f"Output directory: {output_dir}")
    logger.info("=" * 80)


def main():
    parser = argparse.ArgumentParser(
        description='Split STRING database file by species',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--input', '-i', type=str, default=DEFAULT_INPUT,
        help=f'Input protein.links.detailed file (default: {DEFAULT_INPUT})'
    )
    parser.add_argument(
        '--output', '-o', type=str, default=DEFAULT_OUTPUT,
        help=f'Output directory for cache files (default: {DEFAULT_OUTPUT})'
    )
    parser.add_argument(
        '--no-skip', action='store_true',
        help='Do not skip species that already have cache files'
    )
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)
    
    split_string_file(
        input_file=args.input,
        output_dir=args.output,
        skip_existing=not args.no_skip
    )


if __name__ == "__main__":
    main()
