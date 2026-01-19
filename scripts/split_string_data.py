#!/usr/bin/env python
"""
Split STRING Database File by Species (Single-Pass Streaming)

Uses a single-pass streaming approach with file handles management:
- Read input line by line
- Keep output file handles open (up to a limit)
- Use LRU cache to manage file handles for memory efficiency

Usage:
    python split_string_data.py
"""

import os
import sys
import gzip
import hashlib
import argparse
import logging
from pathlib import Path
from time import time
from collections import OrderedDict

# Configuration
DEFAULT_INPUT = "/GenSIvePFS/users/data/STRING/protein.links.detailed.v12.0.txt.gz"
DEFAULT_OUTPUT = "/GenSIvePFS/users/data/pypath_cache/string"
STRING_VERSION = "v12.0"
MAX_OPEN_FILES = 1000  # Maximum number of open file handles

# Setup logging with forced flush
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
# Force unbuffered output
for handler in logging.root.handlers:
    handler.flush = sys.stdout.flush
logger = logging.getLogger(__name__)


def generate_pypath_hash(tax_id: str, link_type: str = "links") -> str:
    """Generate a hash similar to what pypath uses for cache filenames."""
    url = f"https://stringdb-downloads.org/download/protein.{link_type}.detailed.{STRING_VERSION}/{tax_id}.protein.{link_type}.detailed.{STRING_VERSION}.txt.gz"
    return hashlib.md5(url.encode()).hexdigest()


def get_tax_id(protein_id: str) -> str:
    """Extract tax_id from protein ID (format: taxid.protein_name)."""
    return protein_id.split('.')[0]


class LRUFileHandleCache:
    """
    LRU cache for file handles with automatic closing of old handles.
    """
    
    def __init__(self, max_handles: int, output_dir: Path, header: str, buffer_size: int = 10000):
        self.max_handles = max_handles
        self.output_dir = output_dir
        self.header = header
        self.handles = OrderedDict()  # tax_id -> file handle
        self.line_counts = {}  # tax_id -> number of lines written
        self.buffers = {} # tax_id -> list of lines
        self.buffer_size = buffer_size
        
    def get_handle(self, tax_id: str):
        """Get or create file handle for a tax_id."""
        if tax_id in self.handles:
            # Move to end (most recently used)
            self.handles.move_to_end(tax_id)
            return self.handles[tax_id]
        
        # Create new handle
        file_hash = generate_pypath_hash(tax_id)
        filename = f"{file_hash}-{tax_id}.protein.links.detailed.{STRING_VERSION}.txt.gz"
        filepath = self.output_dir / filename
        
        # Close oldest handle if at limit
        if len(self.handles) >= self.max_handles:
            oldest_tax_id, oldest_handle = self.handles.popitem(last=False)
            self._flush_buffer(oldest_tax_id) # Write remaining lines before closing
            oldest_handle.close()
            # Clean up buffer for closed handle if exists (should be empty after flush)
            if oldest_tax_id in self.buffers:
                del self.buffers[oldest_tax_id]
        
        # Open new handle with fast compression
        handle = gzip.open(str(filepath), 'wt', encoding='utf-8', compresslevel=1)
        handle.write(self.header + '\n')
        self.handles[tax_id] = handle
        self.line_counts[tax_id] = 0
        self.buffers[tax_id] = []
        
        return handle
    
    def write_line(self, tax_id: str, line: str):
        """Buffer line and write if buffer full."""
        # Ensure buffer exists
        if tax_id not in self.buffers:
            self.buffers[tax_id] = []
            
        self.buffers[tax_id].append(line)
        self.line_counts[tax_id] = self.line_counts.get(tax_id, 0) + 1
        
        # Flush if buffer is large enough
        if len(self.buffers[tax_id]) >= self.buffer_size:
            self._flush_buffer(tax_id)
            
    def _flush_buffer(self, tax_id: str):
        """Write buffered lines to file."""
        if tax_id in self.buffers and self.buffers[tax_id]:
            handle = self.get_handle(tax_id)
            handle.write('\n'.join(self.buffers[tax_id]) + '\n')
            self.buffers[tax_id] = []
    
    def close_all(self):
        """Flush all buffers and close all open handles."""
        # Flush all pending buffers first
        for tax_id in list(self.buffers.keys()):
            self._flush_buffer(tax_id)
            
        for handle in self.handles.values():
            try:
                handle.close()
            except:
                pass
        self.handles.clear()
        self.buffers.clear()
    
    def get_stats(self) -> dict:
        """Get statistics."""
        return {
            'total_species': len(self.line_counts),
            'total_lines': sum(self.line_counts.values()),
            'open_handles': len(self.handles)
        }


def split_string_file(
    input_file: str,
    output_dir: str,
    skip_existing: bool = True,
    max_handles: int = MAX_OPEN_FILES,
    start_line: int = 0
):
    """
    Split the large STRING file into per-species files using streaming.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Find existing cache files to skip
    existing_tax_ids = set()
    if skip_existing:
        for f in output_path.glob("*-*.protein.links.detailed.*.txt.gz"):
            parts = f.name.split('-', 1)
            if len(parts) == 2:
                tax_part = parts[1].split('.')[0]
                existing_tax_ids.add(tax_part)
        logger.info(f"Found {len(existing_tax_ids)} existing species in cache, will skip them")
    
    logger.info(f"Input file: {input_file}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Max open file handles: {max_handles}")
    if start_line > 0:
        logger.info(f"Resuming from line: {start_line:,}")
    logger.info("=" * 80)
    logger.info("Streaming processing...")
    logger.info("=" * 80)
    
    start_time = time()
    last_log_time = start_time
    line_count = 0
    skipped_count = 0
    new_species = set()
    
    with gzip.open(input_file, 'rt', encoding='utf-8') as f:
        header = f.readline().strip()
        logger.info(f"Header: {header[:60]}...")
        
        # Fast-forward if needed
        if start_line > 0:
            logger.info(f"Fast-forwarding {start_line:,} lines...")
            from itertools import islice
            # Consume lines efficiently
            # We subtract 1 because header is already read
            lines_to_skip = start_line
            
            # Using loop with islice is memory efficient
            # Skipping in chunks to allow progress logging
            chunk_size = 1_000_000
            skipped_so_far = 0
            
            while skipped_so_far < lines_to_skip:
                remaining = lines_to_skip - skipped_so_far
                current_chunk = min(remaining, chunk_size)
                
                # Consume chunk
                for _ in islice(f, current_chunk):
                    pass
                    
                skipped_so_far += current_chunk
                line_count += current_chunk
                
                # Log progress
                if time() - last_log_time >= 10:
                    rate = skipped_so_far / (time() - start_time)
                    logger.info(
                        f"Skipping: {skipped_so_far/1e6:.1f}M / {lines_to_skip/1e6:.1f}M lines, "
                        f"{rate/1e6:.1f}M lines/sec"
                    )
                    sys.stdout.flush()
                    last_log_time = time()
            
            logger.info("Fast-forward complete. Resuming processing...")
            last_log_time = time()

        # Initialize file handle cache
        cache = LRUFileHandleCache(max_handles, output_path, header)
        
        try:
            for line in f:
                line_count += 1
                
                parts = line.split()
                if len(parts) < 2:
                    continue
                
                tax_id = get_tax_id(parts[0])
                
                # Log progress every 30 seconds (before skip check to ensure logging)
                current_time = time()
                if current_time - last_log_time >= 30:
                    elapsed = current_time - start_time
                    rate = line_count / elapsed
                    logger.info(
                        f"Progress: {line_count/1e6:.1f}M lines, "
                        f"{len(new_species)} new species, "
                        f"{skipped_count/1e6:.1f}M skipped, "
                        f"{rate/1000:.1f}K lines/sec"
                    )
                    sys.stdout.flush()
                    last_log_time = current_time
                
                # Skip if already cached
                if tax_id in existing_tax_ids:
                    skipped_count += 1
                    continue
                
                # Track new species
                new_species.add(tax_id)
                
                # Write to output
                cache.write_line(tax_id, line)
        
        finally:
            cache.close_all()
    
    total_time = time() - start_time
    stats = cache.get_stats()
    
    logger.info("=" * 80)
    logger.info("Summary")
    logger.info("=" * 80)
    logger.info(f"Total time: {total_time/60:.1f} minutes")
    logger.info(f"Lines processed: {line_count:,}")
    logger.info(f"Lines skipped (cached): {skipped_count:,}")
    logger.info(f"New species written: {len(new_species)}")
    logger.info(f"Total interactions written: {stats['total_lines']:,}")
    logger.info(f"Output directory: {output_dir}")
    logger.info("=" * 80)
    
    # Update progress file
    progress_file = output_path / "string_download_progress.json"
    try:
        import json
        progress = {'completed': [], 'total_completed': 0}
        if progress_file.exists():
            with open(progress_file, 'r') as pf:
                progress = json.load(pf)
        
        completed_set = set(progress.get('completed', []))
        completed_set.update(new_species)
        
        progress['completed'] = list(completed_set)
        progress['total_completed'] = len(completed_set)
        progress['timestamp'] = time()
        
        with open(progress_file, 'w') as pf:
            json.dump(progress, pf, indent=2)
        
        logger.info(f"Updated progress file: {len(completed_set)} total species")
    except Exception as e:
        logger.warning(f"Could not update progress file: {e}")


def main():
    parser = argparse.ArgumentParser(
        description='Split STRING database file by species (streaming)'
    )
    parser.add_argument(
        '--input', '-i', type=str, default=DEFAULT_INPUT,
        help=f'Input file (default: {DEFAULT_INPUT})'
    )
    parser.add_argument(
        '--output', '-o', type=str, default=DEFAULT_OUTPUT,
        help=f'Output directory (default: {DEFAULT_OUTPUT})'
    )
    parser.add_argument(
        '--no-skip', action='store_true',
        help='Do not skip existing species'
    )
    parser.add_argument(
        '--max-handles', type=int, default=MAX_OPEN_FILES,
        help=f'Max open file handles (default: {MAX_OPEN_FILES})'
    )
    parser.add_argument(
        '--start-line', type=int, default=0,
        help='Line number to resume from (skips first N lines)'
    )
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)
    
    split_string_file(
        input_file=args.input,
        output_dir=args.output,
        skip_existing=not args.no_skip,
        max_handles=args.max_handles,
        start_line=args.start_line
    )


if __name__ == "__main__":
    main()
