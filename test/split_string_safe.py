#!/usr/bin/env python
"""
Split STRING Database File by Species (SAFE Mode)

Combines pigz speed with streaming writes - WILL NOT OOM:
- Uses pigz for parallel gzip decompression (fast)
- LRU file handle cache (memory efficient)
- Writes to disk immediately (no full memory storage)
- Parallel chunk processing

Usage:
    python split_string_safe.py --workers 32
"""

import os
import sys
import gzip
import hashlib
import argparse
import logging
import json
import subprocess
from pathlib import Path
from time import time
from itertools import islice
from multiprocessing import Pool, cpu_count
from collections import defaultdict, OrderedDict

# Memory monitoring
try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False

def get_memory_usage():
    """Get current memory usage in GB."""
    if HAS_PSUTIL:
        process = psutil.Process(os.getpid())
        mem_gb = process.memory_info().rss / (1024**3)
        return f"{mem_gb:.1f}GB"
    return "N/A"

# Configuration
DEFAULT_INPUT = "/GenSIvePFS/users/data/STRING/protein.links.detailed.v12.0.txt.gz"
DEFAULT_OUTPUT = "/GenSIvePFS/users/data/pypath_cache/string"
STRING_VERSION = "v12.0"
CHUNK_SIZE = 1_000_000   # Lines per chunk (default; keep memory bounded)
MAX_OPEN_FILES = 1000    # Maximum open file handles
BUFFER_SIZE = 10000      # Lines to buffer per species before writing

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
for handler in logging.root.handlers:
    handler.flush = sys.stdout.flush
logger = logging.getLogger(__name__)


def generate_pypath_hash(tax_id: str) -> str:
    """Generate hash matching pypath's cache naming."""
    url = f"https://stringdb-downloads.org/download/protein.links.detailed.{STRING_VERSION}/{tax_id}.protein.links.detailed.{STRING_VERSION}.txt.gz"
    return hashlib.md5(url.encode()).hexdigest()


def get_tax_id(protein_id: str) -> str:
    """Extract tax_id from protein ID."""
    return protein_id.split('.')[0]


def check_pigz():
    """Check if pigz is available."""
    try:
        result = subprocess.run(['pigz', '--version'], capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False


class LRUFileHandleCache:
    """
    LRU cache for file handles with automatic closing of old handles.
    Memory-efficient - writes to disk immediately instead of storing in memory.
    """
    
    def __init__(self, max_handles: int, output_dir: Path, header: str, buffer_size: int = BUFFER_SIZE):
        self.max_handles = max_handles
        self.output_dir = output_dir
        self.header = header
        self.handles = OrderedDict()  # tax_id -> file handle
        self.line_counts = {}  # tax_id -> number of lines written
        self.buffers = {}  # tax_id -> list of lines (small buffer)
        self.buffer_size = buffer_size
        
    def _open_handle(self, tax_id: str):
        """Open a handle in append mode, writing header if file is new/empty."""
        file_hash = generate_pypath_hash(tax_id)
        filename = f"{file_hash}-{tax_id}.protein.links.detailed.{STRING_VERSION}.txt.gz"
        filepath = self.output_dir / filename

        is_new = (not filepath.exists()) or filepath.stat().st_size == 0
        handle = gzip.open(str(filepath), 'at', encoding='utf-8', compresslevel=1)
        if is_new:
            handle.write(self.header + '\n')
        return handle

    def get_handle(self, tax_id: str):
        """Get or create file handle for a tax_id."""
        if tax_id in self.handles:
            self.handles.move_to_end(tax_id)
            return self.handles[tax_id]
        
        # Close oldest handle if at limit
        if len(self.handles) >= self.max_handles:
            oldest_tax_id, oldest_handle = self.handles.popitem(last=False)
            self._flush_buffer(oldest_tax_id, handle=oldest_handle, allow_open=False)
            oldest_handle.close()
            if oldest_tax_id in self.buffers:
                del self.buffers[oldest_tax_id]
        
        # Open new handle with fast compression (append; avoid truncation)
        handle = self._open_handle(tax_id)
        self.handles[tax_id] = handle
        if tax_id not in self.line_counts:
            self.line_counts[tax_id] = 0
        self.buffers[tax_id] = []
        
        return handle
    
    def write_lines(self, tax_id: str, lines: list):
        """Write multiple lines - buffers then flushes to disk."""
        if tax_id not in self.buffers:
            self.buffers[tax_id] = []
        
        self.buffers[tax_id].extend(lines)
        self.line_counts[tax_id] = self.line_counts.get(tax_id, 0) + len(lines)
        
        # Flush if buffer is large enough - KEY: writes to disk, frees memory
        if len(self.buffers[tax_id]) >= self.buffer_size:
            self._flush_buffer(tax_id)
            
    def _flush_buffer(self, tax_id: str, handle=None, allow_open: bool = True):
        """Write buffered lines to file and FREE memory."""
        if tax_id in self.buffers and self.buffers[tax_id]:
            if handle is None:
                if tax_id in self.handles:
                    handle = self.handles[tax_id]
                elif allow_open:
                    handle = self._open_handle(tax_id)
                    self.handles[tax_id] = handle
                else:
                    return
            handle.write(''.join(self.buffers[tax_id]))
            self.buffers[tax_id] = []  # Free memory!
    
    def flush_all(self):
        """Flush all buffers."""
        for tax_id in list(self.buffers.keys()):
            self._flush_buffer(tax_id)
    
    def close_all(self):
        """Flush all buffers and close all handles."""
        for tax_id in list(self.buffers.keys()):
            self._flush_buffer(tax_id)

        for tax_id, handle in list(self.handles.items()):
            try:
                handle.close()
                # Delete empty files (only header, no data)
                if self.line_counts.get(tax_id, 0) == 0:
                    file_hash = generate_pypath_hash(tax_id)
                    filename = f"{file_hash}-{tax_id}.protein.links.detailed.{STRING_VERSION}.txt.gz"
                    filepath = self.output_dir / filename
                    if filepath.exists():
                        filepath.unlink()
                        logger.debug(f"Deleted empty file for tax_id {tax_id}")
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


def process_chunk(args):
    """Process a chunk of lines and group by tax_id."""
    chunk_lines, chunk_id, skip_tax_ids = args

    grouped = defaultdict(list)
    processed = 0
    skipped = 0

    for line in chunk_lines:
        parts = line.split(None, 1)
        if not parts:
            continue

        tax_id = parts[0].split('.', 1)[0]

        if tax_id in skip_tax_ids:
            skipped += 1
            continue

        grouped[tax_id].append(line)
        processed += 1

    return chunk_id, dict(grouped), processed, skipped


def iter_chunks(lines_iter, chunk_size: int):
    """Yield lists of lines with a bounded size."""
    chunk = []
    for line in lines_iter:
        chunk.append(line)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def is_complete_file(path: Path) -> bool:
    """Return True if gz file has at least header + 1 data line."""
    if not path.exists():
        return False
    try:
        with gzip.open(str(path), 'rt', encoding='utf-8') as fh:
            header = fh.readline()
            if not header:
                return False
            data_line = fh.readline()
            return bool(data_line)
    except Exception:
        return False


def split_string_safe(
    input_file: str,
    output_dir: str,
    workers: int = None,
    skip_existing: bool = True,
    max_handles: int = MAX_OPEN_FILES,
    start_line: int = 0,
    chunk_size: int = CHUNK_SIZE
):
    """
    Split STRING file using SAFE mode:
    - pigz for fast decompression
    - LRU cache for memory efficiency
    - Streaming writes (NO OOM)
    """
    if workers is None:
        workers = cpu_count()

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Check for pigz
    has_pigz = check_pigz()
    if has_pigz:
        logger.info("✓ pigz detected - using parallel decompression")
    else:
        logger.warning("✗ pigz not found - using Python gzip (slower)")

    # Find existing files
    existing_tax_ids = set()
    if skip_existing:
        for f in output_path.glob("*-*.protein.links.detailed.*.txt.gz"):
            parts = f.name.split('-', 1)
            if len(parts) == 2:
                tax_part = parts[1].split('.')[0]
                if is_complete_file(f):
                    existing_tax_ids.add(tax_part)
                else:
                    try:
                        f.unlink()
                        logger.warning(f"Removed incomplete file: {f.name}")
                    except Exception:
                        logger.warning(f"Could not remove incomplete file: {f.name}")
        logger.info(f"Total species to skip: {len(existing_tax_ids)}")

    logger.info(f"Input: {input_file}")
    logger.info(f"Output: {output_dir}")
    logger.info(f"Workers: {workers}")
    logger.info(f"Chunk size: {chunk_size:,} lines")
    logger.info(f"Max open files: {max_handles}")
    logger.info(f"Mode: SAFE (streaming writes, NO OOM risk)")
    if start_line > 0:
        logger.info(f"Fast-forward to line: {start_line:,}")
    logger.info("=" * 80)

    start_time = time()
    last_log_time = start_time

    # Open input with pigz or Python gzip
    if has_pigz:
        try:
            proc = subprocess.Popen(
                ['pigz', '-dc', '-p', str(min(workers, 16)), input_file],
                stdout=subprocess.PIPE,
                bufsize=1024*1024*256
            )
            lines_iter = (line.decode('utf-8') for line in proc.stdout)
        except Exception as e:
            logger.warning(f"pigz failed ({e}); falling back to Python gzip")
            has_pigz = False
    if not has_pigz:
        f = gzip.open(input_file, 'rt', encoding='utf-8')
        lines_iter = iter(f)

    # Read header
    header = next(lines_iter).strip()
    logger.info(f"Header: {header[:60]}...")

    # Fast-forward if needed
    if start_line > 0:
        logger.info(f"Fast-forwarding {start_line:,} lines...")
        ff_start = time()
        ff_last_log = ff_start
        skipped_so_far = 0
        chunk_size_ff = 10_000_000
        
        eof_reached = False
        while skipped_so_far < start_line:
            remaining = start_line - skipped_so_far
            current_chunk = min(remaining, chunk_size_ff)
            
            # Consume lines efficiently
            consumed = 0
            for _ in islice(lines_iter, current_chunk):
                consumed += 1
            if consumed < current_chunk:
                eof_reached = True
                skipped_so_far += consumed
                break
            
            skipped_so_far += current_chunk
            
            # Log progress every 30 seconds
            if time() - ff_last_log >= 30:
                elapsed = time() - ff_start
                rate = skipped_so_far / elapsed if elapsed > 0 else 0
                eta = (start_line - skipped_so_far) / rate if rate > 0 else 0
                logger.info(
                    f"Fast-forward: {skipped_so_far/1e9:.2f}B / {start_line/1e9:.2f}B lines, "
                    f"{rate/1e6:.1f}M lines/sec, ETA: {eta/60:.1f} min"
                )
                sys.stdout.flush()
                ff_last_log = time()
        
        ff_time = time() - ff_start
        if eof_reached and skipped_so_far < start_line:
            logger.warning(
                f"Reached EOF at line {skipped_so_far:,} before start_line {start_line:,}"
            )
            # Close input and exit early
            if has_pigz:
                proc.wait()
            else:
                try:
                    f.close()
                except Exception:
                    pass
            return
        logger.info(f"Fast-forward complete in {ff_time/60:.1f} min")

    # Initialize LRU file handle cache - KEY for memory efficiency
    cache = LRUFileHandleCache(max_handles, output_path, header)
    
    total_lines = 0
    total_skipped = 0
    chunks_processed = 0
    new_species = set()

    logger.info("Processing with streaming writes (memory safe)...")
    logger.info("=" * 80)

    try:
        with Pool(workers) as pool:
            for cid, grouped, processed, skipped in pool.imap_unordered(
                process_chunk,
                ((chunk, idx, existing_tax_ids)
                 for idx, chunk in enumerate(iter_chunks(lines_iter, chunk_size))),
                chunksize=1
            ):
                for tax_id, lines in grouped.items():
                    cache.write_lines(tax_id, lines)
                    new_species.add(tax_id)
                total_lines += processed
                total_skipped += skipped
                chunks_processed += 1

                if chunks_processed % 10 == 0:
                    cache.flush_all()

                if time() - last_log_time >= 30:
                    mem_usage = get_memory_usage()
                    stats = cache.get_stats()
                    logger.info(
                        f"Chunks: {chunks_processed}, Lines: {total_lines/1e6:.1f}M, "
                        f"Skipped: {total_skipped/1e6:.1f}M, "
                        f"Species: {len(new_species)}, "
                        f"Written: {stats['total_lines']/1e6:.1f}M, "
                        f"Memory: {mem_usage}"
                    )
                    sys.stdout.flush()
                    last_log_time = time()

    finally:
        cache.close_all()
        
        # Close input
        if has_pigz:
            proc.wait()
        else:
            try:
                f.close()
            except:
                pass

    total_time = time() - start_time
    stats = cache.get_stats()

    # Save progress
    progress_file = output_path / "split_progress.json"
    try:
        progress_data = {
            'completed': list(new_species),
            'total_completed': len(new_species),
            'total_lines': total_lines,
            'timestamp': time(),
            'final': True
        }
        with open(progress_file, 'w') as pf:
            json.dump(progress_data, pf, indent=2)
    except Exception as e:
        logger.warning(f"Could not save progress: {e}")

    # Summary
    logger.info("=" * 80)
    logger.info("✅ SAFE Mode Complete (No OOM)")
    logger.info("=" * 80)
    logger.info(f"Total time: {total_time/60:.1f} min")
    logger.info(f"Lines processed: {total_lines:,}")
    logger.info(f"Lines skipped: {total_skipped:,}")
    logger.info(f"Species written: {len(new_species)}")
    logger.info(f"Total interactions: {stats['total_lines']:,}")
    logger.info(f"Throughput: {(total_lines+total_skipped)/total_time/1e6:.2f}M lines/sec")
    logger.info(f"Output directory: {output_dir}")
    logger.info("=" * 80)


def main():
    parser = argparse.ArgumentParser(
        description='Split STRING file (SAFE mode - no OOM risk)'
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
        '--workers', '-w', type=int, default=None,
        help=f'Number of workers (default: CPU count = {cpu_count()})'
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
        help='Line number to resume from (default: 0)'
    )
    parser.add_argument(
        '--chunk-size', type=int, default=CHUNK_SIZE,
        help=f'Lines per chunk (default: {CHUNK_SIZE:,})'
    )

    args = parser.parse_args()

    if not os.path.exists(args.input):
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)

    split_string_safe(
        input_file=args.input,
        output_dir=args.output,
        workers=args.workers,
        skip_existing=not args.no_skip,
        max_handles=args.max_handles,
        start_line=args.start_line,
        chunk_size=args.chunk_size
    )


if __name__ == "__main__":
    main()
