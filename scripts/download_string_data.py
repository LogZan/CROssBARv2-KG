#!/usr/bin/env python
"""
STRING Database Downloader with Enhanced Resilience

This script downloads STRING protein interaction data for all species
with robust error handling, intelligent retry mechanism, and resume capability.

Key improvements:
- Connection pooling with requests Session
- Configurable timeouts (connect + read)
- Adaptive rate limiting to avoid server throttling
- Automatic cleanup of corrupted cache files
- Smarter exponential backoff with jitter
- Periodic heartbeat saves
- Better error classification and handling

Usage:
    python download_string_data.py [--workers 4] [--max-retries 5]
"""

import os
import sys
import json
import argparse
import logging
import signal
import random
import socket
import glob
import gzip
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed, TimeoutError as FuturesTimeoutError
from threading import Lock, Event
from time import time, sleep
from typing import Set, Dict, Optional, List
from contextlib import contextmanager
from functools import wraps
import urllib3

# Suppress SSL warnings if needed
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# Setup paths
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pypath.inputs import string
from pypath.share import curl, settings
from tqdm import tqdm

# Configuration
CACHE_DIR = "/GenSIvePFS/users/data/pypath_cache"
PROGRESS_FILE = CACHE_DIR + "/string_download_progress.json"
FAILED_FILE = CACHE_DIR + "/string_download_failed.json"
LOG_FILE = CACHE_DIR + "/string_download.log"

# Network settings
DEFAULT_CONNECT_TIMEOUT = 30  # seconds
DEFAULT_READ_TIMEOUT = 300    # 5 minutes for large files
DEFAULT_WORKERS = 4           # Conservative default to avoid throttling
DEFAULT_MAX_RETRIES = 5
DEFAULT_RETRY_BASE_DELAY = 5  # seconds
MAX_RETRY_DELAY = 300         # 5 minutes max delay

# Ensure cache directory exists
os.makedirs(CACHE_DIR, exist_ok=True)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - [%(threadName)s] - %(message)s',
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Thread-safe locks
progress_lock = Lock()
failed_lock = Lock()
save_lock = Lock()

# Graceful shutdown flag
shutdown_event = Event()


def signal_handler(signum, frame):
    """Handle shutdown signals gracefully."""
    logger.warning("Received shutdown signal. Finishing current downloads...")
    shutdown_event.set()


# Register signal handlers
signal.signal(signal.SIGINT, signal_handler)
signal.signal(signal.SIGTERM, signal_handler)


class NetworkError(Exception):
    """Network-related errors that might be recoverable."""
    pass


class ServerError(Exception):
    """Server-side errors (5xx) that might be recoverable with backoff."""
    pass


class ClientError(Exception):
    """Client-side errors (4xx) that are unlikely to be recoverable."""
    pass


class DownloadStats:
    """Thread-safe statistics tracker."""
    
    def __init__(self):
        self._lock = Lock()
        self.success_count = 0
        self.failed_count = 0
        self.retry_count = 0
        self.total_interactions = 0
        self.start_time = time()
    
    def record_success(self, interactions: int = 0):
        with self._lock:
            self.success_count += 1
            self.total_interactions += interactions
    
    def record_failure(self):
        with self._lock:
            self.failed_count += 1
    
    def record_retry(self):
        with self._lock:
            self.retry_count += 1
    
    def get_summary(self) -> dict:
        with self._lock:
            elapsed = time() - self.start_time
            return {
                'success': self.success_count,
                'failed': self.failed_count,
                'retries': self.retry_count,
                'total_interactions': self.total_interactions,
                'elapsed_seconds': elapsed,
                'rate_per_minute': self.success_count / elapsed * 60 if elapsed > 0 else 0
            }


def calculate_backoff(retry_count: int, base_delay: float = DEFAULT_RETRY_BASE_DELAY) -> float:
    """
    Calculate exponential backoff with jitter.
    
    Uses decorrelated jitter algorithm for better distribution.
    """
    # Exponential backoff
    delay = base_delay * (2 ** retry_count)
    
    # Add jitter (up to 25% variation)
    jitter = delay * 0.25 * random.random()
    delay = delay + jitter
    
    # Cap at maximum delay
    return min(delay, MAX_RETRY_DELAY)


def retry_with_backoff(max_retries: int = DEFAULT_MAX_RETRIES, 
                       base_delay: float = DEFAULT_RETRY_BASE_DELAY,
                       recoverable_exceptions: tuple = (NetworkError, ServerError, TimeoutError)):
    """
    Decorator for retry with exponential backoff.
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            last_exception = None
            
            for attempt in range(max_retries + 1):
                if shutdown_event.is_set():
                    raise InterruptedError("Shutdown requested")
                
                try:
                    return func(*args, **kwargs)
                except recoverable_exceptions as e:
                    last_exception = e
                    
                    if attempt < max_retries:
                        delay = calculate_backoff(attempt, base_delay)
                        logger.warning(
                            f"Attempt {attempt + 1}/{max_retries + 1} failed: {e}. "
                            f"Retrying in {delay:.1f}s..."
                        )
                        sleep(delay)
                    else:
                        logger.error(f"All {max_retries + 1} attempts failed: {e}")
                        raise
                except ClientError:
                    # Don't retry client errors
                    raise
            
            raise last_exception
        return wrapper
    return decorator


class StringDownloader:
    """STRING database downloader with enhanced resilience."""
    
    def __init__(
        self, 
        cache_dir: str, 
        max_workers: int = DEFAULT_WORKERS, 
        max_retries: int = DEFAULT_MAX_RETRIES,
        connect_timeout: int = DEFAULT_CONNECT_TIMEOUT,
        read_timeout: int = DEFAULT_READ_TIMEOUT,
        rate_limit_delay: float = 0.5
    ):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.max_workers = max_workers
        self.max_retries = max_retries
        self.connect_timeout = connect_timeout
        self.read_timeout = read_timeout
        self.rate_limit_delay = rate_limit_delay
        
        # Setup pypath cache and timeouts
        settings.setup(cachedir=str(self.cache_dir))
        
        # Configure curl module for better resilience
        self._configure_curl()
        
        # Load progress
        self.completed: Set[str] = self._load_progress()
        self.failed: Dict[str, dict] = self._load_failed()
        
        # Statistics
        self.stats = DownloadStats()
        
        # Last save time for periodic saves
        self.last_save_time = time()
        self.save_interval = 60  # Save every 60 seconds
        
        logger.info(f"Initialized downloader with cache dir: {self.cache_dir}")
        logger.info(f"Already completed: {len(self.completed)} species")
        logger.info(f"Previously failed: {len(self.failed)} species")
        logger.info(f"Settings: workers={max_workers}, retries={max_retries}, "
                   f"connect_timeout={connect_timeout}s, read_timeout={read_timeout}s")
    
    def _configure_curl(self):
        """Configure curl/requests settings for better resilience."""
        try:
            # Set default socket timeout
            socket.setdefaulttimeout(self.read_timeout)
            
            # Try to configure pypath's curl settings if available
            if hasattr(curl, 'Curl'):
                # These are common settings that might help
                curl.Curl.timeout = self.read_timeout
                curl.Curl.connect_timeout = self.connect_timeout
            
            logger.debug("Configured network settings")
        except Exception as e:
            logger.warning(f"Could not configure some curl settings: {e}")
    
    def _load_progress(self) -> Set[str]:
        """Load list of already completed downloads."""
        if os.path.exists(PROGRESS_FILE):
            try:
                with open(PROGRESS_FILE, 'r') as f:
                    data = json.load(f)
                    return set(data.get('completed', []))
            except Exception as e:
                logger.warning(f"Failed to load progress file: {e}")
        return set()
    
    def _save_progress(self, force: bool = False):
        """Save current progress with rate limiting."""
        current_time = time()
        
        # Only save if forced or enough time has passed
        if not force and (current_time - self.last_save_time) < self.save_interval:
            return
        
        with save_lock:
            try:
                self.last_save_time = current_time
                
                # Save progress
                with open(PROGRESS_FILE, 'w') as f:
                    json.dump({
                        'completed': list(self.completed),
                        'total_completed': len(self.completed),
                        'timestamp': current_time,
                        'stats': self.stats.get_summary()
                    }, f, indent=2)
                
                logger.debug(f"Progress saved: {len(self.completed)} completed")
            except Exception as e:
                logger.error(f"Failed to save progress: {e}")
    
    def _load_failed(self) -> Dict[str, dict]:
        """Load list of failed downloads with error details."""
        if os.path.exists(FAILED_FILE):
            try:
                with open(FAILED_FILE, 'r') as f:
                    data = json.load(f)
                    # Convert old format if needed
                    if data and isinstance(next(iter(data.values())), str):
                        return {k: {'error': v, 'retry_count': 0} for k, v in data.items()}
                    return data
            except Exception as e:
                logger.warning(f"Failed to load failed file: {e}")
        return {}
    
    def _save_failed(self):
        """Save failed downloads."""
        with failed_lock:
            try:
                with open(FAILED_FILE, 'w') as f:
                    json.dump(self.failed, f, indent=2)
            except Exception as e:
                logger.error(f"Failed to save failed list: {e}")
    
    def _mark_completed(self, tax_id: str, interactions_count: int = 0):
        """Mark a tax_id as completed."""
        with progress_lock:
            self.completed.add(tax_id)
            # Remove from failed if it was there
            if tax_id in self.failed:
                with failed_lock:
                    del self.failed[tax_id]
        
        self.stats.record_success(interactions_count)
    
    def _mark_failed(self, tax_id: str, error: str, retry_count: int = 0):
        """Mark a tax_id as failed with details."""
        with failed_lock:
            self.failed[tax_id] = {
                'error': error,
                'retry_count': retry_count,
                'timestamp': time(),
                'error_type': self._classify_error(error)
            }
        self.stats.record_failure()
    
    def _classify_error(self, error_msg: str) -> str:
        """Classify error type for better retry strategies."""
        error_lower = error_msg.lower()
        
        # Corrupted cache/compressed file errors - need to clean cache
        if any(x in error_lower for x in [
            'compressed file ended', 
            'decompressing data',
            'invalid compressed data',
            'crc check failed',
            'truncated',
            'end-of-stream marker'
        ]):
            return 'corrupted_cache'
        elif any(x in error_lower for x in ['timeout', 'timed out']):
            return 'timeout'
        elif any(x in error_lower for x in ['connection', 'connect', 'network']):
            return 'network'
        elif any(x in error_lower for x in ['404', 'not found']):
            return 'not_found'
        elif any(x in error_lower for x in ['429', 'rate limit', 'too many']):
            return 'rate_limit'
        elif any(x in error_lower for x in ['500', '502', '503', '504']):
            return 'server_error'
        else:
            return 'unknown'
    
    def _clear_corrupted_cache(self, tax_id: str) -> int:
        """
        Clear cache files for a specific tax_id.
        
        Since we're only called when a corruption error has already been detected,
        we delete all cache files for this tax_id without further verification.
        
        Args:
            tax_id: NCBI taxonomy ID
            
        Returns:
            Number of files deleted
        """
        deleted = 0
        
        # Scan all gz files in the cache directory and look for ones matching this tax_id
        # Pattern in filename is: {hash}-{taxid}.protein.*.txt.gz
        cache_files = list(self.cache_dir.glob("*.gz"))
        
        for cache_file in cache_files:
            filename = cache_file.name
            # Check if filename contains the tax_id in the expected pattern
            # Example: 3d0563c961eccfbe1a888ad5881aa85a-1648.protein.links.detailed.v12.0.txt.gz
            if f"-{tax_id}.protein." in filename:
                try:
                    os.remove(str(cache_file))
                    deleted += 1
                    logger.info(f"Deleted cache file: {filename}")
                except Exception as e:
                    logger.warning(f"Could not delete cache file {filename}: {e}")
        
        return deleted
        
        return deleted
    
    def _should_retry_error(self, error_type: str, retry_count: int) -> bool:
        """Determine if an error should be retried based on type."""
        # Never retry not_found errors
        if error_type == 'not_found':
            return False
        
        # Rate limit errors need longer backoff
        if error_type == 'rate_limit':
            return retry_count < 3  # Fewer retries for rate limits
        
        # Corrupted cache should definitely be retried after cache cleanup
        if error_type == 'corrupted_cache':
            return retry_count < self.max_retries
        
        return retry_count < self.max_retries
    
    def download_single_species(self, tax_id: str) -> tuple:
        """
        Download STRING data for a single species with robust error handling.
        
        Args:
            tax_id: NCBI taxonomy ID
            
        Returns:
            Tuple of (tax_id, success, error_message, interactions_count)
        """
        if shutdown_event.is_set():
            return (tax_id, False, "Shutdown requested", 0)
        
        retry_count = 0
        last_error = None
        
        while retry_count <= self.max_retries:
            try:
                if shutdown_event.is_set():
                    return (tax_id, False, "Shutdown requested", 0)
                
                logger.debug(f"Downloading tax_id {tax_id} (attempt {retry_count + 1})")
                
                # Add rate limiting between downloads
                if self.rate_limit_delay > 0:
                    sleep(self.rate_limit_delay * random.uniform(0.5, 1.5))
                
                # Call STRING API with timeout awareness
                interactions = list(string.string_links_interactions(
                    ncbi_tax_id=int(tax_id),
                    score_threshold="high_confidence"
                ))
                
                count = len(interactions)
                logger.info(f"✓ Downloaded {count:,} interactions for tax_id {tax_id}")
                return (tax_id, True, None, count)
                
            except KeyboardInterrupt:
                shutdown_event.set()
                return (tax_id, False, "Keyboard interrupt", 0)
                
            except Exception as e:
                error_msg = str(e)
                error_type = self._classify_error(error_msg)
                last_error = error_msg
                
                # For corrupted cache errors, clean the cache before retrying
                if error_type == 'corrupted_cache':
                    deleted = self._clear_corrupted_cache(tax_id)
                    if deleted > 0:
                        logger.info(f"Cleaned {deleted} corrupted cache file(s) for tax_id {tax_id}")
                
                # Check if we should retry
                if self._should_retry_error(error_type, retry_count):
                    retry_count += 1
                    self.stats.record_retry()
                    
                    # Calculate appropriate backoff
                    if error_type == 'rate_limit':
                        delay = calculate_backoff(retry_count, base_delay=30)  # Longer for rate limits
                    elif error_type == 'server_error':
                        delay = calculate_backoff(retry_count, base_delay=10)
                    elif error_type == 'corrupted_cache':
                        delay = calculate_backoff(retry_count, base_delay=2)  # Faster retry after cache cleanup
                    else:
                        delay = calculate_backoff(retry_count)
                    
                    logger.warning(
                        f"✗ Failed tax_id {tax_id} ({error_type}): {error_msg}. "
                        f"Retry {retry_count}/{self.max_retries} in {delay:.1f}s..."
                    )
                    sleep(delay)
                else:
                    logger.error(
                        f"✗ Failed tax_id {tax_id} ({error_type}): {error_msg}. "
                        f"Not retrying ({retry_count} attempts made)"
                    )
                    break
        
        return (tax_id, False, last_error, 0)
    
    def download_all(
        self, 
        skip_completed: bool = True, 
        retry_failed: bool = True,
        priority_list: Optional[list] = None
    ):
        """
        Download STRING data for all species.
        
        Args:
            skip_completed: Skip already downloaded species
            retry_failed: Retry previously failed downloads
            priority_list: Optional list of tax_ids to prioritize
        """
        logger.info("=" * 80)
        logger.info("Starting STRING data download")
        logger.info("=" * 80)
        
        # Get list of all species
        logger.info("Fetching list of all species...")
        try:
            string_species = string.string_species()
            all_tax_ids = list(string_species.keys())
        except Exception as e:
            logger.error(f"Failed to fetch species list: {e}")
            return
        
        # Filter out known problematic tax IDs
        tax_ids_to_skip = {'4565', '8032'}  # Known problematic IDs
        
        # Determine which tax_ids to download
        pending_tax_ids = []
        
        if skip_completed:
            pending_tax_ids = [
                tid for tid in all_tax_ids 
                if tid not in self.completed and tid not in tax_ids_to_skip
            ]
        else:
            pending_tax_ids = [tid for tid in all_tax_ids if tid not in tax_ids_to_skip]
        
        # Add failed ones if retry_failed is True (with better filtering)
        if retry_failed and self.failed:
            # Only retry failed ones that might succeed
            failed_to_retry = [
                tid for tid, info in self.failed.items() 
                if tid not in pending_tax_ids 
                and tid not in tax_ids_to_skip
                and info.get('error_type') not in ['not_found', 'client_error']
            ]
            pending_tax_ids.extend(failed_to_retry)
            logger.info(f"Will retry {len(failed_to_retry)} previously failed downloads")
        
        # Apply priority list if provided
        if priority_list:
            priority_set = set(priority_list)
            pending_tax_ids.sort(key=lambda x: (x not in priority_set, pending_tax_ids.index(x)))
        
        total_species = len(all_tax_ids)
        total_pending = len(pending_tax_ids)
        
        logger.info(f"Total species: {total_species}")
        logger.info(f"Already completed: {len(self.completed)}")
        logger.info(f"Pending downloads: {total_pending}")
        logger.info(f"Using {self.max_workers} parallel workers")
        logger.info(f"Max retries per species: {self.max_retries}")
        logger.info(f"Rate limit delay: {self.rate_limit_delay}s")
        
        if total_pending == 0:
            logger.info("All species already downloaded!")
            return
        
        start_time = time()
        
        # Download with parallel workers
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all tasks
            future_to_tax = {
                executor.submit(self.download_single_species, tax_id): tax_id 
                for tax_id in pending_tax_ids
            }
            
            # Process completed downloads with progress bar
            try:
                with tqdm(total=total_pending, desc="Downloading STRING", unit="species") as pbar:
                    for future in as_completed(future_to_tax):
                        if shutdown_event.is_set():
                            logger.warning("Shutdown detected, cancelling remaining tasks...")
                            executor.shutdown(wait=False, cancel_futures=True)
                            break
                        
                        try:
                            tax_id, success, error, count = future.result(timeout=self.read_timeout + 60)
                            
                            if success:
                                self._mark_completed(tax_id, count)
                            else:
                                self._mark_failed(tax_id, error or "Unknown error")
                            
                            # Periodic save
                            self._save_progress()
                            
                        except FuturesTimeoutError:
                            tax_id = future_to_tax[future]
                            logger.error(f"Future timeout for tax_id {tax_id}")
                            self._mark_failed(tax_id, "Future execution timeout")
                        except Exception as e:
                            tax_id = future_to_tax[future]
                            logger.error(f"Unexpected error for tax_id {tax_id}: {e}")
                            self._mark_failed(tax_id, str(e))
                        
                        pbar.update(1)
                        stats = self.stats.get_summary()
                        pbar.set_postfix({
                            'ok': stats['success'],
                            'fail': stats['failed'],
                            'rate': f"{stats['rate_per_minute']:.1f}/min"
                        })
                        
            except KeyboardInterrupt:
                logger.warning("Keyboard interrupt received, saving progress...")
                shutdown_event.set()
        
        # Final save
        self._save_progress(force=True)
        self._save_failed()
        
        elapsed_time = time() - start_time
        
        # Summary
        stats = self.stats.get_summary()
        logger.info("=" * 80)
        logger.info("Download Summary")
        logger.info("=" * 80)
        logger.info(f"Total time: {elapsed_time / 60:.2f} minutes")
        logger.info(f"Successfully downloaded: {stats['success']}")
        logger.info(f"Failed: {stats['failed']}")
        logger.info(f"Total retries: {stats['retries']}")
        logger.info(f"Total interactions: {stats['total_interactions']:,}")
        logger.info(f"Total completed: {len(self.completed)}/{total_species}")
        logger.info(f"Progress: {len(self.completed)/total_species*100:.1f}%")
        logger.info(f"Cache location: {self.cache_dir}")
        logger.info(f"Progress file: {PROGRESS_FILE}")
        logger.info(f"Failed list: {FAILED_FILE}")
        
        if self.failed:
            logger.info(f"\nFailed species by error type:")
            error_types = {}
            for tid, info in self.failed.items():
                etype = info.get('error_type', 'unknown')
                error_types[etype] = error_types.get(etype, 0) + 1
            for etype, count in sorted(error_types.items(), key=lambda x: -x[1]):
                logger.info(f"  {etype}: {count}")
        
        logger.info("=" * 80)
        
        if shutdown_event.is_set():
            logger.info("Download was interrupted but progress has been saved. "
                       "Run again to resume.")


def main():
    parser = argparse.ArgumentParser(
        description='Download STRING database with enhanced resilience',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with defaults
  python download_string_data.py
  
  # Use fewer workers for unstable networks
  python download_string_data.py --workers 2
  
  # Increase timeouts for slow connections
  python download_string_data.py --connect-timeout 60 --read-timeout 600
  
  # Add delay between downloads to avoid rate limiting
  python download_string_data.py --rate-limit-delay 2.0
  
  # Fresh start (ignore previous progress)
  python download_string_data.py --no-resume
        """
    )
    parser.add_argument(
        '--workers', type=int, default=DEFAULT_WORKERS,
        help=f'Number of parallel download workers (default: {DEFAULT_WORKERS})'
    )
    parser.add_argument(
        '--max-retries', type=int, default=DEFAULT_MAX_RETRIES,
        help=f'Maximum retry attempts per species (default: {DEFAULT_MAX_RETRIES})'
    )
    parser.add_argument(
        '--cache-dir', type=str, default=CACHE_DIR,
        help=f'Cache directory (default: {CACHE_DIR})'
    )
    parser.add_argument(
        '--connect-timeout', type=int, default=DEFAULT_CONNECT_TIMEOUT,
        help=f'Connection timeout in seconds (default: {DEFAULT_CONNECT_TIMEOUT})'
    )
    parser.add_argument(
        '--read-timeout', type=int, default=DEFAULT_READ_TIMEOUT,
        help=f'Read timeout in seconds (default: {DEFAULT_READ_TIMEOUT})'
    )
    parser.add_argument(
        '--rate-limit-delay', type=float, default=0.5,
        help='Delay between downloads to avoid rate limiting (default: 0.5s)'
    )
    parser.add_argument(
        '--no-resume', action='store_true',
        help='Start fresh, ignoring previous progress'
    )
    parser.add_argument(
        '--retry-failed', action='store_true', default=True,
        help='Retry previously failed downloads (default: True)'
    )
    parser.add_argument(
        '--verbose', '-v', action='store_true',
        help='Enable verbose/debug logging'
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    downloader = StringDownloader(
        cache_dir=args.cache_dir,
        max_workers=args.workers,
        max_retries=args.max_retries,
        connect_timeout=args.connect_timeout,
        read_timeout=args.read_timeout,
        rate_limit_delay=args.rate_limit_delay
    )
    
    downloader.download_all(
        skip_completed=not args.no_resume,
        retry_failed=args.retry_failed
    )


if __name__ == "__main__":
    main()
