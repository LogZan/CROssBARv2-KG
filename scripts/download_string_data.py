#!/usr/bin/env python
"""
STRING Database Downloader with Resume and Retry Support

This script downloads STRING protein interaction data for all species
with robust error handling, retry mechanism, and resume capability.

Usage:
    python download_string_data.py [--workers 8] [--max-retries 3]
"""

import os
import sys
import json
import argparse
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
from time import time, sleep
from typing import Set, Dict

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

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Thread-safe locks
progress_lock = Lock()
failed_lock = Lock()


class StringDownloader:
    """STRING database downloader with resume and retry capabilities."""
    
    def __init__(self, cache_dir: str, max_workers: int = 8, max_retries: int = 3):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.max_workers = max_workers
        self.max_retries = max_retries
        
        # Setup pypath cache
        settings.setup(cachedir=str(self.cache_dir))
        
        # Load progress
        self.completed: Set[str] = self._load_progress()
        self.failed: Dict[str, str] = self._load_failed()
        
        logger.info(f"Initialized downloader with cache dir: {self.cache_dir}")
        logger.info(f"Already completed: {len(self.completed)} species")
        logger.info(f"Previously failed: {len(self.failed)} species")
    
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
    
    def _save_progress(self):
        """Save current progress."""
        with progress_lock:
            try:
                with open(PROGRESS_FILE, 'w') as f:
                    json.dump({
                        'completed': list(self.completed),
                        'total_completed': len(self.completed),
                        'timestamp': time()
                    }, f, indent=2)
            except Exception as e:
                logger.error(f"Failed to save progress: {e}")
    
    def _load_failed(self) -> Dict[str, str]:
        """Load list of failed downloads with error messages."""
        if os.path.exists(FAILED_FILE):
            try:
                with open(FAILED_FILE, 'r') as f:
                    return json.load(f)
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
    
    def _mark_completed(self, tax_id: str):
        """Mark a tax_id as completed."""
        with progress_lock:
            self.completed.add(tax_id)
            # Remove from failed if it was there
            if tax_id in self.failed:
                with failed_lock:
                    del self.failed[tax_id]
    
    def _mark_failed(self, tax_id: str, error: str):
        """Mark a tax_id as failed."""
        with failed_lock:
            self.failed[tax_id] = error
    
    def download_single_species(self, tax_id: str, retry_count: int = 0) -> tuple:
        """
        Download STRING data for a single species.
        
        Args:
            tax_id: NCBI taxonomy ID
            retry_count: Current retry attempt number
            
        Returns:
            Tuple of (tax_id, success, error_message)
        """
        try:
            # Attempt download
            logger.debug(f"Downloading tax_id {tax_id} (attempt {retry_count + 1})")
            
            # Call STRING API - this will use pypath's caching automatically
            interactions = list(string.string_links_interactions(
                ncbi_tax_id=int(tax_id),
                score_threshold="high_confidence"
            ))
            
            logger.info(f"✓ Downloaded {len(interactions)} interactions for tax_id {tax_id}")
            return (tax_id, True, None)
            
        except Exception as e:
            error_msg = str(e)
            
            # Check if we should retry
            if retry_count < self.max_retries:
                # Exponential backoff
                wait_time = 2 ** retry_count
                logger.warning(f"✗ Failed tax_id {tax_id}: {error_msg}. Retrying in {wait_time}s...")
                sleep(wait_time)
                return self.download_single_species(tax_id, retry_count + 1)
            else:
                logger.error(f"✗ Failed tax_id {tax_id} after {self.max_retries} retries: {error_msg}")
                return (tax_id, False, error_msg)
    
    def download_all(self, skip_completed: bool = True, retry_failed: bool = True):
        """
        Download STRING data for all species.
        
        Args:
            skip_completed: Skip already downloaded species
            retry_failed: Retry previously failed downloads
        """
        logger.info("=" * 80)
        logger.info("Starting STRING data download")
        logger.info("=" * 80)
        
        # Get list of all species
        logger.info("Fetching list of all species...")
        string_species = string.string_species()
        all_tax_ids = list(string_species.keys())
        
        # Filter out problematic tax IDs
        tax_ids_to_skip = {'4565', '8032'}  # Known problematic IDs
        
        # Determine which tax_ids to download
        if skip_completed:
            pending_tax_ids = [
                tid for tid in all_tax_ids 
                if tid not in self.completed and tid not in tax_ids_to_skip
            ]
        else:
            pending_tax_ids = [tid for tid in all_tax_ids if tid not in tax_ids_to_skip]
        
        # Add failed ones if retry_failed is True
        if retry_failed and self.failed:
            failed_to_retry = [tid for tid in self.failed.keys() if tid not in pending_tax_ids]
            pending_tax_ids.extend(failed_to_retry)
            logger.info(f"Will retry {len(failed_to_retry)} previously failed downloads")
        
        total_species = len(all_tax_ids)
        total_pending = len(pending_tax_ids)
        
        logger.info(f"Total species: {total_species}")
        logger.info(f"Already completed: {len(self.completed)}")
        logger.info(f"Pending downloads: {total_pending}")
        logger.info(f"Using {self.max_workers} parallel workers")
        logger.info(f"Max retries per species: {self.max_retries}")
        
        if total_pending == 0:
            logger.info("All species already downloaded!")
            return
        
        start_time = time()
        success_count = 0
        failed_count = 0
        
        # Download with parallel workers
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all tasks
            future_to_tax = {
                executor.submit(self.download_single_species, tax_id): tax_id 
                for tax_id in pending_tax_ids
            }
            
            # Process completed downloads with progress bar
            with tqdm(total=total_pending, desc="Downloading STRING", unit="species") as pbar:
                for future in as_completed(future_to_tax):
                    tax_id, success, error = future.result()
                    
                    if success:
                        self._mark_completed(tax_id)
                        success_count += 1
                        # Save progress every 10 successful downloads
                        if success_count % 10 == 0:
                            self._save_progress()
                    else:
                        self._mark_failed(tax_id, error)
                        failed_count += 1
                        # Save failed list
                        if failed_count % 5 == 0:
                            self._save_failed()
                    
                    pbar.update(1)
                    pbar.set_postfix({
                        'success': success_count,
                        'failed': failed_count,
                        'rate': f"{success_count / (time() - start_time) * 60:.1f} species/min"
                    })
        
        # Final save
        self._save_progress()
        self._save_failed()
        
        elapsed_time = time() - start_time
        
        # Summary
        logger.info("=" * 80)
        logger.info("Download Summary")
        logger.info("=" * 80)
        logger.info(f"Total time: {elapsed_time / 60:.2f} minutes")
        logger.info(f"Successfully downloaded: {success_count}")
        logger.info(f"Failed: {failed_count}")
        logger.info(f"Total completed: {len(self.completed)}/{total_species}")
        logger.info(f"Progress: {len(self.completed)/total_species*100:.1f}%")
        logger.info(f"Cache location: {self.cache_dir}")
        logger.info(f"Progress file: {PROGRESS_FILE}")
        logger.info(f"Failed list: {FAILED_FILE}")
        logger.info("=" * 80)


def main():
    parser = argparse.ArgumentParser(description='Download STRING database with resume support')
    parser.add_argument('--workers', type=int, default=8,
                       help='Number of parallel download workers (default: 8)')
    parser.add_argument('--max-retries', type=int, default=3,
                       help='Maximum retry attempts per species (default: 3)')
    parser.add_argument('--cache-dir', type=str, default=CACHE_DIR,
                       help=f'Cache directory (default: {CACHE_DIR})')
    parser.add_argument('--no-resume', action='store_true',
                       help='Start fresh, ignoring previous progress')
    parser.add_argument('--retry-failed', action='store_true', default=True,
                       help='Retry previously failed downloads (default: True)')
    
    args = parser.parse_args()
    
    downloader = StringDownloader(
        cache_dir=args.cache_dir,
        max_workers=args.workers,
        max_retries=args.max_retries
    )
    
    downloader.download_all(
        skip_completed=not args.no_resume,
        retry_failed=args.retry_failed
    )


if __name__ == "__main__":
    main()
