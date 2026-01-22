#!/usr/bin/env python
"""
Test splitting STRING into per-species cache files and reading via PPI adapter cache logic.

- Extracts 10 species from the large STRING file into a mini gzip input.
- Runs split_string_safe.split_string_safe to create cache files (overwrite if exist).
- Uses PPI.download_string_data(cache=True) with patched species list and uniprot mapping.
"""

import argparse
import collections
import gzip
import os
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Set, Tuple

# Ensure repo root and test dir are importable
CURRENT_DIR = Path(__file__).resolve().parent
REPO_ROOT = CURRENT_DIR.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(CURRENT_DIR))

from split_string_safe import (
    DEFAULT_INPUT,
    DEFAULT_OUTPUT,
    generate_pypath_hash,
    is_complete_file,
    split_string_safe,
)

from bccb.ppi_adapter import PPI
from bccb import cache_config

from pypath.inputs import string as pypath_string
from pypath.inputs import uniprot as pypath_uniprot


def _collect_species_lines(
    input_path: str,
    species_count: int,
    min_lines_per_species: int,
    max_lines_total: int,
    min_tax_id: int,
    only_tax_ids: Set[str],
) -> Tuple[List[str], Dict[str, List[str]], Set[str]]:
    """Collect header and lines for the first N species with at least min lines each."""
    species_lines: Dict[str, List[str]] = collections.defaultdict(list)
    species_set: Set[str] = set()
    total_lines = 0

    with gzip.open(input_path, "rt", encoding="utf-8") as fh:
        header = fh.readline().strip()
        for line in fh:
            if max_lines_total > 0 and total_lines >= max_lines_total:
                break
            parts = line.split(None, 1)
            if not parts:
                continue
            tax_id = parts[0].split(".", 1)[0]
            if only_tax_ids and tax_id not in only_tax_ids:
                continue
            try:
                tax_id_int = int(tax_id)
            except ValueError:
                continue
            if tax_id_int < min_tax_id:
                continue

            if tax_id not in species_set and len(species_set) >= species_count:
                # Only keep collecting if species already selected
                if tax_id not in species_set:
                    continue

            species_set.add(tax_id)
            species_lines[tax_id].append(line)
            total_lines += 1

            # Stop once all species have enough lines
            if len(species_set) >= species_count:
                ready = all(
                    len(species_lines[t]) >= min_lines_per_species
                    for t in species_set
                )
                if ready:
                    break

    if len(species_set) < species_count:
        raise RuntimeError(
            f"Only found {len(species_set)} species; expected {species_count}."
        )

    return header, species_lines, species_set


def _write_mini_input(tmp_dir: Path, header: str, species_lines: Dict[str, List[str]]) -> str:
    """Write a mini gzip input file for splitting."""
    tmp_path = tmp_dir / "string_subset.txt.gz"
    with gzip.open(tmp_path, "wt", encoding="utf-8") as fh:
        fh.write(header + "\n")
        for tax_id, lines in species_lines.items():
            for line in lines:
                fh.write(line)
    return str(tmp_path)


def _remove_existing_cache_files(output_dir: Path, tax_ids: Set[str]) -> None:
    """Remove existing cache files for specified tax IDs to allow overwrite."""
    for tax_id in tax_ids:
        file_hash = generate_pypath_hash(tax_id)
        pattern = f"{file_hash}-{tax_id}.protein.links.detailed.*.txt.gz"
        for f in output_dir.glob(pattern):
            try:
                f.unlink()
            except Exception:
                pass


def _build_fake_uniprot_mapping(species_lines: Dict[str, List[str]]) -> Dict[str, str]:
    """Build a minimal uniprot->string mapping from sampled STRING lines."""
    mapping: Dict[str, str] = {}
    for tax_id, lines in species_lines.items():
        for line in lines:
            parts = line.split()
            if len(parts) < 2:
                continue
            prot_a = parts[0]
            prot_b = parts[1]
            for prot in (prot_a, prot_b):
                # Use a stable fake uniprot ID per protein
                up_id = f"UP_{tax_id}_{prot.replace('.', '_')}"
                mapping[up_id] = prot
    return mapping


def _patch_pypath_sources(tax_ids: Set[str], uniprot_mapping: Dict[str, str]):
    """Patch pypath string species + uniprot mapping to avoid network access."""
    original_string_species = pypath_string.string_species
    original_uniprot_data = pypath_uniprot.uniprot_data
    original_all_uniprots = pypath_uniprot._all_uniprots

    def _fake_string_species():
        return {tax_id: f"species_{tax_id}" for tax_id in tax_ids}

    def _fake_uniprot_data(*args, **kwargs):
        return dict(uniprot_mapping)

    pypath_string.string_species = _fake_string_species
    pypath_uniprot.uniprot_data = _fake_uniprot_data
    pypath_uniprot._all_uniprots = lambda *args, **kwargs: []

    return original_string_species, original_uniprot_data, original_all_uniprots


def main():
    parser = argparse.ArgumentParser(
        description="Test STRING split cache files are readable by PPI adapter."
    )
    parser.add_argument(
        "--input",
        type=str,
        default=DEFAULT_INPUT,
        help=f"Input STRING file (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=DEFAULT_OUTPUT,
        help=f"Output cache dir (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--species",
        type=int,
        default=10,
        help="Number of species to include (default: 10)",
    )
    parser.add_argument(
        "--min-lines",
        type=int,
        default=5,
        help="Min lines per species in subset (default: 5)",
    )
    parser.add_argument(
        "--max-lines",
        type=int,
        default=50000000,
        help="Max total lines to scan, 0 means no limit (default: 50000000)",
    )
    parser.add_argument(
        "--min-tax-id",
        type=int,
        default=1001,
        help="Only include species with tax_id >= this value (default: 1001)",
    )
    parser.add_argument(
        "--tax-ids",
        type=str,
        default="",
        help="Comma-separated tax_ids to include (overrides min-tax-id/species count)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Workers for split (default: 1)",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=10000,
        help="Split chunk size (default: 10000)",
    )

    args = parser.parse_args()

    if not os.path.exists(args.input):
        raise SystemExit(f"Input file not found: {args.input}")

    # Ensure cache dirs set up
    cache_config.setup_pypath_cache()
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Collect species and build mini input
    only_tax_ids: Set[str] = set()
    if args.tax_ids:
        only_tax_ids = {t.strip() for t in args.tax_ids.split(",") if t.strip()}
        args.species = len(only_tax_ids)

    header, species_lines, tax_ids = _collect_species_lines(
        args.input,
        args.species,
        args.min_lines,
        args.max_lines,
        args.min_tax_id,
        only_tax_ids,
    )

    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        mini_input = _write_mini_input(tmp_dir, header, species_lines)

        # Overwrite existing cache files for these species
        _remove_existing_cache_files(output_dir, tax_ids)

        # Run split
        split_string_safe(
            input_file=mini_input,
            output_dir=str(output_dir),
            workers=args.workers,
            skip_existing=False,
            max_handles=256,
            start_line=0,
            chunk_size=args.chunk_size,
        )

    # Verify split files exist and are valid
    missing = []
    invalid = []
    for tax_id in tax_ids:
        file_hash = generate_pypath_hash(tax_id)
        filename = f"{file_hash}-{tax_id}.protein.links.detailed.v12.0.txt.gz"
        path = output_dir / filename
        if not path.exists():
            missing.append(str(path))
        elif not is_complete_file(path):
            invalid.append(str(path))

    if missing or invalid:
        raise SystemExit(
            f"Split output check failed. Missing: {len(missing)} Invalid: {len(invalid)}"
        )

    # Build fake uniprot mapping for filtering
    uniprot_mapping = _build_fake_uniprot_mapping(species_lines)

    # Patch pypath to avoid full species list / network
    original_string_species, original_uniprot_data, original_all_uniprots = _patch_pypath_sources(
        tax_ids, uniprot_mapping
    )

    try:
        ppi = PPI(organism=None, test_mode=False)
        ppi.download_string_data(cache=True)
        if not ppi.string_ints:
            raise SystemExit("PPI adapter returned 0 STRING interactions.")
        print(
            f"OK: Read {len(ppi.string_ints)} STRING interactions from cache for {len(tax_ids)} species."
        )
    finally:
        pypath_string.string_species = original_string_species
        pypath_uniprot.uniprot_data = original_uniprot_data
        pypath_uniprot._all_uniprots = original_all_uniprots


if __name__ == "__main__":
    main()
