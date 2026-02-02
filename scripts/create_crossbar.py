import os
import sys
import gc
import csv
import json
import re
import yaml
import ctypes
import argparse
import logging
from collections import defaultdict
from datetime import datetime, timezone, timedelta
from pathlib import Path

# Setup timezone (UTC+8 for Asia/Shanghai)
TZ = timezone(timedelta(hours=8))

project_root = Path(__file__).resolve().parent.parent
# Add parent directory to Python path
sys.path.insert(0, str(project_root))

# CRITICAL: Setup pypath cache directory FIRST before any pypath imports
# This ensures all temp files go to the large filesystem, not /tmp or /root
from bccb import cache_config
cache_config.setup_pypath_cache()

# Load configuration first
with open(project_root / "config/crossbar_config.yaml", 'r') as f:
    crossbar_config = yaml.safe_load(f)

with open(project_root / "config/biocypher_config.yaml", 'r') as f:
    biocypher_config = yaml.safe_load(f)

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Build CROssBARv2 Knowledge Graph')
parser.add_argument('--cache', type=lambda x: x.lower() == 'true', default=crossbar_config['settings']['cache'], help='Use cached data')
parser.add_argument('--test-mode', type=lambda x: x.lower() == 'true', default=crossbar_config['settings']['test_mode'], help='Enable test mode')
parser.add_argument('--no-checkpoint', action='store_true', help='Ignore checkpoint and start from beginning')
parser.add_argument('--skip-until', type=str, help='Skip adapters until reaching specified adapter')
parser.add_argument('--only', type=str, help='Only run specified adapters (comma-separated)')
parser.add_argument('--reset-checkpoint', action='store_true', help='Clear checkpoint and restart')
parser.add_argument('--stats-version', type=str, help='Stats output version (defaults to output_dir basename)')
parser.add_argument('--stats-only', action='store_true', help='Only compute stats from existing output, do not run adapters')
args = parser.parse_args()

# Calculate output_dir_path from config
output_dir_path = str(project_root / crossbar_config['settings'].get('output_dir', 'biocypher-out'))
csv_output_dir = Path(output_dir_path).parent / "csv"
stats_version = args.stats_version or Path(output_dir_path).name
stats_dir = project_root / "stats" / stats_version
stats_dir.mkdir(parents=True, exist_ok=True)

# Initialize checkpoint manager
checkpoint_enabled = crossbar_config.get('checkpoint', {}).get('enabled', True) and not args.no_checkpoint
checkpoint_file = str(Path(output_dir_path) / '.checkpoint.json')

from bccb.checkpoint import CheckpointManager
checkpoint = CheckpointManager(checkpoint_file) if checkpoint_enabled else None

# Handle checkpoint reset
if args.reset_checkpoint and checkpoint:
    checkpoint.reset()
    print("Checkpoint reset. Starting fresh.")

# Parse 'only' argument
only_adapters = set(args.only.split(',')) if args.only else None

# Use checkpoint flag
use_checkpoint = checkpoint_enabled


def aggressive_memory_cleanup(adapter_name: str = ""):
    """Perform aggressive memory cleanup between adapters."""
    # Force garbage collection multiple times
    gc.collect()
    gc.collect()
    gc.collect()
    
    # Try to release memory back to OS (Linux specific)
    try:
        libc = ctypes.CDLL("libc.so.6")
        libc.malloc_trim(0)
    except Exception:
        pass
    
    # Print memory status
    try:
        import resource
        mem_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024  # MB
        print(f"[{adapter_name}] Memory cleanup done. Peak RSS: {mem_usage:.0f} MB")
    except Exception:
        print(f"[{adapter_name}] Memory cleanup done.")


# Patch pypath settings for compatibility
try:
    from pypath.share import settings
    if not hasattr(settings, 'context') and hasattr(settings, 'settings'):
        settings.context = settings.settings.context
except ImportError:
    pass

from bccb.uniprot_swissprot_adapter import (
    UniprotSwissprot,
    UniprotNodeType,
    UniprotNodeField,
    UniprotEdgeType,
    UniprotIDField,
)

from bccb.uniprot_keywords_adapter import (
    UniprotKeywords,
    KeywordNodeType,
    KeywordEdgeType,
)

from bccb.ppi_adapter import (
    PPI
)

from bccb.interpro_adapter import (
    InterPro
)

from bccb.go_adapter import (
    GO
)

from bccb.drug_adapter import (
    Drug,
    DrugNodeField,
)

from bccb.compound_adapter import (
    Compound
)

from bccb.orthology_adapter import (
    Orthology
)

from bccb.disease_adapter import (
    Disease
)

from bccb.phenotype_adapter import (
    HPO
)

from bccb.pathway_adapter import (
    Pathway
)

from bccb.side_effect_adapter import (
    SideEffect
)

from bccb.ec_adapter import (
    EC
)

from bccb.tfgen_adapter import (
    TFGene
)

from biocypher import BioCypher

# Load configuration
with open(project_root / "config/crossbar_config.yaml", 'r') as f:
    config = yaml.safe_load(f)

# Extract config values
timestamp = datetime.now(TZ).strftime("%Y%m%d%H%M%S")
embeddings_dir = config['data_paths']['embeddings_dir']
malacards_dir_path = config['data_paths']['malacards_dir']
uniprot_json_path = config['data_paths']['uniprot_json']
interpro_dir_path = config['data_paths']['interpro_dir']

# Embedding file paths
emb_cfg = config.get("embeddings", {}) or {}
prott5_name = emb_cfg.get("prott5")
esm2_name = emb_cfg.get("esm2")
nt_name = emb_cfg.get("nucleotide_transformer")

prott5_embedding_path = f"{embeddings_dir}/{prott5_name}" if prott5_name else None
esm2_embedding_path = f"{embeddings_dir}/{esm2_name}" if esm2_name else None
nt_embedding_path = f"{embeddings_dir}/{nt_name}" if nt_name else None
selformer_drug_embedding_path = f"{embeddings_dir}/{config['embeddings']['selformer_drug']}"
selformer_compound_embedding_path = f"{embeddings_dir}/{config['embeddings']['selformer_compound']}"
doc2vec_disease_embedding_path = f"{embeddings_dir}/{config['embeddings']['doc2vec_disease']}"
biokeen_pathway_embedding_path = f"{embeddings_dir}/{config['embeddings']['biokeen_pathway']}"
rxnfp_ec_embedding_path = f"{embeddings_dir}/{config['embeddings']['rxnfp_ec']}"
anc2vec_go_embedding_path = f"{embeddings_dir}/{config['embeddings']['anc2vec_go']}"
cada_phenotype_embedding_path = f"{embeddings_dir}/{config['embeddings']['cada_phenotype']}"
dom2vec_domain_embedding_path = f"{embeddings_dir}/{config['embeddings']['dom2vec_domain']}"

# MalaCards file paths
malacards_json_path = f"{malacards_dir_path}/{config['malacards']['diseases']}"
malacards_related_diseases_json_path = f"{malacards_dir_path}/{config['malacards']['related_diseases']}"


# Organism parameter conversion helpers for pypath compatibility
# TODO: Unused, now setting organism=None in majority of uniprot.uniprot_data() calls
def organism_for_uniprot_data(organism):
    """Convert organism value for uniprot.uniprot_data() calls.

    In pypath v0.16+, uniprot_data() uses None for all organisms, not '*'.
    """
    return None if organism == "*" else organism

def organism_for_all_uniprots(organism):
    """Convert organism value for uniprot._all_uniprots() calls.

    _all_uniprots() uses '*' for all organisms.
    """
    return organism


# Helper for chunked edge writing to avoid OOM with large generators
def write_edges_chunked(bc, edge_generator, chunk_size=1000000):
    """Write edges in chunks to avoid BioCypher's list(edges) OOM issue.
    
    BioCypher's _batch_writer.py does `edges = list(edges)` which loads
    entire generator into memory. For large datasets (e.g., 77GB InterPro),
    this causes OOM. 
    
    WORKAROUND: We patch BioCypher's batch_writer to not call list() on 
    the entire generator. Instead, we pass a streaming generator directly.
    BioCypher will then batch by its internal batch_size (chunk_size param).
    """
    from biocypher.output.write._batch_writer import _BatchWriter
    from biocypher._create import BioCypherRelAsNode
    from biocypher._logger import logger
    from more_itertools import peekable
    
    # Monkey-patch write_edges to skip list() conversion
    original_write_edges = _BatchWriter.write_edges
    
    def patched_write_edges(self, edges, batch_size=int(1e6)):
        """Patched version that doesn't convert generator to list."""
        edges = peekable(edges)
        try:
            edges.peek()
        except StopIteration:
            logger.debug("No edges to write.")
            return True
        
        nodes_flat = []
        
        def edge_stream():
            for edge in edges:
                if isinstance(edge, BioCypherRelAsNode):
                    if self.deduplicator.rel_as_node_seen(edge):
                        continue
                    nodes_flat.append(edge.get_node())
                    yield edge.get_source_edge()
                    yield edge.get_target_edge()
                else:
                    if self.deduplicator.edge_seen(edge):
                        continue
                    yield edge
        
        passed = self._write_edge_data(edge_stream(), batch_size)
        
        if nodes_flat:
            self.write_nodes(nodes_flat)
        
        if not passed:
            logger.error("Error while writing edge data.")
            return False
            
        passed = self._write_edge_headers()
        if not passed:
            logger.error("Error while writing edge headers.")
            return False
        return True
    
    # Apply patch and write edges
    _BatchWriter.write_edges = patched_write_edges
    try:
        bc.write_edges(edge_generator, batch_size=chunk_size)
    finally:
        _BatchWriter.write_edges = original_write_edges


# Helper for consistent logging
def log_adapter_boundary(adapter_name: str, phase: str):
    """Log adapter execution boundary."""
    separator = "=" * 60
    timestamp = datetime.now(TZ).strftime("%Y-%m-%d %H:%M:%S")
    if phase == "start":
        print(f"\n{separator}")
        print(f"[{timestamp}] [{timestamp}] [START] {adapter_name} adapter")
        print(f"[{timestamp}] {separator}")
    else:
        print(f"[{timestamp}] [END] {adapter_name} adapter")
        print(f"{separator}\n")


STATS_LOGGER = logging.getLogger("crossbar_stats")
if not STATS_LOGGER.handlers:
    logging.basicConfig(
        level=logging.INFO,
        format="[STATS] {message}",
        style="{",
    )
STATS_LOGGER.disabled = True

ADAPTER_ORDER = [
    "uniprot",
    "keywords",
    "ppi",
    "interpro",
    "go",
    "drug",
    "compound",
    "orthology",
    "disease",
    "phenotype",
    "pathway",
    "side_effect",
    "ec",
    "tfgene",
]

ADAPTER_DISPLAY_NAMES = {
    "uniprot": "UniProtKB-SwissProt",
    "keywords": "UniProtKB-Keywords",
    "ppi": "PPI",
    "interpro": "InterPro",
    "go": "Gene Ontology",
    "drug": "Drug",
    "compound": "Compound",
    "orthology": "Orthology",
    "disease": "Disease",
    "phenotype": "Phenotype",
    "pathway": "Pathway",
    "side_effect": "Side Effect",
    "ec": "EC",
    "tfgene": "TFgene",
}

HEADER_SUFFIX = "-header.csv"
PART_RE = re.compile(r"^(?P<base>.+)-part\d+\.csv$")


def _snapshot_csv_files(dir_path: str) -> dict:
    snapshot = {}
    for path in Path(dir_path).glob("*.csv"):
        try:
            stat = path.stat()
        except FileNotFoundError:
            continue
        snapshot[path.name] = (stat.st_size, stat.st_mtime)
    return snapshot


def _diff_csv_files(before: dict, after: dict) -> list:
    changed = []
    for name, meta in after.items():
        if name not in before or before[name] != meta:
            changed.append(name)
    return sorted(changed)


def _derive_type_name(filename: str) -> str:
    if filename.endswith(HEADER_SUFFIX):
        return filename[:-len(HEADER_SUFFIX)]
    match = PART_RE.match(filename)
    if match:
        return match.group("base")
    if filename.endswith(".csv"):
        return filename[:-4]
    return filename


def _is_header_file(filename: str) -> bool:
    return filename.endswith(HEADER_SUFFIX)


def _is_part_file(filename: str) -> bool:
    return PART_RE.match(filename) is not None


def _is_single_data_file(filename: str) -> bool:
    return filename.endswith(".csv") and not _is_header_file(filename) and not _is_part_file(filename)


def _read_header_columns(path: Path) -> list:
    with path.open("r", encoding="utf-8") as f:
        line = f.readline().rstrip("\n")
    return line.split("\t") if line else []


def _count_lines(path: Path) -> int:
    count = 0
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            count += chunk.count(b"\n")
    return count


def _infer_type_kind(columns: list) -> str:
    edge_markers = (":START_ID", ":END_ID", ":TYPE")
    for col in columns:
        if col.startswith(edge_markers):
            return "edge"
    return "node"


def _parse_sources(raw: str) -> list:
    raw = raw.strip()
    if not raw:
        return []
    raw = raw.strip("\"' ")
    if ";" in raw:
        return [part.strip().strip("\"' ") for part in raw.split(";") if part.strip()]
    return [raw]


def _normalize_ppi_source(source: str) -> str:
    normalized = source.strip().strip("\"' ")
    upper = normalized.upper()
    if upper == "INTACT":
        return "PPI-IntAct"
    if upper in ("BIOGRID", "BIOGRID "):
        return "PPI-BioGrid"
    if upper == "STRING":
        return "PPI-STRING"
    return ""


class StatsManager:
    def __init__(self, output_dir: str, stats_dir_path: Path, test_mode: bool):
        self.output_dir = Path(output_dir)
        self.stats_dir = stats_dir_path
        self.mapping_path = self.stats_dir / "adapter_files.json"
        self.adapter_files = {}
        self.mapping_mode = False
        self._snapshot_before = None

        if self.mapping_path.exists():
            self._load_mapping()
        elif test_mode:
            self.mapping_mode = True
        else:
            STATS_LOGGER.warning(
                "adapter_files.json not found. Run once in test mode to generate mapping."
            )

    def _load_mapping(self):
        with self.mapping_path.open("r", encoding="utf-8") as f:
            self.adapter_files = json.load(f)

    def _save_mapping(self):
        with self.mapping_path.open("w", encoding="utf-8") as f:
            json.dump(self.adapter_files, f, indent=2, sort_keys=True)

    def on_adapter_start(self, adapter_name: str):
        if self.mapping_mode:
            self._snapshot_before = _snapshot_csv_files(str(self.output_dir))

    def on_adapter_end(self, adapter_name: str):
        if self.mapping_mode:
            snapshot_after = _snapshot_csv_files(str(self.output_dir))
            changed = _diff_csv_files(self._snapshot_before or {}, snapshot_after)
            if changed:
                self.adapter_files[adapter_name] = changed
                self._save_mapping()
        self.write_stats(current_adapter=adapter_name)

    def _collect_type_info(self, adapter_name: str) -> dict:
        files = self.adapter_files.get(adapter_name, [])
        type_info = {}
        for filename in files:
            if not filename.endswith(".csv"):
                continue
            type_name = _derive_type_name(filename)
            info = type_info.setdefault(
                type_name,
                {"all_files": [], "header_files": [], "data_files": []},
            )
            info["all_files"].append(filename)
            if _is_header_file(filename):
                info["header_files"].append(filename)
            elif _is_part_file(filename) or _is_single_data_file(filename):
                info["data_files"].append(filename)
        return type_info

    def _get_header_columns(self, type_name: str, info: dict) -> list:
        if info["header_files"]:
            header_path = self.output_dir / info["header_files"][0]
            if header_path.exists():
                return _read_header_columns(header_path)
        for data_file in info["data_files"]:
            if _is_single_data_file(data_file):
                data_path = self.output_dir / data_file
                if data_path.exists():
                    return _read_header_columns(data_path)
        STATS_LOGGER.warning(f"No header found for type {type_name}")
        return []

    def _count_rows_and_size(self, info: dict) -> tuple:
        row_count = 0
        size_bytes = 0
        for filename in info["data_files"]:
            path = self.output_dir / filename
            if not path.exists():
                continue
            size_bytes += path.stat().st_size
            lines = _count_lines(path)
            if _is_single_data_file(filename):
                row_count += max(0, lines - 1)
            else:
                row_count += lines
        for filename in info["header_files"]:
            path = self.output_dir / filename
            if path.exists():
                size_bytes += path.stat().st_size
        return row_count, size_bytes

    def _compute_schema_rows(self, adapter_name: str) -> list:
        rows = []
        type_info = self._collect_type_info(adapter_name)
        for type_name, info in sorted(type_info.items()):
            columns = self._get_header_columns(type_name, info)
            type_kind = _infer_type_kind(columns)
            row_count, size_bytes = self._count_rows_and_size(info)
            rows.append(
                {
                    "adapter": ADAPTER_DISPLAY_NAMES.get(adapter_name, adapter_name),
                    "type_name": type_name,
                    "type_kind": type_kind,
                    "count": row_count,
                    "property_count": len(columns),
                    "size_bytes": size_bytes,
                }
            )
        return rows

    def _compute_adapter_stats(self, adapter_name: str) -> dict:
        type_info = self._collect_type_info(adapter_name)
        nodes_count = 0
        edges_count = 0
        nodes_size = 0
        edges_size = 0
        for type_name, info in type_info.items():
            columns = self._get_header_columns(type_name, info)
            type_kind = _infer_type_kind(columns)
            row_count, size_bytes = self._count_rows_and_size(info)
            if type_kind == "edge":
                edges_count += row_count
                edges_size += size_bytes
            else:
                nodes_count += row_count
                nodes_size += size_bytes
        return {
            "adapter": ADAPTER_DISPLAY_NAMES.get(adapter_name, adapter_name),
            "nodes_count": nodes_count,
            "edges_count": edges_count,
            "nodes_size_bytes": nodes_size,
            "edges_size_bytes": edges_size,
            "total_size_bytes": nodes_size + edges_size,
        }

    def _compute_ppi_sources(self) -> list:
        files = self.adapter_files.get("ppi", [])
        if not files:
            return []
        type_info = self._collect_type_info("ppi")
        ppi_info = type_info.get("Protein_interacts_with_protein")
        if not ppi_info:
            return []
        columns = self._get_header_columns("Protein_interacts_with_protein", ppi_info)
        source_idx = None
        for idx, col in enumerate(columns):
            if col == "source" or col.startswith("source"):
                source_idx = idx
                break
        if source_idx is None:
            STATS_LOGGER.warning("PPI source column not found.")
            return []

        counts = defaultdict(int)
        sizes = defaultdict(int)
        for filename in ppi_info["data_files"]:
            path = self.output_dir / filename
            with path.open("rb") as f:
                for line in f:
                    if not line:
                        continue
                    parts = line.rstrip(b"\n").split(b"\t")
                    if source_idx >= len(parts):
                        continue
                    raw = parts[source_idx].decode("utf-8", "ignore")
                    sources = [_normalize_ppi_source(s) for s in _parse_sources(raw)]
                    sources = [s for s in sources if s]
                    if not sources:
                        continue
                    size = len(line)
                    size_base = size // len(sources)
                    remainder = size % len(sources)
                    for i, src in enumerate(sources):
                        counts[src] += 1
                        sizes[src] += size_base + (1 if i < remainder else 0)

        rows = []
        for src in ("PPI-IntAct", "PPI-BioGrid", "PPI-STRING"):
            if src in counts:
                rows.append(
                    {
                        "adapter": src,
                        "nodes_count": 0,
                        "edges_count": counts[src],
                        "nodes_size_bytes": 0,
                        "edges_size_bytes": sizes[src],
                        "total_size_bytes": sizes[src],
                    }
                )
        return rows

    def _write_adapter_stats(self, adapter_rows: list):
        path = self.stats_dir / "adapter_stats.csv"
        with path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "adapter",
                    "nodes_count",
                    "edges_count",
                    "nodes_size_bytes",
                    "edges_size_bytes",
                    "total_size_bytes",
                ],
            )
            writer.writeheader()
            writer.writerows(adapter_rows)

    def _write_schema_stats(self, schema_rows: list):
        path = self.stats_dir / "schema_stats.csv"
        with path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "adapter",
                    "type_name",
                    "type_kind",
                    "count",
                    "property_count",
                    "size_bytes",
                ],
            )
            writer.writeheader()
            writer.writerows(schema_rows)

    def write_stats(self, current_adapter: str = ""):
        if not self.adapter_files:
            return

        adapter_rows = []
        schema_rows = []

        for adapter_name in ADAPTER_ORDER:
            if adapter_name not in self.adapter_files:
                continue
            if adapter_name == "ppi":
                adapter_rows.extend(self._compute_ppi_sources())
            else:
                adapter_rows.append(self._compute_adapter_stats(adapter_name))
            schema_rows.extend(self._compute_schema_rows(adapter_name))

        # Append global schema summary
        global_schema = {}
        for row in schema_rows:
            key = row["type_name"]
            entry = global_schema.setdefault(
                key,
                {
                    "adapter": "ALL",
                    "type_name": key,
                    "type_kind": row["type_kind"],
                    "count": 0,
                    "property_count": row["property_count"],
                    "size_bytes": 0,
                },
            )
            entry["count"] += row["count"]
            entry["size_bytes"] += row["size_bytes"]
            entry["property_count"] = max(entry["property_count"], row["property_count"])

        schema_rows.extend(sorted(global_schema.values(), key=lambda r: r["type_name"]))

        # Append ALL row for adapter stats
        if adapter_rows:
            total_nodes = sum(r["nodes_count"] for r in adapter_rows)
            total_edges = sum(r["edges_count"] for r in adapter_rows)
            total_nodes_size = sum(r["nodes_size_bytes"] for r in adapter_rows)
            total_edges_size = sum(r["edges_size_bytes"] for r in adapter_rows)
            adapter_rows.append(
                {
                    "adapter": "ALL",
                    "nodes_count": total_nodes,
                    "edges_count": total_edges,
                    "nodes_size_bytes": total_nodes_size,
                    "edges_size_bytes": total_edges_size,
                    "total_size_bytes": total_nodes_size + total_edges_size,
                }
            )

        self._write_adapter_stats(adapter_rows)
        self._write_schema_stats(schema_rows)

        if current_adapter:
            if current_adapter == "ppi":
                for row in adapter_rows:
                    if row["adapter"].startswith("PPI-"):
                        STATS_LOGGER.info(
                            f"{row['adapter']} edges={row['edges_count']} "
                            f"edges_size={row['edges_size_bytes']} total_size={row['total_size_bytes']}"
                        )
            else:
                for row in adapter_rows:
                    if row["adapter"] == ADAPTER_DISPLAY_NAMES.get(current_adapter, current_adapter):
                        STATS_LOGGER.info(
                            f"{row['adapter']} nodes={row['nodes_count']} edges={row['edges_count']} "
                            f"nodes_size={row['nodes_size_bytes']} edges_size={row['edges_size_bytes']} "
                            f"total_size={row['total_size_bytes']}"
                        )

        total_nodes = 0
        total_edges = 0
        total_size = 0
        for row in global_schema.values():
            if row["type_kind"] == "edge":
                total_edges += row["count"]
            else:
                total_nodes += row["count"]
            total_size += row["size_bytes"]
        STATS_LOGGER.info(
            f"ALL nodes={total_nodes} edges={total_edges} total_size={total_size}"
        )


def run_adapter(adapter_name, adapter_func):
    """Run adapter with checkpoint support."""
    # Check if should run based on checkpoint and CLI args
    if use_checkpoint and not checkpoint.should_run(adapter_name, args.skip_until, only_adapters):
        print(f"⏭️  Skipping {adapter_name} (already completed)")
        stats_manager.on_adapter_end(adapter_name)
        return True

    log_adapter_boundary(adapter_name, "start")
    stats_manager.on_adapter_start(adapter_name)
    try:
        adapter_func()
        if use_checkpoint:
            checkpoint.mark_completed(adapter_name)
        print(f"✓ {adapter_name} adapter completed successfully")
        stats_manager.on_adapter_end(adapter_name)
        return True
    except Exception as e:
        print(f"WARNING: {adapter_name} adapter failed: {e}")
        if use_checkpoint:
            checkpoint.mark_failed(adapter_name, e)
        import traceback
        traceback.print_exc()
        stats_manager.on_adapter_end(adapter_name)
        return False
    finally:
        log_adapter_boundary(adapter_name, "end")
        aggressive_memory_cleanup(adapter_name)

bc = BioCypher(
    biocypher_config_path=str(project_root / "config/biocypher_config.yaml"),
    schema_config_path=str(project_root / "config/schema_config.yaml"),
    output_directory=output_dir_path,
)

# Load settings from config and CLI args
CACHE = args.cache
export_as_csv = crossbar_config['settings']['export_csv']
if export_as_csv:
    csv_output_dir.mkdir(parents=True, exist_ok=True)
TEST_MODE = args.test_mode
UPDATE_SCHEMA_DYNAMICALLY = crossbar_config['settings']['update_schema_dynamically']
ORGANISM = crossbar_config['settings']['organism']
USE_EMBEDDINGS = crossbar_config['settings']['use_embeddings']

stats_manager = StatsManager(output_dir_path, stats_dir, TEST_MODE)

if args.stats_only:
    stats_manager.write_stats()
    print(f"Stats written to {stats_dir}")
    sys.exit(0)


def update_schema_with_dynamic_types(schema_path: str, annotation_types: set, feature_types: set):
    """
    Dynamically update schema_config.yaml with discovered annotation and feature types.
    Only needed on first run or when new annotation/feature types are discovered.
    
    Args:
        schema_path: Path to schema_config.yaml
        annotation_types: Set of annotation type names (e.g., 'biophysicochemical_properties_annotation')
        feature_types: Set of feature type names (e.g., 'zinc_finger_feature')
    """
    # Load existing schema
    with open(schema_path, 'r', encoding='utf-8') as f:
        schema = yaml.safe_load(f)
    
    added_types = []
    
    # Add missing annotation types (inherit from functional annotation)
    for ann_type in annotation_types:
        # Convert to schema key format (underscore to space)
        schema_key = ann_type.replace('_', ' ')
        if schema_key not in schema:
            schema[schema_key] = {
                'is_a': 'functional annotation',
                'represented_as': 'node',
                'label_in_input': ann_type,
            }
            added_types.append(schema_key)
    
    # Add missing feature types (inherit from sequence feature)
    for feat_type in feature_types:
        # Convert to schema key format (underscore to space)
        schema_key = feat_type.replace('_', ' ')
        if schema_key not in schema:
            schema[schema_key] = {
                'is_a': 'sequence feature',
                'represented_as': 'node',
                'label_in_input': feat_type,
            }
            added_types.append(schema_key)
    
    if added_types:
        # Write updated schema with blank lines between entries
        yaml_str = yaml.dump(schema, default_flow_style=False, allow_unicode=True, sort_keys=False)
        lines = yaml_str.split('\n')
        formatted_lines = []
        for i, line in enumerate(lines):
            if i > 0 and line and not line[0].isspace() and formatted_lines and formatted_lines[-1]:
                formatted_lines.append('')
            formatted_lines.append(line)

        with open(schema_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(formatted_lines))

        print(f"Schema updated with {len(added_types)} new types:")
        for t in sorted(added_types):
            print(f"  + {t}")
    else:
        print("No new schema types needed.")
    
    return added_types


# Load credentials and adapter configurations from config
drugbank_user = crossbar_config['drugbank']['user']
drugbank_passwd = crossbar_config['drugbank']['password']

# UniProt configuration
uniprot_node_types = [getattr(UniprotNodeType, nt) for nt in crossbar_config['uniprot']['node_types']]
uniprot_node_fields = [getattr(UniprotNodeField, nf) for nf in crossbar_config['uniprot']['node_fields']]
uniprot_edge_types = [getattr(UniprotEdgeType, et) for et in crossbar_config['uniprot']['edge_types']]
uniprot_id_type = [getattr(UniprotIDField, idf) for idf in crossbar_config['uniprot']['id_fields']]


# UniProt SwissProt adapter
def run_uniprot():
    uniprot_adapter = UniprotSwissprot(
        json_path=uniprot_json_path,
        organism=ORGANISM,
        node_types=uniprot_node_types,
        node_fields=uniprot_node_fields,
        edge_types=uniprot_edge_types,
        id_fields=uniprot_id_type,
        test_mode=TEST_MODE,
    )
    if TEST_MODE:
        uniprot_adapter.test_limit = crossbar_config['settings'].get('test_limit', uniprot_adapter.test_limit)

    timestamp = datetime.now(TZ).strftime("%Y-%m-%d %H:%M:%S")
    print(
        f"[{timestamp}] UniprotSwissprot init: organism={ORGANISM}, "
        f"test_mode={TEST_MODE}, test_limit={uniprot_adapter.test_limit}, "
        f"nodes={len(uniprot_node_types)}, edges={len(uniprot_edge_types)}, "
        f"node_fields={len(uniprot_node_fields)}, id_fields={len(uniprot_id_type)}"
    )

    uniprot_adapter.download_uniprot_data(cache=CACHE,
        prott5_embedding_output_path=prott5_embedding_path if USE_EMBEDDINGS else None,
        esm2_embedding_path=esm2_embedding_path if USE_EMBEDDINGS else None,
        nucleotide_transformer_embedding_path=nt_embedding_path if USE_EMBEDDINGS else None)

    # Optionally update schema with dynamically discovered types
    if UPDATE_SCHEMA_DYNAMICALLY:
        print("Updating schema with dynamically discovered annotation/feature types...")
        schema_path = str(project_root / "config/schema_config.yaml")
        uniprot_nodes = list(uniprot_adapter.get_nodes())
        annotation_types = {ntype for _, ntype, _ in uniprot_nodes if ntype.endswith("_annotation")}
        feature_types = {ntype for _, ntype, _ in uniprot_nodes if ntype.endswith("_feature")}
        update_schema_with_dynamic_types(schema_path, annotation_types, feature_types)
        print("Schema updated. Please restart the script for changes to take effect.")
        sys.exit(0)

    uniprot_nodes = list(uniprot_adapter.get_nodes())
    uniprot_edges = list(uniprot_adapter.get_edges())
    # Normalize edge properties per label to avoid BioCypher schema mismatch
    if uniprot_edges:
        cleaned_edges = []
        id_prop_count = 0
        label_keys: dict[str, set[str]] = {}
        for edge in uniprot_edges:
            if edge and len(edge) == 5:
                props = edge[4] if isinstance(edge[4], dict) else {}
                props = dict(props)
                if "id" in props:
                    id_prop_count += 1
                    props.pop("id", None)
                label = edge[3]
                label_keys.setdefault(label, set()).update(props.keys())
                cleaned_edges.append((edge[0], edge[1], edge[2], label, props))
            else:
                cleaned_edges.append(edge)

        normalized_edges = []
        for edge in cleaned_edges:
            if edge and len(edge) == 5:
                label = edge[3]
                props = dict(edge[4]) if isinstance(edge[4], dict) else {}
                for key in label_keys.get(label, set()):
                    if key not in props:
                        props[key] = ""
                normalized_edges.append((edge[0], edge[1], edge[2], label, props))
            else:
                normalized_edges.append(edge)

        if id_prop_count:
            print(f"[WARN] Removed 'id' property from {id_prop_count} UniProt edges.")
        uniprot_edges = normalized_edges

    # Write nodes and edges (includes extended types)
    bc.write_nodes(uniprot_nodes)
    bc.write_edges(uniprot_edges)

    if export_as_csv:
        uniprot_adapter.export_data_to_csv(path=str(csv_output_dir),
                                            node_data=uniprot_nodes,
                                            edge_data=uniprot_edges)
    print(f"SwissProt data exported to CSV successfully.")
    annotation_types = {ntype for _, ntype, _ in uniprot_nodes if ntype.endswith("_annotation")}
    feature_types = {ntype for _, ntype, _ in uniprot_nodes if ntype.endswith("_feature")}
    disease_nodes = [n for n in uniprot_nodes if n[1] == "uniprot_disease"]
    protein_with_keywords = {
        source for _, source, target, label, _ in uniprot_edges if label == "protein_has_keyword"
    }
    print(f"  Annotation types: {len(annotation_types)}")
    print(f"  Feature types: {len(feature_types)}")
    print(f"  Disease nodes: {len(disease_nodes)}")
    print(f"  Proteins with keywords: {len(protein_with_keywords)}")

run_adapter("uniprot", run_uniprot)


# UniProt Keywords (vocabulary with hierarchy and GO mappings)
def run_keywords():
    print("\n" + "="*60)
    print("Loading UniProt Keywords...")
    print("="*60)

    keywords_adapter = UniprotKeywords(
        json_path=crossbar_config['data_paths']['uniprot_keywords_json'],
        node_types=[getattr(KeywordNodeType, nt).value for nt in crossbar_config['keywords']['node_types']],
        edge_types=[getattr(KeywordEdgeType, et).value for et in crossbar_config['keywords']['edge_types']],
        test_mode=TEST_MODE,
    )
    if TEST_MODE:
        keywords_adapter.test_limit = crossbar_config['settings'].get('test_limit', keywords_adapter.test_limit)

    # Write keyword nodes
    keyword_nodes = keywords_adapter.get_nodes()
    bc.write_nodes(keyword_nodes)

    # Write keyword edges (hierarchy and GO mappings)
    keyword_edges = list(keywords_adapter.get_edges())
    if keyword_edges:
        bc.write_edges(keyword_edges)

    print(f"UniProt Keywords written successfully.")

run_adapter("keywords", run_keywords)


# PPI
def run_ppi():
    ppi_adapter = PPI(organism=ORGANISM,
                      output_dir=output_dir_path,
                      export_csv=export_as_csv,
                      test_mode=TEST_MODE)

    ppi_adapter.download_ppi_data(cache=CACHE, n_workers=96)
    ppi_adapter.process_ppi_data()
    bc.write_edges(ppi_adapter.get_ppi_edges())

run_adapter("ppi", run_ppi)


# protein domain
# Note: pypath 0.16.28+ uses rescued data from OmniPath, currently only human (9606) available.
# Now using new API to download all species interpro data
def run_interpro():
    interpro_adapter = InterPro(
        organism=ORGANISM,
        test_mode=TEST_MODE,
        data_dir=interpro_dir_path
    )

    interpro_adapter.download_interpro_data(cache=CACHE)
    interpro_adapter.download_domain_node_data(dom2vec_embedding_path=dom2vec_domain_embedding_path if USE_EMBEDDINGS else None, cache=CACHE)

    bc.write_nodes(interpro_adapter.get_interpro_nodes())
    write_edges_chunked(bc, interpro_adapter.get_interpro_edges())

    if export_as_csv:
        interpro_adapter.export_as_csv(path=output_dir_path)

run_adapter("interpro", run_interpro)


# gene ontology
def run_go():
    go_adapter = GO(
        organism=ORGANISM,
        test_mode=TEST_MODE
    )
    go_adapter.download_go_data(cache=CACHE, anc2vec_embedding_path=anc2vec_go_embedding_path if USE_EMBEDDINGS else None)
    bc.write_nodes(go_adapter.get_go_nodes())
    bc.write_edges(go_adapter.get_go_edges())
    if export_as_csv:
        go_adapter.export_as_csv(path=output_dir_path)

run_adapter("go", run_go)

# drug - with memory optimization
def run_drug():
    drug_adapter = Drug(
        drugbank_user=drugbank_user,
        drugbank_passwd=drugbank_passwd,
        organism=ORGANISM,
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        test_mode=TEST_MODE,
    )
    drug_adapter.download_drug_data(cache=CACHE, selformer_embedding_path=selformer_drug_embedding_path if USE_EMBEDDINGS else None)
    drug_adapter.process_drug_data()
    bc.write_nodes(drug_adapter.get_drug_nodes())
    bc.write_edges(drug_adapter.get_edges())

run_adapter("drug", run_drug)

def run_compound():
    compound_adapter = Compound(
        stitch_organism=ORGANISM,
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        test_mode=TEST_MODE,
    )
    compound_adapter.download_compound_data(cache=CACHE, selformer_embedding_path=selformer_compound_embedding_path if USE_EMBEDDINGS else None)
    compound_adapter.process_compound_data()
    bc.write_nodes(compound_adapter.get_compound_nodes())
    bc.write_edges(compound_adapter.get_cti_edges())

run_adapter("compound", run_compound)

def run_orthology():
    orthology_adapter = Orthology(
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        test_mode=TEST_MODE,
        organism=ORGANISM,
    )
    orthology_adapter.download_orthology_data(cache=CACHE)
    # Note: Orthology only has edges, no nodes
    bc.write_edges(orthology_adapter.get_orthology_edges())

run_adapter("orthology", run_orthology)


def run_disease():
    disease_adapter = Disease(
        drugbank_user=drugbank_user,
        drugbank_passwd=drugbank_passwd,
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        test_mode=TEST_MODE
    )
    disease_adapter.download_disease_data(cache=CACHE,
        doc2vec_embedding_path=doc2vec_disease_embedding_path if USE_EMBEDDINGS else None,
        malacards_json_path=malacards_json_path,
        malacards_related_diseases_json_path=malacards_related_diseases_json_path
        )
    bc.write_nodes(disease_adapter.get_nodes())
    bc.write_edges(disease_adapter.get_edges())

run_adapter("disease", run_disease)


def run_phenotype():
    phenotype_adapter = HPO(
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        test_mode=TEST_MODE
    )
    phenotype_adapter.download_hpo_data(cache=CACHE, cada_embedding_path=cada_phenotype_embedding_path if USE_EMBEDDINGS else None)
    bc.write_nodes(phenotype_adapter.get_nodes())
    bc.write_edges(phenotype_adapter.get_edges())

run_adapter("phenotype", run_phenotype)


def run_pathway():
    # Use crossbar_config kegg_organism, or limit to human in test mode
    kegg_organism = crossbar_config['pathway']['kegg_organism'] if not TEST_MODE else ["hsa"]
    pathway_adapter = Pathway(
        drugbank_user=drugbank_user,
        drugbank_passwd=drugbank_passwd,
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        test_mode=TEST_MODE,
        kegg_organism=kegg_organism,
    )
    pathway_adapter.download_pathway_data(cache=CACHE, biokeen_embedding_path=biokeen_pathway_embedding_path if USE_EMBEDDINGS else None)
    bc.write_nodes(pathway_adapter.get_nodes())
    bc.write_edges(pathway_adapter.get_edges())

run_adapter("pathway", run_pathway)


def run_side_effect():
    side_effect_adapter = SideEffect(
        drugbank_user=drugbank_user,
        drugbank_passwd=drugbank_passwd,
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        test_mode=TEST_MODE
    )
    side_effect_adapter.download_side_effect_data(cache=CACHE)
    bc.write_nodes(side_effect_adapter.get_nodes())
    bc.write_edges(side_effect_adapter.get_edges())

run_adapter("side_effect", run_side_effect)


def run_ec():
    ec_adapter = EC(
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        test_mode=TEST_MODE,
        organism=ORGANISM,
    )
    ec_adapter.download_ec_data(cache=CACHE, rxnfp_embedding_path=rxnfp_ec_embedding_path if USE_EMBEDDINGS else None)
    bc.write_nodes(ec_adapter.get_nodes())
    bc.write_edges(ec_adapter.get_edges())

run_adapter("ec", run_ec)


def run_tfgene():
    tfgene_adapter = TFGene(
        export_csv=export_as_csv,
        output_dir=output_dir_path,
        organism=ORGANISM,
        test_mode=TEST_MODE,
    )
    tfgene_adapter.download_tfgen_data(cache=CACHE)
    bc.write_edges(tfgene_adapter.get_edges())

run_adapter("tfgene", run_tfgene)


# Write import call and other post-processing
bc.write_import_call()
# bc.summary()  # Disabled: ontology summary can hang on multiple inheritance

print("="*80)
print("Script completed successfully!")
print("="*80)
