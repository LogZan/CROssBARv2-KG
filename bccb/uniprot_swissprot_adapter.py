"""Thin wrapper to use the shared UniProt SwissProt adapter via symlink."""

from __future__ import annotations

import os
import sys

SHARED_DIR = os.path.join(os.path.dirname(__file__), "_shared_uniprot_swissprot")
if SHARED_DIR not in sys.path:
    sys.path.insert(0, SHARED_DIR)

from adapter import UniprotSwissprot  # noqa: E402
from enums import (  # noqa: E402
    UniprotEnumMeta,
    UniprotNodeType,
    UniprotNodeField,
    UniprotEdgeType,
    UniprotIDField,
    UniprotDocField,
)

__all__ = [
    "UniprotSwissprot",
    "UniprotEnumMeta",
    "UniprotNodeType",
    "UniprotNodeField",
    "UniprotEdgeType",
    "UniprotIDField",
    "UniprotDocField",
]
