"""Thin wrapper to use the shared UniProt Keywords adapter via symlink."""

from __future__ import annotations

import os
import sys

SHARED_DIR = os.path.join(os.path.dirname(__file__), "_shared_uniprot_keywords")
if SHARED_DIR not in sys.path:
    sys.path.insert(0, SHARED_DIR)

from adapter import UniprotKeywords  # noqa: E402
from enums import (  # noqa: E402
    UniprotKeywordEnumMeta,
    UniprotKeywordNodeType,
    UniprotKeywordEdgeType,
    UniprotKeywordDocField,
)

__all__ = [
    "UniprotKeywords",
    "UniprotKeywordEnumMeta",
    "UniprotKeywordNodeType",
    "UniprotKeywordEdgeType",
    "UniprotKeywordDocField",
]
