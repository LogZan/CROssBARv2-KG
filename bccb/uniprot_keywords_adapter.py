"""Thin wrapper to use the shared UniProt Keywords adapter via symlink."""

from __future__ import annotations

import importlib.util
import os
import sys
import types

SHARED_DIR = os.path.join(os.path.dirname(__file__), "_shared_uniprot_keywords")
PKG_NAME = "bccb._shared_uniprot_keywords"

def _ensure_pkg(name: str, path: str) -> None:
    if name in sys.modules:
        return
    pkg = types.ModuleType(name)
    pkg.__path__ = [path]
    sys.modules[name] = pkg


def _load_module(name: str, path: str, package: str | None = None):
    spec = importlib.util.spec_from_file_location(
        name, path, submodule_search_locations=[SHARED_DIR]
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load module {name} from {path}")
    module = importlib.util.module_from_spec(spec)
    if package is not None:
        module.__package__ = package
    spec.loader.exec_module(module)
    return module


_ensure_pkg(PKG_NAME, SHARED_DIR)
_enums_mod = _load_module(
    f"{PKG_NAME}.enums",
    os.path.join(SHARED_DIR, "enums.py"),
    package=PKG_NAME,
)
_prev_enums_mod = sys.modules.get("enums")
sys.modules["enums"] = _enums_mod
try:
    _adapter_mod = _load_module(
        f"{PKG_NAME}.adapter",
        os.path.join(SHARED_DIR, "adapter.py"),
        package=PKG_NAME,
    )
finally:
    if _prev_enums_mod is None:
        sys.modules.pop("enums", None)
    else:
        sys.modules["enums"] = _prev_enums_mod

UniprotKeywords = _adapter_mod.UniprotKeywords  # noqa: E402
UniprotKeywordEnumMeta = _enums_mod.UniprotKeywordEnumMeta  # noqa: E402
UniprotKeywordNodeType = _enums_mod.UniprotKeywordNodeType  # noqa: E402
UniprotKeywordEdgeType = _enums_mod.UniprotKeywordEdgeType  # noqa: E402
UniprotKeywordDocField = _enums_mod.UniprotKeywordDocField  # noqa: E402
KeywordNodeType = UniprotKeywordNodeType  # noqa: E402
KeywordEdgeType = UniprotKeywordEdgeType  # noqa: E402

__all__ = [
    "UniprotKeywords",
    "UniprotKeywordEnumMeta",
    "UniprotKeywordNodeType",
    "UniprotKeywordEdgeType",
    "UniprotKeywordDocField",
    "KeywordNodeType",
    "KeywordEdgeType",
]
