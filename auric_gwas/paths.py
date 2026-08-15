"""Frozen-panel locations, made configurable for distribution, plus the read-only guard.

The AURIC reference artifacts (profile HMMs, coordinate table, frozen column order, DuckDB catalogue)
are large and are distributed separately from this package (a Zenodo bundle; see the README). Their
location is resolved, in order:
  1. the AURIC_HOME environment variable, or a call to configure(auric_home=...);
  2. the built-in default /fastpool/sausnp (the development mount).
Under AURIC_HOME the layout is db/ (sausnp.duckdb, db_version.json, coordinate_table.parquet,
genotype_matrix/) and profiles/ + profiles_masked/.

Every path here is an INPUT, opened read-only. The frozen panel must never receive a write; two panel
aggregators derive the callable-mask denominator by counting subdirectories under the predict dir, and
one does CREATE OR REPLACE TABLE on the catalogue. assert_not_panel() is the runtime guard.
"""
from __future__ import annotations

import os
from pathlib import Path

_DEFAULT_HOME = Path("/fastpool/sausnp")


def auric_home() -> Path:
    return Path(os.environ.get("AURIC_HOME", str(_DEFAULT_HOME)))


def configure(auric_home: str | os.PathLike) -> None:
    """Point the tool at a local AURIC reference directory (overrides AURIC_HOME for this process)."""
    os.environ["AURIC_HOME"] = str(auric_home)


# --- resolved locations (call the accessor, or use the module-level snapshots below) ---------------
def db() -> Path:                return auric_home() / "db" / "sausnp.duckdb"
def db_version_json() -> Path:   return auric_home() / "db" / "db_version.json"
def matrix() -> Path:            return auric_home() / "db" / "genotype_matrix"
def variant_order() -> Path:     return matrix() / "variant_order.parquet"
def genome_order() -> Path:      return matrix() / "genome_order.parquet"
def genotypes_zarr() -> Path:    return matrix() / "genotypes.zarr"
def coordinate_table() -> Path:  return auric_home() / "db" / "coordinate_table.parquet"
def profiles() -> Path:          return auric_home() / "profiles"
def profiles_masked() -> Path:   return auric_home() / "profiles_masked"


# Backwards-compatible module-level snapshots (resolved at import against the current AURIC_HOME).
# Prefer the accessors above in long-lived processes where AURIC_HOME may change.
DB = db()
DB_VERSION_JSON = db_version_json()
MATRIX = matrix()
VARIANT_ORDER = variant_order()
GENOME_ORDER = genome_order()
GENOTYPES_ZARR = genotypes_zarr()
COORDINATE_TABLE = coordinate_table()
PROFILES = profiles()
PROFILES_MASKED = profiles_masked()


def _forbidden_roots() -> tuple[Path, ...]:
    h = auric_home()
    return (h / "db", h / "profiles", h / "profiles_masked")


def assert_not_panel(path) -> Path:
    """Raise if `path` would write under any frozen-panel root. Call before opening any write target."""
    p = Path(path).resolve()
    for root in _forbidden_roots():
        try:
            p.relative_to(root.resolve())
        except (ValueError, OSError):
            continue
        raise PermissionError(f"refusing to write under the frozen panel: {p} is inside {root}")
    return p


def db_version(path: Path | None = None) -> str:
    """Read the catalogue version fresh (never cache — a concurrent writer can bump it)."""
    import json
    return json.loads(Path(path or db_version_json()).read_text()).get("db_version", "unknown")


def set_thread_caps(n: int = 4) -> None:
    """Cap BLAS/OMP threads. MUST run before numpy is imported by the caller for it to bite fully."""
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ.setdefault(var, str(n))
