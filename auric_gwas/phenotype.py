"""Generic phenotype + covariate loading — replaces the per-drug hard-coding in the AMR scripts.

A phenotype table is any TSV/CSV with an id column, a phenotype column, and optional covariate columns.
Binary phenotypes may be given as 0/1, R/S, or a MIC with breakpoints (--r-ge/--s-le). Continuous
phenotypes are used as-is (e.g. log2 MIC). Only samples present in BOTH the table and the genotype set
are analysed.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def _read(path: Path) -> pd.DataFrame:
    sep = "\t" if str(path).endswith((".tsv", ".txt")) else ","
    return pd.read_csv(path, sep=sep, dtype=str).convert_dtypes()


def load(pheno_path, *, id_col: str, pheno_col: str, pheno_type: str,
         covar_cols: list[str] | None = None,
         r_ge: float | None = None, s_le: float | None = None,
         filter_col: str | None = None, filter_val: str | None = None) -> dict:
    """Return {samples, y, covar, dropped_intermediate, n_input}.

    pheno_type='binary': y in {0,1}. Accepts 0/1, R/S (case-insensitive first letter), or a numeric
      MIC with r_ge/s_le breakpoints (values strictly between are intermediate and dropped).
    pheno_type='continuous': y is float(pheno_col) as-is; no rows dropped for being intermediate.
    covar is a DataFrame aligned to samples (may be empty).
    """
    df = _read(Path(pheno_path))
    for c in [id_col, pheno_col, *(covar_cols or [])]:
        if c not in df.columns:
            raise KeyError(f"column {c!r} not in phenotype table (has {list(df.columns)})")
    if filter_col is not None:
        if filter_col not in df.columns:
            raise KeyError(f"filter column {filter_col!r} not in phenotype table")
        df = df[df[filter_col].astype(str) == str(filter_val)].copy()
    df = df[df[pheno_col].notna()].copy()
    df[id_col] = df[id_col].astype(str)
    n_input = len(df)
    dropped_intermediate = 0

    if pheno_type == "continuous":
        y = df[pheno_col].astype(float).to_numpy()
    elif pheno_type == "binary":
        raw = df[pheno_col]
        if r_ge is not None and s_le is not None:
            mic = raw.astype(float).to_numpy()
            code = np.where(mic >= r_ge, 1, np.where(mic <= s_le, 0, -1))
        else:
            up = raw.astype(str).str.strip().str.upper().str[0]
            code = np.where(up.isin(["R"]), 1, np.where(up.isin(["S"]), 0,
                     np.where(up.isin(["1"]), 1, np.where(up.isin(["0"]), 0, -1))))
        keep = code >= 0
        dropped_intermediate = int((~keep).sum())
        df = df[keep].copy()
        y = code[keep].astype(float)
    else:
        raise ValueError("pheno_type must be 'binary' or 'continuous'")

    df = df.drop_duplicates(id_col).reset_index(drop=True)
    # y was aligned to df before dedup; re-derive after dedup to stay aligned
    if pheno_type == "continuous":
        y = df[pheno_col].astype(float).to_numpy()
    else:
        raw = df[pheno_col]
        if r_ge is not None and s_le is not None:
            mic = raw.astype(float).to_numpy()
            y = np.where(mic >= r_ge, 1.0, 0.0)
        else:
            up = raw.astype(str).str.strip().str.upper().str[0]
            y = np.where(up.isin(["R", "1"]), 1.0, 0.0)

    covar = df[covar_cols].copy() if covar_cols else pd.DataFrame(index=df.index)
    return {
        "samples": df[id_col].tolist(),
        "y": y,
        "covar": covar,
        "dropped_intermediate": dropped_intermediate,
        "n_input": n_input,
    }
