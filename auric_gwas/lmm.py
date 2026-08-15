"""EMMAX-style lineage-aware association — the reusable core of AURIC-associate.

Extracted from the four copy-pasted implementations in scripts/analysis/amr/ (the reference is
assoc_grm_recomb_sweep.py, the only copy whose scan() takes its rotation basis as arguments). Pure
NumPy/SciPy on arrays already in memory; no I/O, no globals, no phenotype-specific logic.

The method (Kang et al., Nat Genet 2010): estimate one variance-component ratio (delta) from the
kinship-only null by REML, then hold it fixed and test every variant with GLS weights — the
"eXpedited" shortcut that avoids a per-variant mixed-model refit.

Genotype convention throughout: int8 (n_genomes, m_variants), 0 = called reference OR a different alt,
1 = called this alt, -1 = MISSING (never folded to reference).
"""
from __future__ import annotations

import numpy as np
from scipy import linalg, optimize, stats


def build_grm(G: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (Xc, K): column-centred genotypes and the genomic relatedness matrix.

    Each column is centred on its called-genotype allele frequency; missing cells are set to exactly
    the site mean (0 after centring) so they contribute nothing. K = Xc Xc' normalised to unit mean
    diagonal.
    """
    if G.dtype != np.int8:
        raise ValueError(f"G must be int8 (pattern hashing and encoding depend on it); got {G.dtype}")
    miss = G == -1
    an = (~miss).sum(0)
    ac = (G == 1).sum(0)
    af = ac / np.maximum(an, 1)
    Xc = G.astype(np.float64) - af          # missing rows still hold -1 here...
    Xc[miss] = 0.0                          # ...now exactly the site mean
    K = Xc @ Xc.T
    K /= np.mean(np.diag(K))
    return Xc, K


def fit_reml(K: np.ndarray, y: np.ndarray, Xcov: np.ndarray) -> dict:
    """Fit the variance-component ratio delta under the null (covariates only, no genotype).

    Returns a dict with delta, h2 = 1/(1+delta), the eigenbasis (U, evals) of K, and grid_edge — a
    flag that delta hit the search bound, i.e. h2 saturated (h2->1 means over-correction is likely and
    the aware p-values must not be taken as final; use a PC-fixed-effects fallback).
    """
    n, p_cov = Xcov.shape
    evals, U = linalg.eigh(K)
    evals = np.clip(evals, 0, None)
    Uty = U.T @ y
    UtX = U.T @ Xcov
    logdetXtX = np.linalg.slogdet(Xcov.T @ Xcov)[1]

    def neg_reml(log_delta: float) -> float:
        delta = np.exp(log_delta)
        w = 1.0 / (evals + delta)
        XtWX = UtX.T @ (w[:, None] * UtX)
        beta = np.linalg.solve(XtWX, UtX.T @ (w * Uty))
        r = Uty - UtX @ beta
        rss = float(np.sum(w * r * r))
        dfr = n - p_cov
        return 0.5 * (dfr * np.log(2 * np.pi * rss / dfr) + np.sum(np.log(evals + delta)) + dfr
                      + logdetXtX - np.linalg.slogdet(XtWX)[1])

    grid = np.linspace(-10, 10, 201)
    vals = np.array([neg_reml(g) for g in grid])
    g0 = grid[int(np.argmin(vals))]
    res = optimize.minimize_scalar(neg_reml, bounds=(g0 - 0.5, g0 + 0.5), method="bounded")
    log_delta = float(res.x)
    delta = float(np.exp(log_delta))
    return {
        "delta": delta,
        "h2": 1.0 / (1.0 + delta),
        "U": U,
        "evals": evals,
        "grid_edge": bool(log_delta <= -9.5 or log_delta >= 9.5),
    }


def scan(Xc: np.ndarray, y: np.ndarray, Xcov: np.ndarray,
         U: np.ndarray | None = None, evals: np.ndarray | None = None,
         delta: float | None = None, block: int = 20000) -> dict:
    """Block-vectorised Wald scan of every column of Xc.

    U/evals/delta given -> lineage-aware GLS in the rotated basis (weights 1/(evals+delta)).
    U is None            -> lineage-blind ordinary linear model (unit weights, unrotated).
    A single null fit is reused for every variant. Returns beta, se, p (two-sided t), chi2, and the
    genomic-control lambda.
    """
    n, m = Xc.shape
    p_cov = Xcov.shape[1]
    aware = U is not None
    if aware:
        Ry, RX, w = U.T @ y, U.T @ Xcov, 1.0 / (evals + delta)
    else:
        Ry, RX, w = y, Xcov, np.ones(n)

    XtWX = RX.T @ (w[:, None] * RX)
    Ainv = np.linalg.inv(XtWX)
    XtWy = RX.T @ (w * Ry)
    beta0 = Ainv @ XtWy
    r0 = Ry - RX @ beta0
    rss0 = float(np.sum(w * r0 * r0))
    dfr = n - p_cov - 1

    beta = np.empty(m)
    se = np.empty(m)
    for s in range(0, m, block):
        e = min(s + block, m)
        Xb = Xc[:, s:e]
        Rg = Xb if not aware else (U.T @ Xb)
        wRg = w[:, None] * Rg
        aa = np.einsum("ij,ij->j", Rg, wRg)
        bb = RX.T @ wRg
        cc = (w * Ry) @ Rg
        Ab = Ainv @ bb
        S = aa - np.einsum("ij,ij->j", bb, Ab)
        with np.errstate(divide="ignore", invalid="ignore"):
            bg = (cc - Ab.T @ XtWy) / np.where(S > 0, S, np.nan)
            rss_full = rss0 - bg * bg * S
            se[s:e] = np.sqrt((rss_full / dfr) / S)
        beta[s:e] = bg

    with np.errstate(divide="ignore", invalid="ignore"):
        t = beta / se
    p = 2 * stats.t.sf(np.abs(t), dfr)
    chi2 = t ** 2
    finite = np.isfinite(chi2)
    lam = float(np.median(chi2[finite]) / stats.chi2.ppf(0.5, 1)) if finite.any() else np.nan
    return {"beta": beta, "se": se, "p": p, "chi2": chi2, "lambda_gc": lam}


def n_patterns(G: np.ndarray) -> int:
    """Distinct genotype patterns (missingness included, matching pyseer) for pattern-Bonferroni.

    Requires int8: the void-view hashing assumes exactly one byte per sample.
    """
    if G.dtype != np.int8:
        raise ValueError("pattern hashing requires int8 G (one byte per sample)")
    pat = np.ascontiguousarray(G.T).view(np.dtype((np.void, G.shape[0])))
    return int(np.unique(pat).size)


def build_covariates(cov_df, cols) -> tuple[np.ndarray, list[str]]:
    """Design matrix with an explicit intercept, k-1 dummies per categorical column, dropped constant
    columns, and a rank check. cov_df is a pandas DataFrame aligned to the sample order.

    Returns (Xcov, names). Raises if the design is rank-deficient after expansion.
    """
    import pandas as pd
    n = len(cov_df)
    mats = [np.ones((n, 1))]
    names = ["intercept"]
    for c in cols:
        s = cov_df[c]
        if pd.api.types.is_numeric_dtype(s) and s.nunique(dropna=True) > 2:
            v = s.astype(float).to_numpy()
            if np.nanstd(v) == 0:
                continue  # constant, carries no information
            mats.append(v.reshape(-1, 1))
            names.append(c)
        else:
            levels = sorted(x for x in s.dropna().unique())
            for lv in levels[1:]:  # first level is the reference
                mats.append((s == lv).astype(float).to_numpy().reshape(-1, 1))
                names.append(f"{c}={lv}")
    X = np.column_stack(mats)
    if np.linalg.matrix_rank(X) < X.shape[1]:
        raise ValueError(f"covariate design is rank-deficient ({names}); drop collinear/constant covariates")
    return X, names
