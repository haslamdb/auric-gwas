"""Orchestrate one association: submatrix + phenotype + covariates -> blind & aware ranking.

Ties matrix.extract_submatrix, phenotype.load, lmm.{build_grm,fit_reml,scan}, and annotate together.
Returns a results DataFrame (one row per tested variant) plus a diagnostics dict. Read-only throughout.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import lmm, matrix, phenotype, annotate as annot, paths


def run(pheno_spec: dict, *, source: str = "panel", matrix_dir=paths.MATRIX,
        min_mac: int = 5, max_miss: float = 0.10, annotate_top: int = 200) -> dict:
    """pheno_spec is the kwargs dict for phenotype.load. Returns {results, diagnostics}."""
    paths.assert_bundle_compatible()
    ph = phenotype.load(**pheno_spec)
    samples, y, covar = ph["samples"], ph["y"], ph["covar"]

    # keep only phenotyped samples that are in the genotype source
    G_all, meta = matrix.extract_submatrix(
        samples, source=source, matrix=matrix_dir, min_mac=min_mac, max_miss=max_miss)
    n, m = G_all.shape

    # covariate design (explicit intercept + k-1 dummies, rank-checked)
    covar_cols = list(covar.columns)
    Xcov, cov_names = lmm.build_covariates(covar.reset_index(drop=True), covar_cols)

    Xc, K = lmm.build_grm(G_all)
    reml = lmm.fit_reml(K, y, Xcov)
    blind = lmm.scan(Xc, y, Xcov)
    aware = lmm.scan(Xc, y, Xcov, U=reml["U"], evals=reml["evals"], delta=reml["delta"])

    npat = lmm.n_patterns(G_all)
    bonf = 0.05 / npat

    res = meta.copy()
    res["p_blind"] = blind["p"]
    res["p_aware"] = aware["p"]
    res["beta_aware"] = aware["beta"]
    res["rank_blind"] = pd.Series(blind["p"]).rank(method="min").astype(int)
    res["rank_aware"] = pd.Series(aware["p"]).rank(method="min").astype(int)
    res = res.sort_values("p_aware").reset_index(drop=True)

    # annotate the top hits (read-only join); leave the rest un-annotated for size
    top_ids = res["variant_id"].head(annotate_top).tolist()
    ann = annot.annotate(top_ids)
    keep_cols = ["variant_id", "gene_symbol", "product", "feature_type", "nctc8325_pos",
                 "effect", "aa_change", "amr_symbol", "amr_class", "esm2_llr",
                 "n_origins", "in_recombinant_tract_frac", "AF", "n_CC_fixed"]
    res = res.merge(ann[[c for c in keep_cols if c in ann.columns]], on="variant_id", how="left")

    diagnostics = {
        "db_version": paths.db_version(),
        "n_samples": int(n),
        "n_case": int((y == 1).sum()) if set(np.unique(y)) <= {0.0, 1.0} else None,
        "n_control": int((y == 0).sum()) if set(np.unique(y)) <= {0.0, 1.0} else None,
        "n_variants_tested": int(m),
        "n_patterns": int(npat),
        "bonferroni": bonf,
        "covariates": cov_names,
        "delta": reml["delta"],
        "h2_kinship": reml["h2"],
        "reml_grid_edge": reml["grid_edge"],
        "lambda_blind": blind["lambda_gc"],
        "lambda_aware": aware["lambda_gc"],
        "n_sig_blind": int((blind["p"] < bonf).sum()),
        "n_sig_aware": int((aware["p"] < bonf).sum()),
        "dropped_intermediate": ph["dropped_intermediate"],
        "warnings": _warnings(reml, aware["lambda_gc"], blind["lambda_gc"]),
    }
    return {"results": res, "diagnostics": diagnostics}


def _warnings(reml: dict, lam_aware: float, lam_blind: float) -> list[str]:
    w = []
    if reml["grid_edge"] or reml["h2"] >= 0.999:
        w.append(f"h2 saturated at {reml['h2']:.4f} (delta at the search bound): the mixed model is "
                 "likely over-correcting; treat the aware ranking as conservative, not final, and "
                 "consider a PC-fixed-effects score test (see scripts/analysis/amr/pc_score_test.py).")
    if lam_aware < 0.5:
        w.append(f"lambda_aware={lam_aware:.3f} is deflated (<0.5): over-correction — real signals "
                 "beyond the largest-effect determinant may be attenuated.")
    if lam_blind > 2.0:
        w.append(f"lambda_blind={lam_blind:.2f} is inflated: the unadjusted scan is dominated by "
                 "population structure, as expected; interpret the aware scan.")
    return w
