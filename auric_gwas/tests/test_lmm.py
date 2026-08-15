"""Portable unit tests for the LMM core (no frozen-panel access) + a marked known-answer integration.

Run: python -m pytest auric_gwas/tests/test_lmm.py -v
"""
import numpy as np
import pytest

from auric_gwas import lmm


def _synthetic(n=120, m=200, seed=0):
    rng = np.random.default_rng(seed)
    G = rng.integers(0, 2, size=(n, m)).astype(np.int8)
    # inject a little missingness and one true-effect column correlated with y
    G[rng.random((n, m)) < 0.05] = -1
    y = (G[:, 0] == 1).astype(float) + rng.normal(0, 0.3, n)
    Xcov = np.column_stack([np.ones(n)])
    return G, y, Xcov


def test_grm_unit_mean_diagonal():
    G, _, _ = _synthetic()
    _, K = lmm.build_grm(G)
    assert np.isclose(np.mean(np.diag(K)), 1.0)
    assert np.allclose(K, K.T)


def test_build_grm_rejects_non_int8():
    with pytest.raises(ValueError):
        lmm.build_grm(np.zeros((5, 5), dtype=np.int16))


def test_blind_scan_matches_ols():
    """The lineage-blind Wald scan must equal ordinary least squares on the same design."""
    G, y, Xcov = _synthetic()
    Xc, _ = lmm.build_grm(G)
    out = lmm.scan(Xc, y, Xcov)  # blind
    # reference OLS for a handful of columns: y ~ intercept + g
    for j in [0, 1, 50, 199]:
        X = np.column_stack([Xcov, Xc[:, j]])
        beta_ref, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta_ref
        dof = len(y) - X.shape[1]
        sigma2 = (resid @ resid) / dof
        se_ref = np.sqrt(sigma2 * np.linalg.inv(X.T @ X)[-1, -1])
        assert np.isclose(out["beta"][j], beta_ref[-1], rtol=1e-6, atol=1e-9)
        assert np.isclose(out["se"][j], se_ref, rtol=1e-6, atol=1e-9)


def test_true_effect_is_top_ranked():
    G, y, Xcov = _synthetic()
    Xc, K = lmm.build_grm(G)
    reml = lmm.fit_reml(K, y, Xcov)
    aware = lmm.scan(Xc, y, Xcov, U=reml["U"], evals=reml["evals"], delta=reml["delta"])
    assert int(np.argmin(aware["p"])) == 0  # column 0 is the injected true effect


def test_covariate_design_rank_check():
    import pandas as pd
    df = pd.DataFrame({"cohort": ["a", "a", "b", "b"], "const": ["x", "x", "x", "x"]})
    X, names = lmm.build_covariates(df, ["cohort", "const"])
    assert "intercept" in names and "cohort=b" in names
    assert "const=x" not in names  # constant dropped
    # collinear duplicate must raise
    df2 = pd.DataFrame({"a": ["p", "q", "p", "q"], "b": ["p", "q", "p", "q"]})
    with pytest.raises(ValueError):
        lmm.build_covariates(df2, ["a", "b"])


@pytest.mark.integration
def test_levofloxacin_known_answer():
    """gyrA Ser84Leu must be the rank-1 lineage-aware association (needs the frozen panel)."""
    import os
    pheno = os.environ.get("AURIC_PHENO", "")
    if not (os.path.exists("/fastpool/sausnp/db/genotype_matrix/variant_order.parquet") and pheno and os.path.exists(pheno)):
        pytest.skip("AURIC reference and AURIC_PHENO=<levofloxacin mic table> required")
    from auric_gwas import scan
    out = scan.run(dict(
        pheno_path=pheno, id_col="genome_id",
        pheno_col="mic", pheno_type="binary", r_ge=4, s_le=1, covar_cols=["cohort"],
        filter_col="drug", filter_val="Levofloxacin"), source="panel")
    top = out["results"].sort_values("p_aware").iloc[0]
    assert top["feature_id"] == "group_4352"          # gyrA
    assert out["diagnostics"]["n_case"] == 118
    assert out["diagnostics"]["n_sig_aware"] < 50     # collapse from thousands
