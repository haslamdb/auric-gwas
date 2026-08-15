"""Read-only annotation join over the frozen AURIC catalogue (db 0.6.3+).

annotate(variant_ids)      -> one fully annotated row per variant_id, in request order
annotate_matrix_cols(cols) -> same, addressed by genotype-matrix column index
col_index_for(variant_ids) -> variant_id -> zarr column index (NA = not in the matrix)

Every join is LEFT (coding/noncoding/features/callable_mask are each incomplete over the full variant
set); loci are addressed by feature_id, never gene_symbol (NULL for 63.5% of features, non-unique).
Opens the catalogue with read_only=True and writes nothing.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from . import paths

_ANNOT_SQL = """
SELECT
    v.variant_id, v.feature_id,
    f.gene_symbol, f.product, f.nctc8325_locus_tag, f.usa300_locus_tag,
    f.amr_gene AS feature_amr_gene, f.amr_symbol, f.amr_class,
    f.essential, f.vfdb_category, f.mge_overlap, f.cog_category, f.kegg_ko,
    f.assignment_unreliable,
    v.match_state, v.codon_pos, v.nctc8325_pos, ct.usa300_pos,
    coalesce(ct.feature_type, 'family_masked') AS feature_type,
    v.ref, v.alt,
    coalesce(c.consequence, n.region_class) AS effect,
    c.consequence AS coding_consequence,
    c.aa_ref, c.aa_position, c.aa_alt,
    CASE WHEN c.aa_ref IS NOT NULL
         THEN c.aa_ref || CAST(c.aa_position AS VARCHAR) || c.aa_alt END AS aa_change,
    c.codon_ref, c.codon_alt, c.is_4fold_degenerate,
    c.esm2_llr, c.plddt, c.rsa, c.secondary_structure,
    n.region_class, n.flanking_gene_up, n.flanking_gene_down,
    n.in_minus10, n.in_minus35, n.in_rbs, n.in_tfbs, n.tfbs_regulator,
    v.AC, v.AN, v.AF, v.AF_weighted, v.MAC, v.singleton, v.n_alts_at_site,
    v.AF_by_CC, v.n_CC_observed, v.n_CC_fixed, v.n_CC_intermediate, v.n_CC_rare,
    v.n_origins, v.consistency_index, v.af_matched_pctile, v.in_recombinant_tract_frac,
    v."pass", v.family_masked, v.high_confidence, v.hc_readval,
    v.homopolymer, v.is_transition, v.amr_gene AS variant_amr_gene,
    cm.callable, cm.call_frac, cm.exclusion_reason
FROM variants v
LEFT JOIN features f                      ON f.feature_id = v.feature_id
LEFT JOIN variant_annotations_coding c    ON c.variant_id = v.variant_id
LEFT JOIN variant_annotations_noncoding n ON n.variant_id = v.variant_id
LEFT JOIN callable_mask cm
       ON cm.feature_id = v.feature_id AND cm.match_state = v.match_state AND cm.codon_pos = v.codon_pos
LEFT JOIN read_parquet('{coord}') ct
       ON ct.feature_id = v.feature_id AND ct.match_state = v.match_state AND ct.codon_pos = v.codon_pos
JOIN want w ON w.variant_id = v.variant_id
"""


def _con(db: Path = paths.DB):
    import duckdb
    con = duckdb.connect(str(db), read_only=True)
    con.execute("SET threads TO 4")
    return con


def annotate(variant_ids, db: Path = paths.DB, con=None) -> pd.DataFrame:
    want = pd.DataFrame({"variant_id": pd.unique(pd.Series(list(variant_ids), dtype=object))})
    own = con is None
    con = con or _con(db)
    try:
        con.register("want", want)
        out = con.execute(_ANNOT_SQL.format(coord=paths.COORDINATE_TABLE)).df()
    finally:
        con.unregister("want")
        if own:
            con.close()
    return want.merge(out, on="variant_id", how="left")


def col_index_for(variant_ids, matrix: Path = paths.MATRIX) -> pd.Series:
    vo = pd.read_parquet(matrix / "variant_order.parquet", columns=["col_index", "variant_id"])
    m = dict(zip(vo.variant_id, vo.col_index))
    return pd.Series([m.get(v) for v in variant_ids], index=list(variant_ids), dtype="Int64")


def annotate_matrix_cols(col_indices, db: Path = paths.DB, matrix: Path = paths.MATRIX) -> pd.DataFrame:
    vo = pd.read_parquet(matrix / "variant_order.parquet", columns=["col_index", "variant_id"])
    sel = vo.set_index("col_index").loc[list(col_indices), "variant_id"]
    out = annotate(sel.tolist(), db=db)
    out.insert(0, "col_index", list(col_indices))
    return out
