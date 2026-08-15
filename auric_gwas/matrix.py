"""Read-only extraction of a genotype submatrix from the frozen panel or a cohort matrix.

Two sources, one interface:
  - panel: stream /fastpool/sausnp/db/genotype_matrix/genotypes.zarr for a list of panel genome_ids
           (used for known-answer validation and for any analysis over catalogue genomes).
  - cohort: read a matrix a user built with `auric-gwas genotype` (same frozen column order).

Column order is authoritative and comes from variant_order.parquet; col_index is asserted to be the
identity 0..n-1. Candidate columns for association are pass & ~family_masked & (MAC-in-sample >= min).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import zarr

from . import paths


def load_variant_order(matrix: Path = paths.MATRIX) -> pd.DataFrame:
    matrix = Path(matrix)
    vo = pd.read_parquet(matrix / "variant_order.parquet")
    if not np.array_equal(vo["col_index"].to_numpy(), np.arange(len(vo))):
        raise ValueError("variant_order.col_index is not the identity 0..n-1; column order is unsafe")
    return vo


def panel_row_map(matrix: Path = paths.MATRIX) -> dict:
    go = pd.read_parquet(Path(matrix) / "genome_order.parquet", columns=["genome_id", "row_index"])
    return dict(zip(go.genome_id.astype(str), go.row_index))


def extract_submatrix(genome_ids, *, source: str = "panel", matrix: Path = paths.MATRIX,
                      min_mac: int = 5, max_miss: float = 0.10,
                      candidate_mask: np.ndarray | None = None):
    matrix = Path(matrix)
    """Return (G, meta) for the requested genomes.

    G: int8 (n_genomes, m_kept) submatrix over the association-candidate columns.
    meta: the variant_order rows for the kept columns, plus AC/AN/MAC/missfrac in this sample.
    source='panel' streams the frozen zarr; source='cohort' reads matrix/genotypes.zarr built by the
    genotype step (same column order). Read-only in both cases.
    """
    vo = load_variant_order(matrix)
    z = zarr.open_group(str(matrix / "genotypes.zarr"), mode="r")["genotypes"]

    if source == "panel":
        row_of = panel_row_map(matrix)
        missing = [g for g in genome_ids if str(g) not in row_of]
        if missing:
            raise KeyError(f"{len(missing)} genome_ids are not panel rows (e.g. {missing[:3]}); "
                           "use source='cohort' for user-genotyped genomes")
        rows = np.array([row_of[str(g)] for g in genome_ids])
    elif source == "cohort":
        go = pd.read_parquet(matrix / "genome_order.parquet", columns=["genome_id", "row_index"])
        row_of = dict(zip(go.genome_id.astype(str), go.row_index))
        rows = np.array([row_of[str(g)] for g in genome_ids])
    else:
        raise ValueError("source must be 'panel' or 'cohort'")

    if candidate_mask is None:
        candidate_mask = (vo["pass"].astype(bool) & ~vo["family_masked"].astype(bool)).to_numpy()
    cand_idx = np.where(candidate_mask)[0]

    n = len(rows)
    G = np.empty((n, len(cand_idx)), dtype=np.int8)
    step = 4096 * 16
    nv = z.shape[1]
    w = 0
    for c in range(0, nv, step):
        e = min(c + step, nv)
        blk = z[:, c:e][rows]
        msk = candidate_mask[c:e]
        k = int(msk.sum())
        if k:
            G[:, w:w + k] = blk[:, msk]
            w += k
    assert w == len(cand_idx)

    miss = G == -1
    an = n - miss.sum(0)
    ac = (G == 1).sum(0)
    mac = np.minimum(ac, an - ac)
    missfrac = miss.sum(0) / n
    keep = (mac >= min_mac) & (missfrac <= max_miss)
    ki = np.where(keep)[0]
    G = np.ascontiguousarray(G[:, ki])
    meta = vo.iloc[cand_idx[ki]].reset_index(drop=True)
    meta["AC_sample"] = ac[ki]
    meta["AN_sample"] = an[ki]
    meta["MAC_sample"] = mac[ki]
    meta["missfrac_sample"] = missfrac[ki]
    return G, meta
