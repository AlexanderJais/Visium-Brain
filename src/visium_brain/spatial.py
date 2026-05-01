"""Spatial-aware analyses on Visium HD bins.

The neighbor graph is built **per sample** and then concatenated, so
edges never cross between sections. Downstream analyses available
through ``run()``:

* Moran's I for spatial autocorrelation of the top HVGs
  (``squidpy.gr.spatial_autocorr``).
* Neighborhood-enrichment z-scores between annotated cell types
  (``squidpy.gr.nhood_enrichment``).

Co-occurrence is intentionally not in ``run()`` because at HD scale it's
slow; call ``squidpy.gr.co_occurrence`` directly if you need it.
"""

from __future__ import annotations

import logging
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _build_graph_per_sample(adata: ad.AnnData, n_neighbors: int = 6) -> None:
    """Build a spatial graph per sample, then concatenate into adata.obsp.

    squidpy's `gr.spatial_neighbors` operates on the whole AnnData. To avoid
    bridging samples, we build per sample and stitch the sparse matrices.

    Raises ``RuntimeError`` if no sample has enough bins to produce a
    graph (previously this would surface as a confusing
    ``np.concatenate([])`` ValueError).
    """
    import scipy.sparse as sp
    import squidpy as sq

    obs_idx = pd.Series(np.arange(adata.n_obs), index=adata.obs_names)
    rows_d: list = []; cols_d: list = []; vals_d: list = []
    rows_w: list = []; cols_w: list = []; vals_w: list = []
    n_samples_used = 0
    for sid in adata.obs["sample_id"].cat.categories:
        mask = (adata.obs["sample_id"] == sid).values
        sub = adata[mask].copy()
        if sub.n_obs < n_neighbors + 1:
            logger.warning(
                "[spatial] %s: %d bins < n_neighbors+1=%d; skipping",
                sid, sub.n_obs, n_neighbors + 1,
            )
            continue
        sq.gr.spatial_neighbors(sub, coord_type="generic", n_neighs=n_neighbors)
        local = obs_idx.loc[sub.obs_names].values
        d = sub.obsp["spatial_distances"].tocoo()
        w = sub.obsp["spatial_connectivities"].tocoo()
        rows_d.append(local[d.row]); cols_d.append(local[d.col]); vals_d.append(d.data)
        rows_w.append(local[w.row]); cols_w.append(local[w.col]); vals_w.append(w.data)
        n_samples_used += 1

    if n_samples_used == 0:
        raise RuntimeError(
            f"No sample had enough bins (>= n_neighbors+1 = {n_neighbors + 1}) "
            f"to build a spatial graph. Lower spatial.n_neighbors or relax QC."
        )

    n = adata.n_obs
    distances = sp.csr_matrix(
        (np.concatenate(vals_d), (np.concatenate(rows_d), np.concatenate(cols_d))),
        shape=(n, n),
    )
    connectivities = sp.csr_matrix(
        (np.concatenate(vals_w), (np.concatenate(rows_w), np.concatenate(cols_w))),
        shape=(n, n),
    )
    adata.obsp["spatial_distances"] = distances
    adata.obsp["spatial_connectivities"] = connectivities
    adata.uns["spatial_neighbors"] = {
        "connectivities_key": "spatial_connectivities",
        "distances_key": "spatial_distances",
        "params": {"n_neighbors": n_neighbors, "coord_type": "generic"},
    }


def morans_i(adata: ad.AnnData, n_genes: int = 200) -> pd.DataFrame:
    """Compute Moran's I for the top HVGs."""
    import squidpy as sq

    if "highly_variable" in adata.var:
        # Filter to HVGs first, then sort by rank if available. Previously
        # the seurat_v3 path sorted the *full* var by highly_variable_rank
        # (NaN on non-HVGs) and relied on NaN-last sort order to push the
        # right rows to the top -- correct in practice, but accidental.
        hv = adata.var[adata.var["highly_variable"]]
        if "highly_variable_rank" in hv.columns:
            hv = hv.sort_values("highly_variable_rank")
        genes = hv.index.tolist()[:n_genes]
    else:
        genes = adata.var_names[:n_genes].tolist()

    sq.gr.spatial_autocorr(adata, genes=genes, mode="moran")
    df = adata.uns["moranI"].copy()
    # squidpy stores the result with gene names as the index but no
    # index.name; set it so morans_i.csv has a real "gene" header
    # instead of an empty leading column.
    df.index.name = "gene"
    return df


def neighborhood_enrichment(adata: ad.AnnData, cluster_key: str = "cell_type") -> pd.DataFrame:
    import squidpy as sq

    sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)
    z = adata.uns[f"{cluster_key}_nhood_enrichment"]["zscore"]
    cats = adata.obs[cluster_key].cat.categories
    df = pd.DataFrame(z, index=cats, columns=cats)
    df.index.name = cluster_key
    return df


def run(adata: ad.AnnData, cfg: dict[str, Any]) -> dict[str, pd.DataFrame]:
    sp_cfg = cfg["spatial"]
    _build_graph_per_sample(adata, n_neighbors=sp_cfg.get("n_neighbors", 6))
    out: dict[str, pd.DataFrame] = {}
    out["morans_i"] = morans_i(adata, n_genes=sp_cfg.get("morans_i_n_genes", 200))
    if sp_cfg.get("neighborhood_enrichment", True) and "cell_type" in adata.obs:
        out["neighborhood_enrichment"] = neighborhood_enrichment(adata, "cell_type")
    return out
