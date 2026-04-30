"""Neighborhood graph, UMAP, and Leiden clustering.

Supports per-level configuration via ``cfg['clustering']['l1']`` and
``cfg['clustering']['l2']`` (number of PCs, neighbors, resolution). If a
flat ``cfg['clustering']`` block is provided (legacy), it is treated as
the L1 config.
"""

from __future__ import annotations

import logging
from typing import Any

import anndata as ad
import scanpy as sc

logger = logging.getLogger(__name__)


def _level_cfg(cfg: dict[str, Any], level: str) -> dict[str, Any]:
    cl = cfg["clustering"]
    if level in cl:
        return cl[level]
    if level == "l1":
        # Backward-compatible: treat the flat block as L1 settings.
        return {k: v for k, v in cl.items() if k not in {"l1", "l2", "use_rep"}}
    return {}


def _slice_rep(adata: ad.AnnData, use_rep: str, n_pcs: int | None) -> str:
    """If ``n_pcs`` is set, materialize a sliced rep so neighbors honors it."""
    if n_pcs is None:
        return use_rep
    X = adata.obsm[use_rep]
    if X.shape[1] <= n_pcs:
        return use_rep
    sliced_key = f"_{use_rep}_top{n_pcs}"
    adata.obsm[sliced_key] = X[:, :n_pcs]
    return sliced_key


def run_level(
    adata: ad.AnnData,
    cfg: dict[str, Any],
    level: str = "l1",
    use_rep: str | None = None,
    do_umap: bool | None = None,
) -> ad.AnnData:
    """Build neighbors / (optional) UMAP / Leiden for a single level."""
    cl_cfg = _level_cfg(cfg, level)
    use_rep = use_rep or cfg["clustering"].get("use_rep", "X_pca_harmony")
    seed = cfg["project"].get("random_seed", 0)
    method = cfg.get("integration", {}).get("method", "harmony").lower()

    n_pcs = cl_cfg.get("n_pcs")
    n_neighbors = cl_cfg.get("n_neighbors", 15)
    resolution = cl_cfg.get("resolution", 0.8)
    leiden_key = f"leiden_{level}"

    if method != "bbknn":  # bbknn already populated obsp
        rep_key = _slice_rep(adata, use_rep, n_pcs)
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=rep_key, random_state=seed)

    if do_umap is None:
        do_umap = level == "l1"
    if do_umap:
        sc.tl.umap(adata, random_state=seed)

    sc.tl.leiden(
        adata,
        resolution=resolution,
        random_state=seed,
        flavor="igraph",
        n_iterations=2,
        directed=False,
        key_added=leiden_key,
    )
    if level == "l1":
        # Keep `leiden` as an alias so legacy code paths and the
        # marker-scoring annotator continue to work without changes.
        adata.obs["leiden"] = adata.obs[leiden_key]
    logger.info("[%s] Leiden produced %d clusters", level, adata.obs[leiden_key].nunique())
    return adata


def run(adata: ad.AnnData, cfg: dict[str, Any], use_rep: str | None = None) -> ad.AnnData:
    """Backwards-compatible single-level entrypoint (L1 only)."""
    return run_level(adata, cfg, level="l1", use_rep=use_rep, do_umap=True)
