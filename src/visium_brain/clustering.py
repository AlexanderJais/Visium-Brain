"""Neighborhood graph, UMAP, and Leiden clustering.

Supports per-level configuration via ``cfg['clustering']['l1']`` and
``cfg['clustering']['l2']`` (number of PCs, neighbors, resolution). If a
flat ``cfg['clustering']`` block is provided (legacy), it is treated as
the L1 config.

The neighbors-graph step can be skipped by passing
``skip_neighbors=True`` to :func:`run_level`. The caller (typically
``pipeline.run_cluster_annotate``) is responsible for setting this when
``obsp`` is already populated by an integration step that builds the
graph directly (BBKNN). This was previously inferred from
``cfg['integration']['method']`` inside ``run_level`` itself, which
silently produced an empty graph on the sketch when both BBKNN and the
sketch path were enabled (BBKNN had only ever run on the full adata).
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
    skip_neighbors: bool = False,
) -> ad.AnnData:
    """Build neighbors / (optional) UMAP / Leiden for a single level.

    ``skip_neighbors`` lets the caller declare that ``adata.obsp`` is
    already populated (e.g. by BBKNN). When ``False`` (default) we
    always call ``sc.pp.neighbors``.
    """
    cl_cfg = _level_cfg(cfg, level)
    use_rep = use_rep or cfg["clustering"].get("use_rep", "X_pca_harmony")
    seed = cfg["project"].get("random_seed", 0)

    n_pcs = cl_cfg.get("n_pcs")
    n_neighbors = cl_cfg.get("n_neighbors", 15)
    resolution = cl_cfg.get("resolution", 0.8)
    leiden_key = f"leiden_{level}"

    if not skip_neighbors:
        rep_key = _slice_rep(adata, use_rep, n_pcs)
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=rep_key, random_state=seed)
        # neighbors() has stashed the kNN graph in obsp; the sliced rep
        # is no longer needed and would otherwise persist in every
        # downstream h5ad as a duplicate of obsm[use_rep].
        if rep_key != use_rep:
            del adata.obsm[rep_key]

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
    """Backwards-compatible single-level entrypoint (L1 only).

    Honors the legacy implicit-skip for BBKNN to keep downstream callers
    that don't yet pass ``skip_neighbors`` working. New code should call
    :func:`run_level` directly with an explicit ``skip_neighbors``.
    """
    method = cfg.get("integration", {}).get("method", "harmony").lower()
    return run_level(
        adata, cfg, level="l1", use_rep=use_rep, do_umap=True,
        skip_neighbors=(method == "bbknn"),
    )
