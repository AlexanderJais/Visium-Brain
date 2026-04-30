"""Neighborhood graph, UMAP, and Leiden clustering."""

from __future__ import annotations

import logging
from typing import Any

import anndata as ad
import scanpy as sc

logger = logging.getLogger(__name__)


def run(adata: ad.AnnData, cfg: dict[str, Any], use_rep: str | None = None) -> ad.AnnData:
    cl = cfg["clustering"]
    use_rep = use_rep or cl.get("use_rep", "X_pca_harmony")
    seed = cfg["project"].get("random_seed", 0)
    method = cfg.get("integration", {}).get("method", "harmony").lower()

    # BBKNN already populated the graph; skip neighbors.
    if method != "bbknn":
        sc.pp.neighbors(
            adata,
            n_neighbors=cl.get("n_neighbors", 15),
            use_rep=use_rep,
            random_state=seed,
        )
    sc.tl.umap(adata, random_state=seed)
    sc.tl.leiden(
        adata,
        resolution=cl.get("resolution", 0.8),
        random_state=seed,
        flavor="igraph",
        n_iterations=2,
        directed=False,
        key_added="leiden",
    )
    logger.info("Leiden produced %d clusters", adata.obs["leiden"].nunique())
    return adata
