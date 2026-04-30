"""Two-level (L1 -> L2) hierarchical clustering and annotation.

Mirrors the workflow from Oh et al., *Nat Genet* 2025
(https://www.nature.com/articles/s41588-025-02193-3): cluster + annotate
broad cell types at L1, then for each L1 label recluster at lower
resolution to refine subtypes (L2) and annotate.

The annotator is the same marker-score panel used at L1
(``annotation.MOUSE_BRAIN_MARKERS``); within each L1 subset every panel
is rescored and clusters are labelled by their best-scoring panel,
prefixed by the parent L1 label so the hierarchy stays visible.

Per-cluster mean marker-score tables are stashed in
``adata.uns['l2_marker_scores'][<L1 label>]`` to support manual review.
"""

from __future__ import annotations

import logging
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from . import annotation

logger = logging.getLogger(__name__)


def _l2_settings(cfg: dict[str, Any]) -> dict[str, Any]:
    cl = cfg["clustering"]
    return cl.get("l2", {"n_pcs": 25, "n_neighbors": 15, "resolution": 0.1})


def _parent_label_key(adata: ad.AnnData) -> str:
    for k in ("cell_type_l1", "cell_type", "leiden_l1", "leiden"):
        if k in adata.obs:
            return k
    raise KeyError("No L1 label found in adata.obs (expected one of cell_type_l1/cell_type/leiden_l1/leiden).")


def run_l2(adata: ad.AnnData, cfg: dict[str, Any]) -> ad.AnnData:
    """Recluster within each L1 label and write ``cell_type_l2``."""
    seed = cfg["project"].get("random_seed", 0)
    l2 = _l2_settings(cfg)
    use_rep = adata.uns.get("use_rep") or cfg["clustering"].get("use_rep", "X_pca_harmony")
    if use_rep not in adata.obsm:
        raise KeyError(f"adata.obsm[{use_rep!r}] missing; run preprocessing/integration first.")
    n_pcs = l2.get("n_pcs", 25)
    n_neighbors = l2.get("n_neighbors", 15)
    resolution = l2.get("resolution", 0.1)
    min_cells = cfg.get("annotation", {}).get("l2_min_cells", 200)

    parent_key = _parent_label_key(adata)
    if not isinstance(adata.obs[parent_key].dtype, pd.CategoricalDtype):
        adata.obs[parent_key] = adata.obs[parent_key].astype("category")

    l2_labels = pd.Series(index=adata.obs_names, dtype=object)
    leiden_l2 = pd.Series(index=adata.obs_names, dtype=object)
    score_tables: dict[str, pd.DataFrame] = {}

    for l1 in adata.obs[parent_key].cat.categories:
        mask = (adata.obs[parent_key] == l1).values
        n = int(mask.sum())
        if n < min_cells:
            l2_labels.iloc[np.where(mask)[0]] = f"{l1}__c0"
            leiden_l2.iloc[np.where(mask)[0]] = "0"
            logger.info("[L2] %s: %d bins < min_cells=%d, leaving as single subgroup",
                        l1, n, min_cells)
            continue

        sub = adata[mask].copy()
        rep = sub.obsm[use_rep]
        if rep.shape[1] > n_pcs:
            rep = rep[:, :n_pcs]
        sub.obsm["_l2rep"] = rep
        sc.pp.neighbors(sub, n_neighbors=min(n_neighbors, sub.n_obs - 1),
                        use_rep="_l2rep", random_state=seed)
        sc.tl.leiden(sub, resolution=resolution, random_state=seed, flavor="igraph",
                     n_iterations=2, directed=False, key_added="leiden_l2")

        score_cols = annotation.score_markers(sub)
        if score_cols:
            annotation.assign_cluster_labels(
                sub, score_cols=score_cols, cluster_key="leiden_l2",
                out_key="cell_type_l2_local",
            )
            sub_label = sub.obs["cell_type_l2_local"].astype(str)
            score_tables[str(l1)] = sub.uns["cluster_marker_scores"]
        else:
            sub_label = pd.Series([f"sub{i}" for i in sub.obs["leiden_l2"]], index=sub.obs_names)

        composed = (str(l1) + "__" + sub_label + "_c" + sub.obs["leiden_l2"].astype(str)).values
        idx = np.where(mask)[0]
        l2_labels.iloc[idx] = composed
        leiden_l2.iloc[idx] = sub.obs["leiden_l2"].astype(str).values
        logger.info("[L2] %s: %d bins -> %d sub-clusters",
                    l1, n, sub.obs["leiden_l2"].nunique())

    adata.obs["leiden_l2"] = pd.Categorical(leiden_l2.fillna("0").astype(str))
    adata.obs["cell_type_l2"] = pd.Categorical(l2_labels.fillna("Unknown").astype(str))
    if score_tables:
        adata.uns["l2_marker_scores"] = score_tables
    return adata
