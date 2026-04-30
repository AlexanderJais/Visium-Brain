"""Normalization, HVG selection, scaling, and PCA.

Run order (`preprocessing.run`):

1. Stash raw counts in ``adata.layers['counts']`` (kept for scVI,
   pseudobulk DE, and ``seurat_v3`` HVG selection).
2. ``sc.pp.normalize_total`` to ``target_sum`` and ``log1p``.
3. **Snapshot log1p into ``adata.raw``** so that downstream tools
   (``sc.tl.score_genes``, ``sc.tl.rank_genes_groups``,
   ``sc.pl.dotplot``) can read log-normalized expression via
   ``use_raw=True``. Without this, the scaling step in (5) would
   overwrite ``.X`` and DE log-fold-changes would be computed on
   z-scores instead of log1p values.
4. Highly-variable gene selection. Default flavor ``seurat_v3`` runs on
   the saved raw counts and is batch-aware via ``integration.batch_key``.
5. ``sc.pp.scale`` on HVGs only, with ``max_value`` clipping.
6. ``sc.tl.pca`` on HVGs producing ``obsm['X_pca']`` and
   ``uns['pca']['variance']`` (consumed by ``sketch``).
"""

from __future__ import annotations

import logging
from typing import Any

import anndata as ad
import scanpy as sc

logger = logging.getLogger(__name__)


def normalize_log(adata: ad.AnnData, target_sum: float = 1e4, log1p: bool = True) -> None:
    # Preserve raw counts for downstream tools (scVI, pseudobulk DE, seurat_v3 HVGs).
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    if log1p:
        sc.pp.log1p(adata)


def snapshot_log1p_raw(adata: ad.AnnData) -> None:
    """Snapshot the log1p ``.X`` into ``adata.raw`` so subsequent scaling
    (which overwrites ``.X``) does not corrupt downstream marker scoring
    and DE. After this, every ``sc.tl.score_genes`` /
    ``sc.tl.rank_genes_groups`` / ``sc.pl.dotplot`` call should pass
    ``use_raw=True`` to read log-normalized values.
    """
    adata.raw = adata


def select_hvgs(
    adata: ad.AnnData,
    n_top: int = 4000,
    flavor: str = "seurat_v3",
    batch_key: str | None = "sample_id",
) -> None:
    if flavor == "seurat_v3":
        # seurat_v3 expects raw counts in .X; use the saved layer.
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top,
            flavor="seurat_v3",
            batch_key=batch_key,
            layer="counts",
        )
    else:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top,
            flavor=flavor,
            batch_key=batch_key,
        )
    logger.info("Selected %d HVGs (flavor=%s)", int(adata.var["highly_variable"].sum()), flavor)


def scale_and_pca(adata: ad.AnnData, max_value: float = 10.0, n_pcs: int = 50, seed: int = 0) -> None:
    sc.pp.scale(adata, max_value=max_value, mask_var="highly_variable")
    sc.tl.pca(adata, n_comps=n_pcs, mask_var="highly_variable", random_state=seed)


def run(adata: ad.AnnData, cfg: dict[str, Any]) -> ad.AnnData:
    pp = cfg["preprocessing"]
    integ = cfg.get("integration", {})
    normalize_log(adata, target_sum=pp.get("target_sum", 1e4), log1p=pp.get("log1p", True))
    # Snapshot log1p before HVG/scale so downstream DE can read it via use_raw=True.
    snapshot_log1p_raw(adata)
    select_hvgs(
        adata,
        n_top=pp.get("n_top_hvgs", 4000),
        flavor=pp.get("hvg_flavor", "seurat_v3"),
        batch_key=integ.get("batch_key", "sample_id"),
    )
    scale_and_pca(
        adata,
        max_value=pp.get("scale_max_value", 10.0),
        n_pcs=pp.get("n_pcs", 50),
        seed=cfg["project"].get("random_seed", 0),
    )
    return adata
