"""Plotting helpers. Each function writes to ``out_dir`` and returns the path."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns

logger = logging.getLogger(__name__)


def qc_violins(adata: ad.AnnData, out_dir: Path, dpi: int = 200) -> Path:
    # Count metrics on log y, percentage metrics on linear y. Visium HD
    # bin total_counts / n_genes_by_counts are heavy-tailed; on linear
    # scale the bulk of typical bins is squashed into a flat slab at the
    # bottom by a handful of high-count outliers, and inter-sample
    # differences disappear. Percentages stay linear so 0-100% reads
    # naturally.
    metrics = ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_hb"]
    log_metrics = {"total_counts", "n_genes_by_counts"}
    fig, axes = plt.subplots(1, len(metrics), figsize=(4 * len(metrics), 4))
    for ax, key in zip(axes, metrics):
        sc.pl.violin(adata, keys=key, groupby="sample_id", rotation=45, ax=ax, show=False)
        if key in log_metrics:
            ax.set_yscale("log")
    fig.tight_layout()
    path = out_dir / "qc_violins.png"
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path


def qc_summary_table(summary, out_dir: Path) -> Path:
    path = out_dir / "qc_summary.csv"
    summary.to_csv(path, index=False)
    return path


def umap_overview(adata: ad.AnnData, out_dir: Path, dpi: int = 200) -> Path:
    keys = [k for k in ["leiden", "cell_type", "condition", "sample_id"] if k in adata.obs]
    sc.pl.umap(adata, color=keys, ncols=2, show=False)
    fig = plt.gcf()
    # When the sketch path was used, X_umap is the sketch UMAP NaN-padded
    # to the full bin set (matplotlib silently drops NaN coords). Stamp
    # the title and the filename so a reader cannot mistake the
    # apparently-sparse plot for the actual bin density, and so the
    # filename distinguishes a sketch UMAP from a full-data UMAP at a
    # glance.
    if "X_umap_sketch" in adata.obsm and "X_umap" in adata.obsm:
        n_shown = int(np.isfinite(np.asarray(adata.obsm["X_umap"])[:, 0]).sum())
        fig.suptitle(
            f"sketch UMAP — {n_shown:,}/{adata.n_obs:,} bins shown",
            fontsize=10,
        )
        path = out_dir / "umap_overview_sketch.png"
    else:
        path = out_dir / "umap_overview.png"
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path


def spatial_per_sample(
    adata: ad.AnnData,
    color: str,
    out_dir: Path,
    spot_size: float = 1.4,
    dpi: int = 200,
) -> list[Path]:
    paths = []
    for sid in adata.obs["sample_id"].cat.categories:
        sub = adata[adata.obs["sample_id"] == sid].copy()
        sub.uns["spatial"] = {sid: adata.uns["spatial"][sid]} if sid in adata.uns.get("spatial", {}) else {}
        sc.pl.spatial(
            sub,
            color=color,
            library_id=sid if sid in sub.uns.get("spatial", {}) else None,
            spot_size=spot_size,
            show=False,
        )
        path = out_dir / f"spatial_{color}_{sid}.png"
        plt.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close()
        paths.append(path)
    return paths


def neighborhood_heatmap(matrix, out_dir: Path, dpi: int = 200, name: str = "nhood_enrichment") -> Path:
    # Scale figure size and tick-label font to category count: with
    # hierarchical L1->L2 the matrix can be 30-50 cell types per side,
    # at which point the previous fixed (8, 7) figure ran tick labels
    # into each other.
    n = matrix.shape[0]
    side = max(7.0, 0.35 * n + 2.0)
    fontsize = max(6, min(10, int(round(200.0 / max(n, 1)))))
    fig, ax = plt.subplots(figsize=(side + 1.0, side))  # extra width for colorbar
    sns.heatmap(matrix, cmap="vlag", center=0, ax=ax, square=True, cbar_kws={"label": "z-score"})
    ax.set_title("Neighborhood enrichment")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=fontsize)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=fontsize)
    path = out_dir / f"{name}.png"
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path


def save_top_de_heatmap(
    adata: ad.AnnData,
    de_df,
    out_dir: Path,
    n_top: int = 5,
    groupby: str = "cell_type",
    dpi: int = 200,
    max_genes: int = 60,
) -> Path | None:
    """Dotplot of the top condition-DE genes per cluster.

    ``max_genes`` caps the number of distinct genes shown. With
    hierarchical L1->L2 annotation the de_df can carry 30-50 L2
    clusters; even at ``n_top=5`` that produces 150+ unique genes and
    the dotplot becomes an unreadable barcode. Genes are taken in
    cluster order (each cluster contributes its top ``n_top`` ranked by
    pval_adj asc / score desc, deduplicated against earlier clusters'
    picks), then the head ``max_genes`` is plotted.
    """
    if de_df.empty:
        return None
    ranked = de_df.sort_values(
        ["cluster", "pval_adj", "score"], ascending=[True, True, False]
    )
    top = ranked.groupby("cluster").head(n_top)["gene"]
    # Deduplicate while preserving rank order (every cluster gets at
    # least its #1 marker before any cluster gets its #2, then so on).
    top = top.drop_duplicates().tolist()
    top = [g for g in top if g in adata.var_names][:max_genes]
    if not top:
        return None
    # use_raw=True so the dotplot reflects log1p expression, not the
    # scaled values left in .X by sc.pp.scale.
    sc.pl.dotplot(
        adata, var_names=top, groupby=groupby, standard_scale="var",
        use_raw=True, show=False,
    )
    path = out_dir / f"de_top{n_top}_dotplot.png"
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close()
    return path
