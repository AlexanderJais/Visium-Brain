"""Plotting helpers. Each function writes to ``out_dir`` and returns the path."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable

import anndata as ad
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns

logger = logging.getLogger(__name__)


def qc_violins(adata: ad.AnnData, out_dir: Path, dpi: int = 200) -> Path:
    metrics = ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_hb"]
    fig = sc.pl.violin(
        adata, keys=metrics, groupby="sample_id", rotation=45, multi_panel=True, show=False
    )
    path = out_dir / "qc_violins.png"
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close()
    return path


def qc_summary_table(summary, out_dir: Path) -> Path:
    path = out_dir / "qc_summary.csv"
    summary.to_csv(path, index=False)
    return path


def umap_overview(adata: ad.AnnData, out_dir: Path, dpi: int = 200) -> Path:
    keys = [k for k in ["leiden", "cell_type", "condition", "sample_id"] if k in adata.obs]
    sc.pl.umap(adata, color=keys, ncols=2, show=False)
    path = out_dir / "umap_overview.png"
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close()
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
    fig, ax = plt.subplots(figsize=(8, 7))
    sns.heatmap(matrix, cmap="vlag", center=0, ax=ax, square=True, cbar_kws={"label": "z-score"})
    ax.set_title("Neighborhood enrichment")
    path = out_dir / f"{name}.png"
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close()
    return path


def save_top_de_heatmap(
    adata: ad.AnnData,
    de_df,
    out_dir: Path,
    n_top: int = 5,
    groupby: str = "cell_type",
    dpi: int = 200,
) -> Path | None:
    if de_df.empty:
        return None
    top = (
        de_df.sort_values(["cluster", "pval_adj", "score"], ascending=[True, True, False])
        .groupby("cluster")
        .head(n_top)["gene"].unique().tolist()
    )
    top = [g for g in top if g in adata.var_names]
    if not top:
        return None
    sc.pl.dotplot(adata, var_names=top, groupby=groupby, standard_scale="var", show=False)
    path = out_dir / f"de_top{n_top}_dotplot.png"
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close()
    return path
