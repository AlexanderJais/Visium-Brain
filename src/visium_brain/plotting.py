"""Plotting helpers. Each function writes to ``out_dir`` and returns the path."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

logger = logging.getLogger(__name__)


# Cardinality threshold above which we pre-bake a HUSL-spaced palette
# into adata.uns. Below this, scanpy's default_20 / default_28 already
# pick visually distinct colors and we leave the choice to scanpy.
_LARGE_PALETTE_THRESHOLD = 20


def _ensure_categorical_palette(adata: ad.AnnData, key: str) -> None:
    """Ensure ``adata.uns[f'{key}_colors']`` is set to a perceptually
    uniform palette when the category count would otherwise rely on
    scanpy's default_102.

    With hierarchical L1->L2 annotation, ``cell_type_l2`` has 30-50
    alphabetically-sorted categories like ``Excitatory_neuron__Cux2_c0``
    / ``...__Cux2_c1``. scanpy's default_102 is curated for max
    distinguishability across the *whole* set of 102 colors, but
    adjacent positions can still look similar -- and when alphabetical
    sort puts related subtypes next to each other in the category
    order, those visually-similar colors land on biologically-similar
    clusters. A HUSL palette spreads hues evenly across N regardless
    of order, so adjacent categories are always maximally apart in hue.

    Idempotent: if the colors slot is already populated and matches the
    cardinality, leave it alone (so a re-plot reuses the same colors).
    """
    if key not in adata.obs:
        return
    col = adata.obs[key]
    if not hasattr(col, "cat"):
        return
    n = len(col.cat.categories)
    if n <= _LARGE_PALETTE_THRESHOLD:
        return
    palette_key = f"{key}_colors"
    existing = adata.uns.get(palette_key)
    if existing is not None and len(existing) == n:
        return
    palette = sns.color_palette("husl", n_colors=n).as_hex()
    adata.uns[palette_key] = list(palette)


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
    summary.to_csv(path, index=False, float_format="%.4g")
    return path


def umap_overview(adata: ad.AnnData, out_dir: Path, dpi: int = 200) -> Path:
    keys = [k for k in ["leiden", "cell_type", "condition", "sample_id"] if k in adata.obs]
    for k in keys:
        _ensure_categorical_palette(adata, k)
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
    # Pre-bake a HUSL palette on the full adata once. Each per-sample
    # sub.copy() inherits adata.uns, so every panel uses the SAME
    # category-to-color mapping -- the same L2 cell type is plotted in
    # the same color across ctrl_m1, ctrl_m2, trt_m3, trt_m4.
    _ensure_categorical_palette(adata, color)
    # Map sample_id -> condition once. Filenames include the condition
    # so figures stay self-describing even when sample_ids are generic
    # (S1, S2, ...) and the user later loses track of which is which.
    if "condition" in adata.obs:
        # dropna here so NaN condition rows (rare, but possible if the
        # user's config is incomplete) don't end up as the literal
        # string "nan" in the filename when astype(str) coerces them.
        cond_series = (
            adata.obs.drop_duplicates("sample_id")
            .set_index("sample_id")["condition"]
            .dropna()
        )
        cond_map = cond_series.astype(str).to_dict()
    else:
        cond_map = {}
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
        cond = cond_map.get(sid, "")
        prefix = f"spatial_{color}_{cond}_{sid}" if cond else f"spatial_{color}_{sid}"
        path = out_dir / f"{prefix}.png"
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


def save_pseudobulk_volcano(
    de_df,
    out_dir: Path,
    max_clusters: int = 12,
    n_label: int = 8,
    dpi: int = 200,
) -> Path | None:
    """Per-cluster volcano-style figure for the pseudobulk DE table.

    Designed to be useful in both regimes:

    * **Powered** (``low_power=False``, p-values populated): standard
      volcano with ``-log10(pval_adj)`` on y.
    * **Underpowered** (``low_power=True``, p-values NaN'd by
      ``pseudobulk_de``): switches y to ``|Wilcoxon U|`` so the panel
      still has a meaningful "separation" axis. The panel title is
      tagged ``(low power)`` so a reader cannot mistake the alternate
      y-axis for a p-value-based volcano.

    Top ``n_label`` genes per cluster (by ``|logfoldchange|``) are
    labelled; cap at ``max_clusters`` panels so the figure stays
    readable for hierarchical L2 (30-50 clusters).
    """
    if de_df.empty or "logfoldchange" not in de_df.columns:
        return None
    # Preserve original cluster order without sorting (dict.fromkeys
    # is order-preserving in Python 3.7+).
    clusters = list(dict.fromkeys(de_df["cluster"].tolist()))[:max_clusters]
    if not clusters:
        return None
    ncols = min(3, len(clusters))
    nrows = (len(clusters) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.5 * nrows), squeeze=False)

    for ax, c in zip(axes.flat, clusters):
        sub = de_df[de_df["cluster"] == c]
        if sub.empty:
            ax.set_visible(False)
            continue
        # Guard against NaN low_power: bool(NaN) is True, which would
        # silently mislabel the panel. pseudobulk_de always populates
        # this with a real bool, but be defensive against hand-built
        # frames or future schema drift.
        low_power_raw = sub["low_power"].iloc[0] if "low_power" in sub.columns else False
        low_power = bool(low_power_raw) if pd.notna(low_power_raw) else False

        x = sub["logfoldchange"].to_numpy()
        if not low_power and "pval_adj" in sub.columns and sub["pval_adj"].notna().any():
            y = -np.log10(sub["pval_adj"].clip(lower=1e-300).to_numpy())
            ylabel = "-log10 adj. p"
        else:
            y = sub["score"].abs().to_numpy()
            ylabel = "|Wilcoxon U|"

        ax.scatter(x, y, s=8, alpha=0.6, c="steelblue", edgecolors="none")
        ax.axvline(0, color="k", lw=0.5, alpha=0.3)
        title = str(c) + (" (low power)" if low_power else "")
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("log2 fold-change")
        ax.set_ylabel(ylabel)

        # Label top |lfc| hits.
        top_idx = sub["logfoldchange"].abs().nlargest(n_label).index
        y_series = pd.Series(y, index=sub.index)
        for i in top_idx:
            ax.annotate(
                str(sub.loc[i, "gene"]),
                (float(sub.loc[i, "logfoldchange"]), float(y_series.loc[i])),
                fontsize=6, alpha=0.8,
            )

    for ax in axes.flat[len(clusters):]:
        ax.set_visible(False)
    fig.tight_layout()
    path = out_dir / "pseudobulk_volcano.png"
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
    # ranked is sorted by cluster first, so head(n_top) emits cluster A's
    # top n_top, then cluster B's, and so on. drop_duplicates keeps each
    # gene at its first appearance, so a gene that ranks well in
    # multiple clusters is attributed to the alphabetically-first one
    # and later clusters effectively contribute fewer unique columns.
    # Acceptable here -- the dotplot is a summary, not a per-cluster
    # exhaustive marker table (those live in cluster_markers*.csv).
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
