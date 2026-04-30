"""End-to-end orchestration of the Visium HD analysis pipeline."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd

from . import (
    annotation,
    clustering,
    differential,
    hierarchical,
    integration,
    io,
    plotting,
    preprocessing,
    qc,
    segmentation as segmentation_mod,
    sketch as sketch_mod,
    spatial,
)
from .utils import get_samples, load_config, output_paths, set_seed, setup_logging

logger = logging.getLogger(__name__)


def run_qc(cfg: dict[str, Any]) -> ad.AnnData:
    paths = output_paths(cfg)
    samples = get_samples(cfg)
    adata = io.read_all_samples(samples, bin_size=cfg["bin_size"])
    adata, summary = qc.run(adata, cfg)
    plotting.qc_violins(adata, paths["qc"], dpi=cfg["plotting"]["dpi"])
    plotting.qc_summary_table(summary, paths["qc"])
    io.write_h5ad(adata, paths["qc"] / "adata_qc.h5ad")
    return adata


def run_preprocess_integrate(cfg: dict[str, Any], adata: ad.AnnData | None = None) -> ad.AnnData:
    paths = output_paths(cfg)
    if adata is None:
        adata = io.read_h5ad(paths["qc"] / "adata_qc.h5ad")
    adata = preprocessing.run(adata, cfg)
    use_rep = integration.run(adata, cfg)
    adata.uns["use_rep"] = use_rep
    io.write_h5ad(adata, paths["integration"] / "adata_integrated.h5ad")
    return adata


def _cluster_annotate_with_sketch(adata: ad.AnnData, cfg: dict[str, Any]) -> ad.AnnData:
    """L1 clustering + annotation done on a leverage-score sketch, then propagated."""
    sk_cfg = cfg["sketch"]
    use_rep = adata.uns.get("use_rep", cfg["clustering"].get("use_rep", "X_pca_harmony"))
    seed = cfg["project"].get("random_seed", 0)
    n_pcs_l1 = cfg["clustering"].get("l1", {}).get("n_pcs")
    integ_method = cfg.get("integration", {}).get("method", "harmony").lower()

    idx = sketch_mod.sketch(
        adata,
        fraction=sk_cfg.get("fraction", 0.15),
        method=sk_cfg.get("method", "leverage_score"),
        n_pcs=n_pcs_l1,
        seed=seed,
    )
    sub = adata[idx].copy()

    # BBKNN builds the kNN graph directly (no corrected embedding). The
    # graph from the full-data run is dropped on slicing, so rebuild it
    # on the sketch; kNN-classifier propagation later runs in X_pca space.
    skip_neighbors = False
    if integ_method == "bbknn":
        integration.run_bbknn(
            sub,
            batch_key=cfg.get("integration", {}).get("batch_key", "sample_id"),
            n_pcs=cfg["preprocessing"].get("n_pcs", 50),
        )
        skip_neighbors = True

    clustering.run_level(sub, cfg, level="l1", use_rep=use_rep, do_umap=True,
                         skip_neighbors=skip_neighbors)
    annotation.run(sub, cfg)

    sketch_mod.propagate_labels(
        adata, idx, sub.obs["leiden_l1"].values, "leiden_l1",
        use_rep=use_rep, n_neighbors=sk_cfg.get("k_propagate", 15),
    )
    sketch_mod.propagate_labels(
        adata, idx, sub.obs["cell_type"].astype(str).values, "cell_type_l1",
        use_rep=use_rep, n_neighbors=sk_cfg.get("k_propagate", 15),
    )
    # Keep aliases pointing at L1 for backwards compatibility.
    adata.obs["leiden"] = adata.obs["leiden_l1"]
    adata.obs["cell_type"] = adata.obs["cell_type_l1"]

    # Carry the sketch UMAP through so it can still be plotted.
    if "X_umap" in sub.obsm:
        import numpy as np
        full_umap = np.full((adata.n_obs, sub.obsm["X_umap"].shape[1]), np.nan)
        full_umap[idx] = sub.obsm["X_umap"]
        adata.obsm["X_umap_sketch"] = full_umap
    if "cluster_marker_scores" in sub.uns:
        adata.uns["l1_marker_scores"] = sub.uns["cluster_marker_scores"]
    return adata


def run_cluster_annotate(cfg: dict[str, Any], adata: ad.AnnData | None = None) -> ad.AnnData:
    paths = output_paths(cfg)
    if adata is None:
        adata = io.read_h5ad(paths["integration"] / "adata_integrated.h5ad")

    if cfg.get("sketch", {}).get("enabled", False):
        adata = _cluster_annotate_with_sketch(adata, cfg)
    else:
        method = cfg.get("integration", {}).get("method", "harmony").lower()
        adata = clustering.run_level(
            adata, cfg, level="l1",
            use_rep=adata.uns.get("use_rep"),
            do_umap=True,
            skip_neighbors=(method == "bbknn"),
        )
        adata = annotation.run(adata, cfg)
        adata.obs["cell_type_l1"] = adata.obs["cell_type"]

    if cfg.get("annotation", {}).get("hierarchical", False):
        adata = hierarchical.run_l2(adata, cfg)
        # Promote L2 to be the primary cell_type label downstream.
        adata.obs["cell_type"] = adata.obs["cell_type_l2"]
        if "l2_marker_scores" in adata.uns:
            for l1, df in adata.uns["l2_marker_scores"].items():
                df.to_csv(paths["annotation"] / f"l2_marker_scores_{l1}.csv", index=False)

    plotting.umap_overview(adata, paths["clusters"], dpi=cfg["plotting"]["dpi"])
    plotting.spatial_per_sample(
        adata, color="cell_type_l1", out_dir=paths["annotation"],
        spot_size=cfg["plotting"]["spot_size"], dpi=cfg["plotting"]["dpi"],
    )
    if "cell_type_l2" in adata.obs:
        plotting.spatial_per_sample(
            adata, color="cell_type_l2", out_dir=paths["annotation"],
            spot_size=cfg["plotting"]["spot_size"], dpi=cfg["plotting"]["dpi"],
        )
    plotting.spatial_per_sample(
        adata, color="leiden_l1", out_dir=paths["clusters"],
        spot_size=cfg["plotting"]["spot_size"], dpi=cfg["plotting"]["dpi"],
    )
    io.write_h5ad(adata, paths["annotation"] / "adata_annotated.h5ad")
    return adata


def run_spatial(cfg: dict[str, Any], adata: ad.AnnData | None = None) -> ad.AnnData:
    paths = output_paths(cfg)
    if adata is None:
        adata = io.read_h5ad(paths["annotation"] / "adata_annotated.h5ad")
    results = spatial.run(adata, cfg)
    for name, df in results.items():
        df.to_csv(paths["spatial"] / f"{name}.csv")
    if "neighborhood_enrichment" in results:
        plotting.neighborhood_heatmap(
            results["neighborhood_enrichment"], paths["spatial"], dpi=cfg["plotting"]["dpi"]
        )
    io.write_h5ad(adata, paths["spatial"] / "adata_spatial.h5ad")
    return adata


def run_differential(cfg: dict[str, Any], adata: ad.AnnData | None = None) -> ad.AnnData:
    paths = output_paths(cfg)
    if adata is None:
        adata = io.read_h5ad(paths["spatial"] / "adata_spatial.h5ad")
    results = differential.run(adata, cfg)
    for name, df in results.items():
        df.to_csv(paths["de"] / f"{name}.csv", index=False)
    if "condition_de_binlevel" in results:
        plotting.save_top_de_heatmap(
            adata,
            results["condition_de_binlevel"],
            paths["de"],
            groupby="cell_type" if "cell_type" in adata.obs else "leiden",
            dpi=cfg["plotting"]["dpi"],
        )
    io.write_h5ad(adata, paths["de"] / "adata_final.h5ad")
    return adata


def run_segmentation(cfg: dict[str, Any]) -> ad.AnnData:
    """Optional: run bin2cell nuclear segmentation on 2 µm bins per sample."""
    return segmentation_mod.run(cfg)


def export_markers_for_review(cfg: dict[str, Any], adata: ad.AnnData | None = None) -> Path:
    """Dump per-cluster marker tables and a cluster_labels.yaml stub for
    manual annotation. Run after stage 3 (clustering) and before
    re-running with annotation.method=manual."""
    paths = output_paths(cfg)
    if adata is None:
        adata = io.read_h5ad(paths["annotation"] / "adata_annotated.h5ad")
    cluster_key = "leiden_l1" if "leiden_l1" in adata.obs else "leiden"
    return annotation.export_marker_tables(
        adata, out_dir=paths["annotation"], cluster_key=cluster_key
    )


def run_all(config_path: str | Path) -> ad.AnnData:
    cfg = load_config(config_path)
    setup_logging()
    set_seed(cfg["project"].get("random_seed", 0))
    adata = run_qc(cfg)
    adata = run_preprocess_integrate(cfg, adata)
    adata = run_cluster_annotate(cfg, adata)
    adata = run_spatial(cfg, adata)
    adata = run_differential(cfg, adata)
    return adata
