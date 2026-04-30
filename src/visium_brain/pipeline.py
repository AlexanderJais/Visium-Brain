"""End-to-end orchestration of the Visium HD analysis pipeline."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import anndata as ad

from . import annotation, clustering, differential, integration, io, plotting, preprocessing, qc, spatial
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


def run_cluster_annotate(cfg: dict[str, Any], adata: ad.AnnData | None = None) -> ad.AnnData:
    paths = output_paths(cfg)
    if adata is None:
        adata = io.read_h5ad(paths["integration"] / "adata_integrated.h5ad")
    adata = clustering.run(adata, cfg, use_rep=adata.uns.get("use_rep"))
    adata = annotation.run(adata, cfg)
    plotting.umap_overview(adata, paths["clusters"], dpi=cfg["plotting"]["dpi"])
    plotting.spatial_per_sample(
        adata, color="cell_type", out_dir=paths["annotation"],
        spot_size=cfg["plotting"]["spot_size"], dpi=cfg["plotting"]["dpi"],
    )
    plotting.spatial_per_sample(
        adata, color="leiden", out_dir=paths["clusters"],
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
