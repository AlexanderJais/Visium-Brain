"""Quality control for Visium HD bins.

Computes standard scanpy QC metrics, flags mitochondrial / hemoglobin /
ribosomal genes, and filters bins/genes per the configured thresholds.
"""

from __future__ import annotations

import logging
import re
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

logger = logging.getLogger(__name__)


def _flag_genes(adata: ad.AnnData, mito_prefix: str, hb_prefix: str) -> None:
    names = adata.var_names.astype(str)
    adata.var["mt"] = names.str.startswith(mito_prefix)
    adata.var["hb"] = names.str.contains(re.compile(hb_prefix))
    adata.var["ribo"] = names.str.startswith(("Rps", "Rpl"))


def compute_qc(adata: ad.AnnData, mito_prefix: str = "mt-", hb_prefix: str = "Hb[ab]-") -> ad.AnnData:
    """Annotate var with gene families and obs with QC metrics."""
    _flag_genes(adata, mito_prefix=mito_prefix, hb_prefix=hb_prefix)
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "hb", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    return adata


def filter_bins_and_genes(
    adata: ad.AnnData,
    min_counts_per_bin: int = 50,
    min_genes_per_bin: int = 20,
    max_pct_mito: float = 25.0,
    min_cells_per_gene: int = 10,
) -> ad.AnnData:
    """Apply per-sample QC thresholds and report counts removed."""
    n0 = adata.n_obs
    sc.pp.filter_cells(adata, min_counts=min_counts_per_bin)
    sc.pp.filter_cells(adata, min_genes=min_genes_per_bin)
    adata = adata[adata.obs["pct_counts_mt"] <= max_pct_mito].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    logger.info(
        "QC: kept %d / %d bins (%.1f%%); %d genes",
        adata.n_obs,
        n0,
        100.0 * adata.n_obs / max(n0, 1),
        adata.n_vars,
    )
    return adata


def per_sample_summary(adata: ad.AnnData) -> pd.DataFrame:
    """Produce a per-sample QC summary table."""
    by = adata.obs.groupby("sample_id", observed=True)
    summary = pd.DataFrame(
        {
            "n_bins": by.size(),
            "median_counts": by["total_counts"].median(),
            "median_genes": by["n_genes_by_counts"].median(),
            "median_pct_mt": by["pct_counts_mt"].median(),
            "median_pct_hb": by["pct_counts_hb"].median(),
            "median_pct_ribo": by["pct_counts_ribo"].median(),
        }
    )
    cond = adata.obs.drop_duplicates("sample_id").set_index("sample_id")["condition"]
    summary["condition"] = cond
    return summary.reset_index()


def run(adata: ad.AnnData, cfg: dict[str, Any]) -> tuple[ad.AnnData, pd.DataFrame]:
    qc_cfg = cfg["qc"]
    adata = compute_qc(
        adata,
        mito_prefix=qc_cfg.get("mito_prefix", "mt-"),
        hb_prefix=qc_cfg.get("hb_prefix", "Hb[ab]-"),
    )
    pre_summary = per_sample_summary(adata)
    adata = filter_bins_and_genes(
        adata,
        min_counts_per_bin=qc_cfg.get("min_counts_per_bin", 50),
        min_genes_per_bin=qc_cfg.get("min_genes_per_bin", 20),
        max_pct_mito=qc_cfg.get("max_pct_mito", 25.0),
        min_cells_per_gene=qc_cfg.get("min_cells_per_gene", 10),
    )
    post_summary = per_sample_summary(adata)
    pre_summary["stage"] = "pre"
    post_summary["stage"] = "post"
    return adata, pd.concat([pre_summary, post_summary], ignore_index=True)
