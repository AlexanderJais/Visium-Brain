"""Differential expression: per cluster, per condition, and pseudobulk.

For Visium HD with biological replicates (2 mice per condition here), the
recommended approach is pseudobulk DE: sum raw counts across bins per
(sample, cluster), then test condition effects with a negative-binomial /
Wilcoxon model over mice. We provide both the bin-level wilcoxon (fast,
exploratory) and the pseudobulk version (statistically valid).
"""

from __future__ import annotations

import logging
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

logger = logging.getLogger(__name__)


def cluster_markers(adata: ad.AnnData, groupby: str = "leiden", method: str = "wilcoxon") -> pd.DataFrame:
    """Find markers for each cluster vs the rest."""
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method, use_raw=False)
    return _rank_genes_to_df(adata)


def condition_de(
    adata: ad.AnnData,
    groupby: str = "condition",
    reference: str = "control",
    method: str = "wilcoxon",
    cluster_key: str | None = "cell_type",
) -> pd.DataFrame:
    """Per-cluster DE between conditions (bin-level; fast exploratory view)."""
    frames = []
    if cluster_key is None or cluster_key not in adata.obs:
        sc.tl.rank_genes_groups(adata, groupby=groupby, reference=reference, method=method, use_raw=False)
        df = _rank_genes_to_df(adata)
        df["cluster"] = "all"
        return df

    for c in adata.obs[cluster_key].cat.categories:
        sub = adata[adata.obs[cluster_key] == c]
        if sub.obs[groupby].nunique() < 2:
            continue
        sc.tl.rank_genes_groups(sub, groupby=groupby, reference=reference, method=method, use_raw=False)
        df = _rank_genes_to_df(sub)
        df["cluster"] = str(c)
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _rank_genes_to_df(adata: ad.AnnData) -> pd.DataFrame:
    res = adata.uns["rank_genes_groups"]
    groups = res["names"].dtype.names
    rows = []
    for g in groups:
        rows.append(
            pd.DataFrame(
                {
                    "group": g,
                    "gene": res["names"][g],
                    "logfoldchange": res["logfoldchanges"][g],
                    "pval": res["pvals"][g],
                    "pval_adj": res["pvals_adj"][g],
                    "score": res["scores"][g],
                }
            )
        )
    return pd.concat(rows, ignore_index=True)


def make_pseudobulk(
    adata: ad.AnnData,
    sample_key: str = "sample_id",
    cluster_key: str = "cell_type",
    condition_key: str = "condition",
    counts_layer: str = "counts",
    min_cells: int = 25,
) -> ad.AnnData:
    """Sum raw counts per (sample, cluster) into a pseudobulk AnnData."""
    if counts_layer not in adata.layers:
        raise KeyError(
            f"Expected raw counts in layer '{counts_layer}'. Run preprocessing.normalize_log first."
        )
    obs = adata.obs[[sample_key, cluster_key, condition_key]].copy()
    obs["_pb"] = obs[sample_key].astype(str) + "::" + obs[cluster_key].astype(str)
    groups = obs.groupby("_pb", observed=True)

    rows_X, rows_obs = [], []
    X = adata.layers[counts_layer]
    is_sparse = sp.issparse(X)
    for key, idx in groups.indices.items():
        if len(idx) < min_cells:
            continue
        block = X[idx]
        summed = np.asarray(block.sum(axis=0)).ravel() if is_sparse else block.sum(axis=0)
        rows_X.append(summed)
        meta = obs.iloc[idx[0]]
        rows_obs.append(
            {
                "pseudobulk_id": key,
                sample_key: meta[sample_key],
                cluster_key: meta[cluster_key],
                condition_key: meta[condition_key],
                "n_bins": len(idx),
            }
        )
    if not rows_X:
        raise ValueError("No pseudobulk groups passed the min_cells threshold.")

    pb = ad.AnnData(
        X=np.vstack(rows_X).astype(np.float32),
        obs=pd.DataFrame(rows_obs).set_index("pseudobulk_id"),
        var=adata.var[[]].copy(),
    )
    return pb


def pseudobulk_de(
    pb: ad.AnnData,
    cluster_key: str = "cell_type",
    condition_key: str = "condition",
    reference: str = "control",
) -> pd.DataFrame:
    """Wilcoxon test on log-CPM pseudobulks per cluster across conditions.

    With 2 mice per condition this is severely underpowered for any single
    test, but it gives effect sizes and ranks that are interpretable across
    the full transcriptome. For publication-grade DE swap in pyDESeq2 or
    edgeR via rpy2.
    """
    # log-CPM normalize the pseudobulks.
    X = pb.X.astype(float)
    libsize = X.sum(axis=1, keepdims=True)
    libsize[libsize == 0] = 1
    logcpm = np.log2(1e6 * X / libsize + 1)
    pb2 = ad.AnnData(X=logcpm, obs=pb.obs.copy(), var=pb.var.copy())

    frames = []
    for c in pb2.obs[cluster_key].unique():
        sub = pb2[pb2.obs[cluster_key] == c]
        if sub.obs[condition_key].nunique() < 2:
            continue
        sc.tl.rank_genes_groups(
            sub, groupby=condition_key, reference=reference, method="wilcoxon", use_raw=False
        )
        df = _rank_genes_to_df(sub)
        df["cluster"] = str(c)
        df["n_pseudobulks"] = sub.n_obs
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def run(adata: ad.AnnData, cfg: dict[str, Any]) -> dict[str, pd.DataFrame]:
    de_cfg = cfg["differential_expression"]
    out: dict[str, pd.DataFrame] = {}

    out["cluster_markers"] = cluster_markers(adata, groupby="leiden", method=de_cfg.get("method", "wilcoxon"))
    out["condition_de_binlevel"] = condition_de(
        adata,
        groupby=de_cfg.get("groupby", "condition"),
        reference=de_cfg.get("reference", "control"),
        method=de_cfg.get("method", "wilcoxon"),
        cluster_key="cell_type" if "cell_type" in adata.obs else "leiden",
    )

    if de_cfg.get("pseudobulk", True):
        cluster_key = "cell_type" if "cell_type" in adata.obs else "leiden"
        pb = make_pseudobulk(
            adata,
            sample_key="sample_id",
            cluster_key=cluster_key,
            condition_key=de_cfg.get("groupby", "condition"),
            min_cells=de_cfg.get("pseudobulk_min_cells", 25),
        )
        out["pseudobulk_obs"] = pb.obs.reset_index()
        out["pseudobulk_de"] = pseudobulk_de(
            pb,
            cluster_key=cluster_key,
            condition_key=de_cfg.get("groupby", "condition"),
            reference=de_cfg.get("reference", "control"),
        )
    return out
