"""Smoke tests against synthetic Visium-HD-shaped AnnData.

These don't require Space Ranger output and exercise the math-y parts
of the pipeline (QC, preprocessing, clustering, pseudobulk DE).
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp


def _make_synthetic_adata(
    n_per_sample: int = 200,
    n_genes: int = 300,
    samples: tuple[tuple[str, str, str], ...] = (
        ("ctrl_m1", "control", "M1"),
        ("ctrl_m2", "control", "M2"),
        ("trt_m3", "treated", "M3"),
        ("trt_m4", "treated", "M4"),
    ),
) -> ad.AnnData:
    rng = np.random.default_rng(0)
    obs_frames = []
    Xs = []
    coords = []
    for sid, cond, mid in samples:
        counts = rng.poisson(0.5, size=(n_per_sample, n_genes)).astype(np.float32)
        # Inject a treatment effect on a handful of genes for treated samples.
        if cond == "treated":
            counts[:, :10] += rng.poisson(2.0, size=(n_per_sample, 10)).astype(np.float32)
        Xs.append(counts)
        xy = rng.uniform(0, 1000, size=(n_per_sample, 2))
        coords.append(xy)
        obs_frames.append(
            pd.DataFrame(
                {
                    "sample_id": sid,
                    "condition": cond,
                    "mouse_id": mid,
                },
                index=[f"{sid}_bc{i}" for i in range(n_per_sample)],
            )
        )
    X = sp.csr_matrix(np.vstack(Xs))
    obs = pd.concat(obs_frames)
    var_names = [f"gene_{i}" for i in range(n_genes)]
    var_names[0] = "mt-Atp6"  # mito for QC
    var_names[1] = "mt-Co1"
    var = pd.DataFrame(index=var_names)
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obsm["spatial"] = np.vstack(coords)
    adata.obs["sample_id"] = adata.obs["sample_id"].astype("category")
    adata.obs["condition"] = adata.obs["condition"].astype("category")
    adata.obs["mouse_id"] = adata.obs["mouse_id"].astype("category")
    return adata


def test_qc_filters_and_summary():
    from visium_brain import qc

    adata = _make_synthetic_adata()
    cfg = {
        "qc": {
            "mito_prefix": "mt-",
            "hb_prefix": "Hb[ab]-",
            "min_counts_per_bin": 1,
            "min_genes_per_bin": 1,
            "max_pct_mito": 100.0,
            "min_cells_per_gene": 1,
        }
    }
    out, summary = qc.run(adata, cfg)
    assert out.n_obs > 0
    assert {"sample_id", "n_bins", "stage"}.issubset(summary.columns)
    assert set(summary["stage"].unique()) == {"pre", "post"}


def test_preprocessing_runs():
    from visium_brain import preprocessing, qc

    adata = _make_synthetic_adata()
    cfg = {
        "project": {"random_seed": 0},
        "qc": {
            "mito_prefix": "mt-", "hb_prefix": "Hb[ab]-",
            "min_counts_per_bin": 1, "min_genes_per_bin": 1,
            "max_pct_mito": 100.0, "min_cells_per_gene": 1,
        },
        "preprocessing": {
            "target_sum": 1e4, "log1p": True,
            "n_top_hvgs": 50, "hvg_flavor": "seurat_v3",
            "scale_max_value": 10.0, "n_pcs": 10,
        },
        "integration": {"batch_key": "sample_id"},
    }
    adata, _ = qc.run(adata, cfg)
    out = preprocessing.run(adata, cfg)
    assert "X_pca" in out.obsm
    assert out.var["highly_variable"].sum() > 0
    assert "counts" in out.layers


def test_pseudobulk_de_shape():
    from visium_brain import differential, preprocessing, qc

    adata = _make_synthetic_adata()
    cfg = {
        "project": {"random_seed": 0},
        "qc": {
            "mito_prefix": "mt-", "hb_prefix": "Hb[ab]-",
            "min_counts_per_bin": 1, "min_genes_per_bin": 1,
            "max_pct_mito": 100.0, "min_cells_per_gene": 1,
        },
        "preprocessing": {
            "target_sum": 1e4, "log1p": True,
            "n_top_hvgs": 50, "hvg_flavor": "seurat_v3",
            "scale_max_value": 10.0, "n_pcs": 10,
        },
        "integration": {"batch_key": "sample_id"},
    }
    adata, _ = qc.run(adata, cfg)
    adata = preprocessing.run(adata, cfg)
    # Fake a single cluster label so pseudobulk works without clustering.
    adata.obs["cell_type"] = pd.Categorical(["clusterA"] * adata.n_obs)
    pb = differential.make_pseudobulk(adata, cluster_key="cell_type", min_cells=10)
    de = differential.pseudobulk_de(pb, cluster_key="cell_type")
    assert pb.n_obs == 4  # one pseudobulk per sample × 1 cluster
    assert not de.empty
    assert {"gene", "pval", "pval_adj", "logfoldchange", "cluster"}.issubset(de.columns)
