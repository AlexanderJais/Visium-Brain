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


def test_sketch_and_propagation():
    """Leverage-score sketch + kNN propagation reproduces labels on a
    synthetic dataset where each sample is also its own latent cluster."""
    from visium_brain import preprocessing, qc, sketch as sketch_mod

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

    scores = sketch_mod.compute_leverage_scores(adata)
    assert scores.shape == (adata.n_obs,)
    assert np.all(scores >= 0)

    idx = sketch_mod.sketch(adata, fraction=0.3, method="leverage_score", seed=0)
    assert 0 < len(idx) < adata.n_obs
    assert adata.obs["sketch"].sum() == len(idx)

    # Use sample_id as a stand-in label and verify kNN propagation recovers it.
    sub_labels = adata.obs["sample_id"].astype(str).values[idx]
    sketch_mod.propagate_labels(adata, idx, sub_labels, "sample_id_propagated",
                                 use_rep="X_pca", n_neighbors=5)
    acc = (adata.obs["sample_id"].astype(str).values
           == adata.obs["sample_id_propagated"].astype(str).values).mean()
    assert acc > 0.9


def test_hierarchical_l2():
    from visium_brain import annotation, clustering, hierarchical, preprocessing, qc

    adata = _make_synthetic_adata(n_per_sample=100, n_genes=200)
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
        "clustering": {
            "use_rep": "X_pca",
            "l1": {"n_neighbors": 10, "n_pcs": 10, "resolution": 0.5},
            "l2": {"n_neighbors": 5, "n_pcs": 10, "resolution": 0.3},
        },
        "annotation": {"method": "markers", "l2_min_cells": 20},
    }
    adata, _ = qc.run(adata, cfg)
    adata = preprocessing.run(adata, cfg)
    adata.uns["use_rep"] = "X_pca"
    clustering.run_level(adata, cfg, level="l1", use_rep="X_pca", do_umap=False)
    # Provide a synthetic L1 label so L2 has something to recluster within.
    adata.obs["cell_type_l1"] = adata.obs["sample_id"].astype(str).astype("category")
    hierarchical.run_l2(adata, cfg)
    assert "cell_type_l2" in adata.obs
    assert adata.obs["cell_type_l2"].nunique() >= adata.obs["cell_type_l1"].nunique()


def test_manual_label_application(tmp_path):
    import yaml as _yaml

    from visium_brain import annotation

    adata = _make_synthetic_adata(n_per_sample=50, n_genes=100)
    adata.obs["leiden_l1"] = pd.Categorical(
        np.array(["0", "1"])[(np.arange(adata.n_obs) % 2)]
    )
    mapping = tmp_path / "labels.yaml"
    with open(mapping, "w") as fh:
        _yaml.safe_dump(
            {"cluster_key": "leiden_l1", "labels": {"0": "TypeA", "1": "TypeB"}}, fh
        )
    annotation.apply_manual_labels(adata, mapping, out_key="cell_type_manual")
    assert set(adata.obs["cell_type_manual"].astype(str)) == {"TypeA", "TypeB"}


def test_export_marker_tables(tmp_path):
    from visium_brain import annotation, clustering, preprocessing, qc

    adata = _make_synthetic_adata(n_per_sample=80, n_genes=150)
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
        "clustering": {
            "use_rep": "X_pca",
            "l1": {"n_neighbors": 10, "n_pcs": 10, "resolution": 0.5},
        },
    }
    adata, _ = qc.run(adata, cfg)
    adata = preprocessing.run(adata, cfg)
    clustering.run_level(adata, cfg, level="l1", use_rep="X_pca", do_umap=False)
    csv_path = annotation.export_marker_tables(adata, tmp_path, cluster_key="leiden_l1", n_top=5)
    assert csv_path.exists()
    yaml_path = tmp_path / "cluster_labels.yaml"
    assert yaml_path.exists()


def test_clustering_skip_neighbors_explicit():
    """Regression for B2: ``clustering.run_level`` must build the
    neighborhood graph based on the explicit ``skip_neighbors``
    argument, not by reading ``integration.method`` from config.

    Previously the function silently skipped ``sc.pp.neighbors`` when
    ``cfg['integration']['method'] == 'bbknn'``; combined with the
    sketch path this produced an empty graph on the sketch (BBKNN had
    only built obsp on the full adata).
    """
    from visium_brain import clustering, preprocessing, qc

    adata = _make_synthetic_adata(n_per_sample=80, n_genes=150)
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
        # Set integration.method to a value the OLD code path would
        # interpret as "skip neighbors" -- we want to confirm run_level
        # ignores this when called with the default skip_neighbors=False.
        "integration": {"method": "bbknn", "batch_key": "sample_id"},
        "clustering": {
            "use_rep": "X_pca",
            "l1": {"n_neighbors": 5, "n_pcs": 10, "resolution": 0.5},
        },
    }
    adata, _ = qc.run(adata, cfg)
    adata = preprocessing.run(adata, cfg)

    # Default skip_neighbors=False: must build the graph regardless of
    # integration.method.
    clustering.run_level(adata, cfg, level="l1", use_rep="X_pca", do_umap=False)
    assert "neighbors" in adata.uns
    assert "connectivities" in adata.obsp
    assert adata.obsp["connectivities"].nnz > 0
    assert "leiden_l1" in adata.obs


def test_log1p_preserved_in_raw_and_de_uses_it():
    """Regression for B1: ``sc.pp.scale`` must not contaminate the values
    DE / marker scoring see. After preprocessing:

    * ``adata.raw.X`` holds log1p (small magnitudes, non-negative).
    * ``adata.X`` for HVGs has been z-scored (clipped to ``scale_max_value``).
    * ``rank_genes_groups(... use_raw=True)`` recovers the planted
      treatment effect (the synthetic adata bumps genes 2..9 up in the
      'treated' samples).
    """
    from visium_brain import differential, preprocessing, qc

    adata = _make_synthetic_adata(n_per_sample=200, n_genes=300)
    cfg = {
        "project": {"random_seed": 0},
        "qc": {
            "mito_prefix": "mt-", "hb_prefix": "Hb[ab]-",
            "min_counts_per_bin": 1, "min_genes_per_bin": 1,
            "max_pct_mito": 100.0, "min_cells_per_gene": 1,
        },
        "preprocessing": {
            "target_sum": 1e4, "log1p": True,
            "n_top_hvgs": 100, "hvg_flavor": "seurat_v3",
            "scale_max_value": 10.0, "n_pcs": 10,
        },
        "integration": {"batch_key": "sample_id"},
    }
    adata, _ = qc.run(adata, cfg)
    adata = preprocessing.run(adata, cfg)

    # 1. .raw populated with log1p (non-negative, modest magnitudes).
    assert adata.raw is not None, "preprocessing must snapshot log1p into adata.raw"
    raw_X = adata.raw.X.toarray() if sp.issparse(adata.raw.X) else np.asarray(adata.raw.X)
    assert raw_X.min() >= 0.0
    assert raw_X.max() < 30.0, f"raw values look unscaled ({raw_X.max():.2f})"

    # 2. .X has been scaled on HVGs (mean ~0, clipped at scale_max_value).
    hvg = adata.var["highly_variable"].values
    X_hvg = adata.X[:, hvg]
    X_hvg = X_hvg.toarray() if sp.issparse(X_hvg) else np.asarray(X_hvg)
    assert X_hvg.max() <= 10.0 + 1e-6
    assert abs(X_hvg.mean()) < 1.0, f"HVG mean {X_hvg.mean():.3f} doesn't look z-scored"

    # 3. DE with use_raw=True recovers the planted up-regulation in
    # treated samples on genes 2..9 (gene_0 / gene_1 are renamed to mt-*
    # by _make_synthetic_adata so the planted set in raw-name space is
    # gene_2..gene_9).
    adata.obs["cell_type"] = pd.Categorical(["A"] * adata.n_obs)
    de = differential.condition_de(
        adata, groupby="condition", reference="control",
        method="wilcoxon", cluster_key="cell_type",
    )
    assert not de.empty
    treated = de[de["group"] == "treated"]
    planted = {f"gene_{i}" for i in range(2, 10)}
    top50 = set(treated.sort_values("score", ascending=False)["gene"].head(50))
    hits = planted & top50
    assert len(hits) >= 4, (
        f"Expected most planted genes in top 50 by score; got {sorted(hits)}"
    )
    # logFC for those hits should be positive (treated > control).
    planted_lfc = treated[treated["gene"].isin(planted)]["logfoldchange"]
    assert (planted_lfc > 0).sum() >= 4


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
