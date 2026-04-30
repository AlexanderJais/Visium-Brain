"""Sketch-based analysis for Visium HD.

At full Visium HD scale (millions of 8 µm bins) full-dataset PCA + Leiden
becomes painful. The standard fix in Seurat v5 is sketch-based analysis:
sample a small, *representative* subset of bins, do clustering and
annotation on the sketch, then propagate the labels back to every bin
via a kNN classifier in the integrated PCA space.

The key sampling primitive is **leverage-score sampling**, which
oversamples rare populations (high row norms in the PCA basis come from
bins that are not well-explained by the bulk of the data, i.e. rare cell
types). This matches the procedure used in Oh et al., *Nat Genet* 2025
(https://www.nature.com/articles/s41588-025-02193-3) and the Seurat
``SketchData`` / ``ProjectData`` workflow.
"""

from __future__ import annotations

import logging
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def compute_leverage_scores(adata: ad.AnnData, n_pcs: int | None = None) -> np.ndarray:
    """Statistical leverage scores from the stored PCA.

    Recovers an orthonormal basis ``U`` from scanpy's PCA — scanpy stores
    ``X_pca[:, k] = sqrt(var_k) * U[:, k]`` and the per-component variances
    in ``adata.uns['pca']['variance']`` — and returns the squared row
    norms of ``U`` (the classical leverage score, ``h_ii``).
    """
    if "X_pca" not in adata.obsm:
        raise KeyError("adata.obsm['X_pca'] not found; run preprocessing.run() first.")
    if "pca" not in adata.uns or "variance" not in adata.uns["pca"]:
        raise KeyError("adata.uns['pca']['variance'] missing; was PCA run via scanpy?")

    X_pca = np.asarray(adata.obsm["X_pca"])
    var = np.asarray(adata.uns["pca"]["variance"])
    if n_pcs is not None:
        X_pca = X_pca[:, :n_pcs]
        var = var[:n_pcs]
    sigma = np.sqrt(np.maximum(var, 1e-12))
    U = X_pca / sigma[None, :]
    return np.einsum("ij,ij->i", U, U)


def sketch(
    adata: ad.AnnData,
    fraction: float = 0.15,
    method: str = "leverage_score",
    n_pcs: int | None = None,
    seed: int = 0,
) -> np.ndarray:
    """Pick a sketch of the data and stash the result in ``adata``.

    Returns the sorted integer indices of selected bins; also writes
    ``adata.obs['sketch']`` (bool) and ``adata.uns['sketch_indices']``.
    """
    n = adata.n_obs
    n_sketch = max(1, int(round(fraction * n)))
    rng = np.random.default_rng(seed)

    if method == "uniform":
        idx = rng.choice(n, size=n_sketch, replace=False)
    elif method == "leverage_score":
        scores = compute_leverage_scores(adata, n_pcs=n_pcs)
        p = scores / scores.sum()
        idx = rng.choice(n, size=n_sketch, replace=False, p=p)
    else:
        raise ValueError(f"Unknown sketch method: {method!r}")

    idx = np.sort(idx)
    mask = np.zeros(n, dtype=bool)
    mask[idx] = True
    adata.obs["sketch"] = mask
    adata.uns["sketch_indices"] = idx
    adata.uns["sketch_method"] = method
    adata.uns["sketch_fraction"] = float(n_sketch) / n
    logger.info(
        "Sketch: %d / %d bins (%.1f%%) via %s",
        n_sketch, n, 100.0 * n_sketch / n, method,
    )
    return idx


def propagate_labels(
    adata: ad.AnnData,
    sketch_idx: np.ndarray,
    sketch_labels: np.ndarray | pd.Series,
    out_key: str,
    use_rep: str = "X_pca_harmony",
    n_neighbors: int = 15,
) -> None:
    """Propagate sketch labels to all bins via kNN in ``use_rep`` space."""
    from sklearn.neighbors import KNeighborsClassifier

    if use_rep not in adata.obsm:
        raise KeyError(f"adata.obsm[{use_rep!r}] not found; run integration first.")
    X = np.asarray(adata.obsm[use_rep])
    y = np.asarray(sketch_labels)

    clf = KNeighborsClassifier(n_neighbors=n_neighbors, n_jobs=-1)
    clf.fit(X[sketch_idx], y)
    preds = clf.predict(X)
    adata.obs[out_key] = pd.Categorical(preds)
    logger.info("Propagated %d labels to %d bins -> obs['%s']", len(set(y)), adata.n_obs, out_key)
