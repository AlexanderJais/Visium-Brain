"""Multi-sample integration / batch correction.

Defaults to Harmony on top of PCA, since it's fast, well-suited to a
small number of samples (4 here), and supported out-of-the-box by
scanpy's external API. Falls back gracefully if optional backends
are not installed.
"""

from __future__ import annotations

import logging
from typing import Any

import anndata as ad

logger = logging.getLogger(__name__)


def run_harmony(
    adata: ad.AnnData,
    batch_key: str = "sample_id",
    max_iter_harmony: int = 20,
    theta: float = 2.0,
) -> None:
    import scanpy.external as sce

    sce.pp.harmony_integrate(
        adata,
        key=batch_key,
        max_iter_harmony=max_iter_harmony,
        theta=theta,
    )
    logger.info("Harmony integration done; embedding stored in obsm['X_pca_harmony']")


def run_bbknn(adata: ad.AnnData, batch_key: str = "sample_id", n_pcs: int = 50) -> None:
    import scanpy.external as sce

    sce.pp.bbknn(adata, batch_key=batch_key, n_pcs=n_pcs)
    logger.info("BBKNN graph built and stored in obsp")


def run_scvi(adata: ad.AnnData, batch_key: str = "sample_id", seed: int = 0) -> None:
    try:
        import scvi  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "scvi-tools is not installed. Install it via `pip install scvi-tools` "
            "to use integration.method=scvi."
        ) from exc

    scvi.settings.seed = seed
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_key)
    model = scvi.model.SCVI(adata)
    model.train()
    adata.obsm["X_scVI"] = model.get_latent_representation()
    logger.info("scVI integration done; latent stored in obsm['X_scVI']")


def run(adata: ad.AnnData, cfg: dict[str, Any]) -> str:
    """Apply the configured integration and return the embedding key to use downstream."""
    integ = cfg.get("integration", {})
    method = integ.get("method", "harmony").lower()
    batch_key = integ.get("batch_key", "sample_id")

    if method == "none":
        logger.info("Integration disabled; using X_pca")
        return "X_pca"
    if method == "harmony":
        h = integ.get("harmony", {})
        run_harmony(
            adata,
            batch_key=batch_key,
            max_iter_harmony=h.get("max_iter_harmony", 20),
            theta=h.get("theta", 2.0),
        )
        return "X_pca_harmony"
    if method == "bbknn":
        run_bbknn(adata, batch_key=batch_key, n_pcs=cfg["preprocessing"].get("n_pcs", 50))
        return "X_pca"  # BBKNN builds a graph directly; neighbors step is skipped later.
    if method == "scvi":
        run_scvi(adata, batch_key=batch_key, seed=cfg["project"].get("random_seed", 0))
        return "X_scVI"
    raise ValueError(f"Unknown integration method: {method}")
