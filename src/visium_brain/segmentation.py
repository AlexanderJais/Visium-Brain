"""Nuclear-segmentation-based aggregation for Visium HD 2 µm bins.

Wraps `bin2cell <https://github.com/Teichlab/bin2cell>`_ (Polanski et al.,
2024) to:

1. Load Space Ranger 2 µm-bin output for a sample, together with the
   full-resolution H&E image.
2. Run StarDist on the H&E to segment nuclei.
3. Aggregate 2 µm bins under each nucleus into a per-cell AnnData.

This is the standard route to recover (near-)single-cell resolution from
Visium HD when 8 µm bins under-sample small cell types (microglia,
endothelial, immune). It is opt-in (config: ``segmentation.enabled``)
because it requires the optional ``bin2cell`` dependency and a path to
each sample's full-resolution H&E TIFF, which Space Ranger does not
write into the ``binned_outputs`` tree.

For each sample provide an extra ``he_image_path`` field in
``config.yaml``:

.. code-block:: yaml

    samples:
      - sample_id: ctrl_m1
        path: data/raw/ctrl_m1
        he_image_path: data/raw/ctrl_m1/he_full.tif
        ...

If you do not have full-resolution H&E images, leave ``segmentation``
disabled and continue analysing 8 µm bins.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import anndata as ad

from .io import _resolve_bin_dir
from .utils import Sample

logger = logging.getLogger(__name__)


def _require_bin2cell():
    try:
        import bin2cell as b2c  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "bin2cell is not installed. Install via "
            "`pip install bin2cell` (and ensure stardist + tensorflow are "
            "available) to use segmentation."
        ) from exc
    return b2c


def segment_sample(
    sample: Sample,
    he_image_path: str | Path,
    out_dir: Path,
    he_mpp: float = 0.5,
    stardist_model: str = "2D_versatile_he",
    prob_thresh: float = 0.01,
    nms_thresh: float = 0.5,
    bin_size: str = "square_002um",
) -> ad.AnnData:
    """Run bin2cell segmentation on one sample and return cell-level AnnData."""
    b2c = _require_bin2cell()

    bin_dir = _resolve_bin_dir(sample, bin_size)
    he_image_path = Path(he_image_path)
    if not he_image_path.exists():
        raise FileNotFoundError(f"H&E image not found for {sample.sample_id}: {he_image_path}")
    sample_out = out_dir / sample.sample_id
    sample_out.mkdir(parents=True, exist_ok=True)

    logger.info("[seg] %s: loading 2 µm bins from %s", sample.sample_id, bin_dir)
    adata = b2c.read_visium(
        path=str(bin_dir),
        source_image_path=str(he_image_path),
    )
    adata.var_names_make_unique()

    he_scaled = sample_out / "he_scaled.tiff"
    labels_npz = sample_out / "labels_he.npz"

    logger.info("[seg] %s: rescaling H&E to mpp=%s", sample.sample_id, he_mpp)
    b2c.scaled_he_image(adata, mpp=he_mpp, save_path=str(he_scaled))

    logger.info("[seg] %s: StarDist (%s, prob=%.3f)", sample.sample_id, stardist_model, prob_thresh)
    b2c.stardist(
        image_path=str(he_scaled),
        labels_npz_path=str(labels_npz),
        stardist_model=stardist_model,
        prob_thresh=prob_thresh,
        nms_thresh=nms_thresh,
    )
    b2c.insert_labels(
        adata,
        labels_npz_path=str(labels_npz),
        basis="spatial",
        spatial_key="spatial",
        mpp=he_mpp,
        labels_key="labels_he",
    )

    logger.info("[seg] %s: aggregating bins under nuclei", sample.sample_id)
    cdata = b2c.bin_to_cell(adata, labels_key="labels_he")

    cdata.obs["sample_id"] = sample.sample_id
    cdata.obs["condition"] = sample.condition
    cdata.obs["mouse_id"] = sample.mouse_id
    cdata.obs_names = [f"{sample.sample_id}_cell_{bc}" for bc in cdata.obs_names]

    out_h5ad = sample_out / "adata_cells.h5ad"
    cdata.write_h5ad(out_h5ad, compression="gzip")
    logger.info("[seg] %s: %d cells -> %s", sample.sample_id, cdata.n_obs, out_h5ad)
    return cdata


def run(cfg: dict[str, Any]) -> ad.AnnData:
    """Run segmentation for every sample with a configured H&E image."""
    seg = cfg.get("segmentation", {})
    if not seg.get("enabled", False):
        raise RuntimeError("segmentation.enabled is false; nothing to do.")

    out_dir = Path(seg.get("output_per_sample", "results/segmentation"))
    out_dir.mkdir(parents=True, exist_ok=True)
    bin_size = seg.get("bin_size", "square_002um")

    cell_adatas: dict[str, ad.AnnData] = {}
    for s in cfg["samples"]:
        sample = Sample.from_dict(s)
        he = s.get("he_image_path")
        if not he:
            logger.warning("[seg] %s: no he_image_path; skipping", sample.sample_id)
            continue
        cell_adatas[sample.sample_id] = segment_sample(
            sample,
            he_image_path=he,
            out_dir=out_dir,
            he_mpp=seg.get("he_mpp", 0.5),
            stardist_model=seg.get("stardist_model", "2D_versatile_he"),
            prob_thresh=seg.get("prob_thresh", 0.01),
            nms_thresh=seg.get("nms_thresh", 0.5),
            bin_size=bin_size,
        )

    if not cell_adatas:
        raise RuntimeError("No samples with he_image_path found in config.")

    # See io.read_all_samples: each cell-level adata already has
    # obs["sample_id"] set, so no `label=` here.
    merged = ad.concat(
        cell_adatas, axis=0, join="outer",
        index_unique=None, merge="unique",
    )
    merged.write_h5ad(out_dir / "adata_cells_merged.h5ad", compression="gzip")
    logger.info("[seg] merged cell-level AnnData: %s", merged.shape)
    return merged
