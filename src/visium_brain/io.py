"""Input/output utilities for Visium HD Space Ranger outputs.

Visium HD Space Ranger writes binned outputs at three resolutions:

    <sample>/binned_outputs/square_002um/
    <sample>/binned_outputs/square_008um/  (recommended starting point)
    <sample>/binned_outputs/square_016um/

Each bin directory contains a `filtered_feature_bc_matrix.h5` and a
`spatial/` folder. We use scanpy's Visium reader against that bin
directory; the read function expects the `tissue_positions.parquet`
or `tissue_positions_list.csv` layout which Space Ranger >= 3.0
emits for HD.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable

import anndata as ad
import scanpy as sc

from .utils import Sample

logger = logging.getLogger(__name__)


def _resolve_bin_dir(sample: Sample, bin_size: str) -> Path:
    """Find the binned-outputs directory for a sample.

    Accepts either a path to the top-level Space Ranger output (containing
    `binned_outputs/`) or directly to a `square_XXXum` directory.
    """
    p = sample.path
    if (p / "binned_outputs" / bin_size).exists():
        return p / "binned_outputs" / bin_size
    if p.name == bin_size:
        return p
    if (p / bin_size).exists():
        return p / bin_size
    raise FileNotFoundError(
        f"Could not locate {bin_size} directory under {p}. "
        f"Expected `<sample>/binned_outputs/{bin_size}/` from Space Ranger."
    )


def read_visium_hd_sample(sample: Sample, bin_size: str = "square_008um") -> ad.AnnData:
    """Load a single Visium HD sample at the given bin resolution."""
    bin_dir = _resolve_bin_dir(sample, bin_size)
    logger.info("Loading %s from %s", sample.sample_id, bin_dir)

    adata = sc.read_visium(
        path=str(bin_dir),
        count_file="filtered_feature_bc_matrix.h5",
        load_images=True,
    )
    adata.var_names_make_unique()

    # Stash sample-level metadata
    adata.obs["sample_id"] = sample.sample_id
    adata.obs["condition"] = sample.condition
    adata.obs["mouse_id"] = sample.mouse_id
    adata.obs["bin_size"] = bin_size
    adata.obs_names = [f"{sample.sample_id}_{bc}" for bc in adata.obs_names]

    # Re-key spatial metadata under the sample id so multi-sample objects keep
    # their per-sample images / scale factors.
    if "spatial" in adata.uns:
        old_keys = list(adata.uns["spatial"].keys())
        if old_keys and sample.sample_id not in adata.uns["spatial"]:
            adata.uns["spatial"][sample.sample_id] = adata.uns["spatial"].pop(old_keys[0])

    return adata


def read_all_samples(samples: Iterable[Sample], bin_size: str) -> ad.AnnData:
    """Load and concatenate all samples into one AnnData."""
    adatas = {s.sample_id: read_visium_hd_sample(s, bin_size=bin_size) for s in samples}
    logger.info("Concatenating %d samples", len(adatas))

    merged = ad.concat(
        adatas,
        axis=0,
        join="outer",
        label="sample_id",
        keys=list(adatas.keys()),
        index_unique=None,
        merge="unique",
        uns_merge="unique",
    )
    # ad.concat drops uns["spatial"]; re-attach manually so plotting works.
    merged.uns["spatial"] = {}
    for sid, a in adatas.items():
        if "spatial" in a.uns and sid in a.uns["spatial"]:
            merged.uns["spatial"][sid] = a.uns["spatial"][sid]

    merged.obs["condition"] = merged.obs["condition"].astype("category")
    merged.obs["sample_id"] = merged.obs["sample_id"].astype("category")
    merged.obs["mouse_id"] = merged.obs["mouse_id"].astype("category")
    return merged


def write_h5ad(adata: ad.AnnData, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Writing %s (%s)", path, adata.shape)
    adata.write_h5ad(path, compression="gzip")


def read_h5ad(path: str | Path) -> ad.AnnData:
    logger.info("Reading %s", path)
    return ad.read_h5ad(path)
