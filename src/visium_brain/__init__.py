"""Analysis pipeline for 10x Genomics Visium HD mouse-brain data.

Public entry points are organised one-stage-per-module:

* :mod:`visium_brain.io`            -- read Space Ranger output
* :mod:`visium_brain.qc`            -- QC metrics and filtering
* :mod:`visium_brain.preprocessing` -- normalize / HVG / PCA
* :mod:`visium_brain.integration`   -- Harmony / BBKNN / scVI
* :mod:`visium_brain.sketch`        -- leverage-score sketching + kNN propagation
* :mod:`visium_brain.clustering`    -- per-level neighbors / UMAP / Leiden
* :mod:`visium_brain.annotation`    -- marker-score / CellTypist / manual labels
* :mod:`visium_brain.hierarchical`  -- L1 -> L2 reclustering and annotation
* :mod:`visium_brain.spatial`       -- Moran's I, neighborhood enrichment
* :mod:`visium_brain.differential`  -- cluster / condition / pseudobulk DE
* :mod:`visium_brain.segmentation`  -- bin2cell 2 µm + H&E nuclear segmentation
* :mod:`visium_brain.plotting`      -- QC, UMAP, spatial, DE figures
* :mod:`visium_brain.pipeline`      -- stage orchestration used by the CLI
* :mod:`visium_brain.cli`           -- the ``visium-brain`` command

See ``README.md`` for the high-level workflow.
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("visium-brain")
except PackageNotFoundError:  # pragma: no cover - package not installed
    __version__ = "0.1.0"

__all__ = ["__version__"]
