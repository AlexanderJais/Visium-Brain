"""Cell-type / region annotation for mouse brain Visium HD bins.

Two strategies are supported:

* `markers`: score canonical mouse brain marker panels with `sc.tl.score_genes`
  and assign each cluster the label of its highest-scoring panel. This is
  robust, has no external model dependency, and works well for region/coarse
  cell-type calls on bin-level Visium HD data.
* `celltypist`: probabilistic per-bin annotation with a pre-trained
  CellTypist mouse brain model. Requires the `celltypist` package and a
  downloaded model file.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

logger = logging.getLogger(__name__)


# Canonical markers for major mouse brain cell types / regions. These are
# intentionally short, well-validated panels — extend as needed.
MOUSE_BRAIN_MARKERS: dict[str, list[str]] = {
    "Excitatory_neuron": ["Slc17a7", "Slc17a6", "Camk2a", "Neurod6", "Satb2"],
    "Inhibitory_neuron": ["Gad1", "Gad2", "Slc32a1", "Pvalb", "Sst", "Vip"],
    "Astrocyte": ["Gfap", "Aqp4", "Slc1a3", "Aldh1l1", "S100b"],
    "Oligodendrocyte": ["Mog", "Mbp", "Plp1", "Mag", "Mobp"],
    "OPC": ["Pdgfra", "Cspg4", "Sox10", "Olig2"],
    "Microglia": ["Cx3cr1", "C1qa", "C1qb", "Tmem119", "P2ry12", "Csf1r"],
    "Endothelial": ["Cldn5", "Pecam1", "Flt1", "Cdh5"],
    "Mural_pericyte": ["Pdgfrb", "Rgs5", "Acta2", "Mcam"],
    "Choroid_plexus": ["Ttr", "Folr1", "Kl"],
    "Ependymal": ["Foxj1", "Ccdc153", "Rarres2"],
    "Granule_cell_DG": ["Prox1", "Calb1", "C1ql2"],
    "CA1_pyramidal": ["Wfs1", "Fibcd1", "Pou3f1"],
    "CA3_pyramidal": ["Cpne7", "Iyd"],
    "Cortical_L2_3": ["Cux2", "Rorb"],
    "Cortical_L5": ["Bcl11b", "Fezf2"],
    "Cortical_L6": ["Foxp2", "Tbr1"],
    "Striatum_MSN": ["Drd1", "Drd2", "Ppp1r1b", "Tac1"],
    "Thalamus": ["Tcf7l2", "Plekhg1"],
}


def score_markers(adata: ad.AnnData, panels: dict[str, list[str]] | None = None) -> list[str]:
    """Score each marker panel; returns list of score column names."""
    panels = panels or MOUSE_BRAIN_MARKERS
    score_cols = []
    for label, genes in panels.items():
        present = [g for g in genes if g in adata.var_names]
        if len(present) < 2:
            logger.warning("Skipping %s (only %d/%d markers present)", label, len(present), len(genes))
            continue
        col = f"score_{label}"
        sc.tl.score_genes(adata, gene_list=present, score_name=col, use_raw=False)
        score_cols.append(col)
    return score_cols


def assign_cluster_labels(
    adata: ad.AnnData,
    score_cols: list[str],
    cluster_key: str = "leiden",
    out_key: str = "cell_type",
) -> ad.AnnData:
    """Assign each cluster the label whose mean score is highest."""
    df = adata.obs[[cluster_key] + score_cols].copy()
    means = df.groupby(cluster_key, observed=True)[score_cols].mean()
    best = means.idxmax(axis=1).str.replace("score_", "", regex=False)
    mapping = best.to_dict()
    adata.obs[out_key] = adata.obs[cluster_key].map(mapping).astype("category")
    # Persist the per-cluster score table for downstream review.
    adata.uns["cluster_marker_scores"] = means.reset_index()
    logger.info("Assigned %d cluster labels", len(mapping))
    return adata


def run_celltypist(adata: ad.AnnData, model: str, out_key: str = "cell_type") -> ad.AnnData:
    try:
        import celltypist  # type: ignore
        from celltypist import models  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "celltypist is not installed. `pip install celltypist` to use this annotator."
        ) from exc

    models.download_models(model=model, force_update=False)
    pred = celltypist.annotate(adata, model=model, majority_voting=True)
    adata.obs[out_key] = pred.predicted_labels["majority_voting"].astype("category")
    return adata


def export_marker_tables(
    adata: ad.AnnData,
    out_dir: str | Path,
    cluster_key: str = "leiden",
    n_top: int = 25,
    method: str = "wilcoxon",
) -> Path:
    """Compute per-cluster top markers and dump them for manual review.

    Writes:
      * ``markers_<cluster_key>.csv`` — long-format table of top markers
      * ``cluster_labels.yaml`` — pre-filled mapping ``cluster_id -> "TBD"``
        for the user to edit
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method=method, use_raw=False)
    res = adata.uns["rank_genes_groups"]
    groups = list(res["names"].dtype.names)
    rows = []
    for g in groups:
        for r in range(min(n_top, len(res["names"][g]))):
            rows.append(
                {
                    "cluster": g,
                    "rank": r + 1,
                    "gene": res["names"][g][r],
                    "logfoldchange": float(res["logfoldchanges"][g][r]),
                    "pval_adj": float(res["pvals_adj"][g][r]),
                    "score": float(res["scores"][g][r]),
                }
            )
    df = pd.DataFrame(rows)
    csv_path = out_dir / f"markers_{cluster_key}.csv"
    df.to_csv(csv_path, index=False)

    label_template = {str(g): "TBD" for g in groups}
    yaml_path = out_dir / "cluster_labels.yaml"
    if not yaml_path.exists():
        with open(yaml_path, "w") as fh:
            yaml.safe_dump(
                {"cluster_key": cluster_key, "labels": label_template},
                fh,
                sort_keys=False,
            )
    logger.info("Exported %d markers x %d clusters -> %s", n_top, len(groups), csv_path)
    return csv_path


def apply_manual_labels(
    adata: ad.AnnData,
    mapping_path: str | Path,
    out_key: str = "cell_type",
) -> ad.AnnData:
    """Apply a user-edited ``cluster_labels.yaml`` to an AnnData.

    The YAML must have the form::

        cluster_key: leiden_l1
        labels:
          "0": Excitatory_neuron
          "1": Astrocyte
          ...
    """
    with open(mapping_path) as fh:
        spec = yaml.safe_load(fh)
    cluster_key = spec["cluster_key"]
    labels = {str(k): str(v) for k, v in spec["labels"].items()}
    if cluster_key not in adata.obs:
        raise KeyError(f"adata.obs[{cluster_key!r}] not found.")

    series = adata.obs[cluster_key].astype(str).map(labels)
    if series.isna().any():
        missing = sorted(set(adata.obs[cluster_key].astype(str)) - set(labels))
        raise ValueError(f"Manual label mapping is incomplete; missing keys: {missing}")
    adata.obs[out_key] = pd.Categorical(series.values)
    logger.info("Applied %d manual labels -> obs['%s']", len(labels), out_key)
    return adata


def run(adata: ad.AnnData, cfg: dict[str, Any]) -> ad.AnnData:
    ann = cfg.get("annotation", {})
    method = ann.get("method", "markers").lower()
    if method == "celltypist":
        return run_celltypist(adata, model=ann.get("celltypist_model", "Mouse_Whole_Brain.pkl"))
    if method == "manual":
        mapping = ann.get("manual_labels_l1") or ann.get("manual_labels")
        if not mapping:
            raise ValueError(
                "annotation.method='manual' requires annotation.manual_labels_l1 "
                "to point at a cluster_labels.yaml."
            )
        return apply_manual_labels(adata, mapping, out_key="cell_type")
    score_cols = score_markers(adata)
    return assign_cluster_labels(adata, score_cols=score_cols, cluster_key="leiden")
