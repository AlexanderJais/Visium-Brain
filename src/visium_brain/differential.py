"""Differential expression: per cluster, per condition, and pseudobulk.

For Visium HD with biological replicates (2 mice per condition in the
default design), the recommended approach is pseudobulk DE: sum raw
counts across bins per (sample, cluster), then test condition effects
*with mice as the unit of replication*. This module provides:

* :func:`cluster_markers` - Wilcoxon markers per cluster vs the rest
  (fast, bin-level).
* :func:`condition_de` - Wilcoxon between conditions per cluster
  (bin-level; **heavily pseudoreplicated** -- use only for
  exploration / ranking, never as a significance test).
* :func:`make_pseudobulk` + :func:`pseudobulk_de` - the
  replication-correct path.

The :func:`run` orchestrator takes a single ``cluster_key`` for every
flavor (``cell_type`` if present, else ``leiden``); when hierarchical
L1->L2 has run, ``cell_type`` is the L2 label, so the bundle reports
DE at L2 by default and additionally emits an L1-keyed
``cluster_markers_l1`` table for side-by-side comparison.

**Important caveat on `pseudobulk_de` p-values.** The Wilcoxon
rank-sum test has a hard lower bound on its achievable p-value that
is determined entirely by group sizes:

    n1 = n2 = 2  -->  two-sided p_min = 2/C(4,2)  = 2/6  ~ 0.333
    n1 = n2 = 3  -->  two-sided p_min = 2/C(6,3)  = 2/20 = 0.100
    n1 = n2 = 4  -->  two-sided p_min = 2/C(8,4)  = 2/70 ~ 0.029

With 2 mice per condition, *no* gene -- not even a perfectly separated
one -- can reach p < 0.05. Reporting p-values in that regime would
mislead readers. We therefore gate :func:`pseudobulk_de`: if the
smaller condition has fewer than ``min_per_group`` (default 3)
pseudobulks for a given cluster, the ``pval`` and ``pval_adj`` columns
are blanked to ``NaN`` for that cluster, a ``low_power=True`` flag is
set, and a warning is logged. The effect-size column
(``logfoldchange``) and the per-cluster pseudobulk counts are always
populated and remain useful for prioritization. For a real
significance test in this regime, use a parametric pseudobulk model
(``pyDESeq2`` or ``edgeR`` via ``rpy2``).
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
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method, use_raw=True)
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
        sc.tl.rank_genes_groups(adata, groupby=groupby, reference=reference, method=method, use_raw=True)
        df = _rank_genes_to_df(adata)
        df["cluster"] = "all"
        return df

    for c in adata.obs[cluster_key].cat.categories:
        # Materialize the slice: sc.tl.rank_genes_groups writes to
        # sub.uns["rank_genes_groups"], and recent scanpy raises
        # ImplicitModificationWarning (and on some versions silently
        # writes to a transient copy) when handed a view.
        sub = adata[adata.obs[cluster_key] == c].copy()
        if sub.obs[groupby].nunique() < 2:
            continue
        sc.tl.rank_genes_groups(sub, groupby=groupby, reference=reference, method=method, use_raw=True)
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
    min_per_group: int = 3,
) -> pd.DataFrame:
    """Wilcoxon test on log-CPM pseudobulks per cluster across conditions.

    Per cluster, builds log-CPM from the (sample, cluster) summed
    counts in ``pb`` and runs ``sc.tl.rank_genes_groups`` between
    conditions using the pseudobulks (mice) as the unit of
    replication.

    The ``min_per_group`` argument gates p-value reporting on
    statistical power. **Read the module-level docstring before
    interpreting the output**:

    * If the smaller condition has at least ``min_per_group``
      pseudobulks for a cluster, the cluster's rows carry usual
      ``pval`` / ``pval_adj`` along with ``logfoldchange``.
    * If it has fewer, ``pval`` and ``pval_adj`` are set to ``NaN``
      and ``low_power=True`` is stamped on every row of that cluster.
      A warning is logged. Effect sizes (``logfoldchange``) and group
      sizes (``min_n_per_group``, ``n_pseudobulks``) are always
      reported.

    The default threshold is 3 because Wilcoxon's two-sided minimum
    achievable p with n1 = n2 = 2 is 2/C(4,2) = 0.333, with n1 = n2 = 3
    is 2/C(6,3) = 0.10, and only with n1 = n2 >= 4 can you cross
    p < 0.05. For a real test in the underpowered regime use
    ``pyDESeq2`` / ``edgeR`` (via ``rpy2``) on the same pseudobulks.

    Parameters
    ----------
    pb
        Pseudobulk AnnData from :func:`make_pseudobulk` (raw counts in
        ``.X``).
    cluster_key, condition_key
        Columns of ``pb.obs`` defining strata and the test factor.
    reference
        Level of ``condition_key`` used as the Wilcoxon reference; the
        other level becomes the test ``group`` and its
        ``logfoldchange`` is signed positive when up in the test
        condition.
    min_per_group
        Minimum smaller-group size required to publish p-values. Set
        to 0 to always report p-values (not recommended for
        publication; useful for diagnostics).

    Returns
    -------
    Long-format DataFrame with columns::

        cluster, group, gene, logfoldchange, pval, pval_adj, score,
        n_pseudobulks, min_n_per_group, low_power
    """
    # log-CPM normalize the pseudobulks. We work on a fresh AnnData so
    # the input ``pb`` (raw counts) is untouched and re-usable.
    X = pb.X.astype(float)
    libsize = X.sum(axis=1, keepdims=True)
    libsize[libsize == 0] = 1
    logcpm = np.log2(1e6 * X / libsize + 1)
    pb2 = ad.AnnData(X=logcpm, obs=pb.obs.copy(), var=pb.var.copy())

    frames = []
    for c in pb2.obs[cluster_key].unique():
        # Materialize the slice before sc.tl.rank_genes_groups: see the
        # matching note in condition_de. Handing rank_genes_groups a
        # view triggers ImplicitModificationWarning on recent scanpy and
        # on older versions silently writes uns into a transient copy.
        sub = pb2[pb2.obs[cluster_key] == c].copy()
        if sub.obs[condition_key].nunique() < 2:
            continue
        # Per-condition pseudobulk counts in this cluster. The Wilcoxon
        # achievable minimum p depends only on these two numbers; see
        # the module docstring for the C(n1+n2, n1) derivation.
        per_group = sub.obs[condition_key].value_counts()
        min_n = int(per_group.min())

        # use_raw=False here: pb2 carries log-CPM directly in .X and
        # has no .raw, so the standard preprocessing-time use_raw=True
        # convention does NOT apply here.
        sc.tl.rank_genes_groups(
            sub, groupby=condition_key, reference=reference,
            method="wilcoxon", use_raw=False,
        )
        df = _rank_genes_to_df(sub)
        df["cluster"] = str(c)
        df["n_pseudobulks"] = sub.n_obs
        df["min_n_per_group"] = min_n
        df["low_power"] = bool(min_n < min_per_group)

        if min_n < min_per_group:
            # Wilcoxon p-values are not informative below the threshold:
            # the smallest two-sided p achievable with n1 = n2 = min_n
            # is 2 / C(2*min_n, min_n), which exceeds 0.05 for min_n < 4.
            # Blank them so a downstream reader cannot mistake the
            # output for a real significance test, but keep
            # logfoldchange (effect size) which remains interpretable.
            df.loc[:, "pval"] = np.nan
            df.loc[:, "pval_adj"] = np.nan
            logger.warning(
                "[pseudobulk_de] cluster=%s: smaller-group n=%d < min_per_group=%d; "
                "p-values blanked (Wilcoxon cannot reach p<0.05 at this n). "
                "logfoldchange retained.",
                c, min_n, min_per_group,
            )
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def run(adata: ad.AnnData, cfg: dict[str, Any]) -> dict[str, pd.DataFrame]:
    de_cfg = cfg["differential_expression"]
    out: dict[str, pd.DataFrame] = {}

    # Single cluster key for every DE flavor. When hierarchical annotation
    # has run, the pipeline promotes cell_type_l2 into cell_type, so this
    # naturally selects L2; otherwise it picks up L1 cell_type or falls
    # back to leiden. Previously cluster_markers was pinned to "leiden"
    # (== L1 alias) while condition_de / pseudobulk_de used cell_type
    # (L2 after hierarchical), producing DE outputs at two different
    # granularities.
    cluster_key = "cell_type" if "cell_type" in adata.obs else "leiden"
    method = de_cfg.get("method", "wilcoxon")

    out["cluster_markers"] = cluster_markers(adata, groupby=cluster_key, method=method)
    # If hierarchical L1->L2 ran, also dump L1 markers so users get both
    # granularities side by side without having to re-run. cell_type_l1
    # is preserved by run_cluster_annotate even after the L2 promotion.
    # rank_genes_groups requires >= 2 groups, so guard the degenerate
    # single-L1 case (e.g. when annotation returned only "Unknown").
    if (
        "cell_type_l1" in adata.obs
        and "cell_type_l2" in adata.obs
        and adata.obs["cell_type_l1"].nunique() >= 2
    ):
        out["cluster_markers_l1"] = cluster_markers(
            adata, groupby="cell_type_l1", method=method
        )

    out["condition_de_binlevel"] = condition_de(
        adata,
        groupby=de_cfg.get("groupby", "condition"),
        reference=de_cfg.get("reference", "control"),
        method=method,
        cluster_key=cluster_key,
    )

    if de_cfg.get("pseudobulk", True):
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
            min_per_group=de_cfg.get("pseudobulk_min_per_group", 3),
        )
    return out
