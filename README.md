# Visium-Brain

A Python pipeline for analyzing 10x Genomics **Visium HD** spatial gene
expression data on mouse brain sections, designed around a 2-condition ×
2-mice study (4 sections total) but configurable for any number of
samples.

Built on `scanpy` + `anndata` + `squidpy`, with optional `harmonypy`,
`celltypist`, `scvi-tools`, and `bin2cell` integrations. The defaults
reproduce the analysis recipe from Oh et al., *Nat Genet* 2025
([10.1038/s41588-025-02193-3](https://www.nature.com/articles/s41588-025-02193-3))
— 8 µm bins, leverage-score sketching, hierarchical L1 → L2 annotation —
and add cross-condition support (Harmony integration, pseudobulk DE)
that the original paper did not need.

---

## Highlights

- **Visium HD aware.** Loads Space Ranger ≥ 3.0 `binned_outputs/`
  directly, preserves per-sample images and scale factors through
  multi-sample concatenation.
- **Sketch-based clustering.** Leverage-score sampling (oversamples
  rare populations) → cluster + annotate on the sketch → kNN
  propagation to all bins. Mirrors Seurat v5 `SketchData` / `ProjectData`.
- **Hierarchical annotation.** L1 broad cell types, then per-L1
  reclustering at lower resolution to refine subtypes (L2). Each level
  has its own `n_pcs` / `resolution`.
- **Replicate-aware DE.** Pseudobulk DE per (sample × cluster) with
  mice as the unit of replication, in addition to fast bin-level
  Wilcoxon tests.
- **Nuclear segmentation.** Optional `bin2cell` route to aggregate
  2 µm bins under StarDist-segmented nuclei for cell-resolution
  analyses of small cell types.
- **Manual annotation hook.** `visium-brain export-markers` dumps top
  markers per cluster + a `cluster_labels.yaml` stub the pipeline
  consumes via `annotation.method: manual`.
- **Spatial statistics.** Per-sample neighborhood graphs, Moran's I on
  HVGs, neighborhood enrichment between cell types (squidpy).

---

## Quick start

```bash
conda env create -f environment.yml && conda activate visium-brain
pip install -e .

# point config/config.yaml at your Space Ranger outputs, then:
visium-brain --config config/config.yaml all
```

Outputs land in `results/`; intermediate `.h5ad` per stage so any step
can be re-run on its own.

---

## Install

```bash
# Recommended: conda for igraph/leidenalg/opencv pinning
conda env create -f environment.yml
conda activate visium-brain
pip install -e .

# Or pure pip (in a virtualenv)
python -m venv .venv && source .venv/bin/activate
pip install -e .
```

Optional extras:

| Extra | Adds | Install |
|---|---|---|
| `segmentation` | `bin2cell` + `stardist` for 2 µm nuclear segmentation | `pip install -e ".[segmentation]"` |
| `celltypist` | Pretrained CellTypist annotators | `pip install -e ".[celltypist]"` |
| `scvi` | scVI-based integration | `pip install -e ".[scvi]"` |
| `spatialdata` | `spatialdata` / `spatialdata-io` HD readers | `pip install -e ".[spatialdata]"` |
| `dev` | `pytest`, `ruff`, `jupyter` | `pip install -e ".[dev]"` |

---

## Configure

Edit `config/config.yaml`. The minimal required edits are the `samples:`
list — every sample points at a Space Ranger top-level directory
containing `binned_outputs/square_008um/`:

```yaml
samples:
  - sample_id: ctrl_m1
    condition: control
    mouse_id: M1
    path: data/raw/ctrl_m1
  - sample_id: ctrl_m2
    condition: control
    mouse_id: M2
    path: data/raw/ctrl_m2
  - sample_id: trt_m3
    condition: treated
    mouse_id: M3
    path: data/raw/trt_m3
  - sample_id: trt_m4
    condition: treated
    mouse_id: M4
    path: data/raw/trt_m4
```

Key knobs in the same file (defaults shown):

| Block | Setting | Default | Notes |
|---|---|---|---|
| top-level | `bin_size` | `square_008um` | use `square_002um` for nuclear segmentation only |
| `qc` | `min_counts_per_bin` / `min_genes_per_bin` / `max_pct_mito` | 50 / 20 / 25 | mouse-tuned (`mito_prefix: mt-`) |
| `preprocessing` | `n_top_hvgs`, `hvg_flavor` | 4000, `seurat_v3` | seurat_v3 needs raw counts (we keep them in `layers['counts']`) |
| `integration` | `method` | `harmony` | `none` / `bbknn` / `scvi` |
| `sketch` | `enabled`, `fraction`, `method` | `true`, 0.15, `leverage_score` | matches Oh et al. 2025 |
| `clustering.l1` | `n_pcs`, `resolution` | 20, 0.8 | L1 broad cell types |
| `clustering.l2` | `n_pcs`, `resolution` | 25, 0.1 | L2 subtypes |
| `annotation` | `method`, `hierarchical`, `l2_min_cells` | `markers`, `true`, 200 | `manual` / `celltypist` also supported |
| `differential_expression` | `pseudobulk` | `true` | mice as replicates |
| `segmentation` | `enabled` | `false` | requires per-sample `he_image_path:` |

---

## Run

End-to-end:

```bash
visium-brain --config config/config.yaml all
# or
python scripts/run_all.py --config config/config.yaml
```

Stage by stage (each reads the previous stage's `.h5ad` from `results/`):

```bash
visium-brain qc               # 1. load + QC
visium-brain preprocess       # 2. normalize, HVG, PCA, integrate
visium-brain cluster          # 3. sketch -> L1 cluster + annotate -> L2
visium-brain spatial          # 4. spatial graph, Moran's I, neighborhood enrichment
visium-brain de               # 5. cluster markers + condition + pseudobulk DE
```

Auxiliary commands:

```bash
visium-brain export-markers   # dump cluster markers + cluster_labels.yaml stub
visium-brain segment          # bin2cell nuclear segmentation on 2 µm bins
```

---

## Pipeline overview

| # | Stage | Module(s) | Primary outputs |
|---|---|---|---|
| 1 | Load + QC | `io`, `qc` | `01_qc/adata_qc.h5ad`, `qc_violins.png`, `qc_summary.csv` |
| 2 | Normalize / HVG / PCA / integrate | `preprocessing`, `integration` | `03_integration/adata_integrated.h5ad` |
| 3 | Cluster + annotate (sketch + hierarchical) | `sketch`, `clustering`, `annotation`, `hierarchical` | `05_annotation/adata_annotated.h5ad`, UMAP, spatial maps for L1 + L2 |
| 4 | Spatial analyses | `spatial` | `06_spatial/morans_i.csv`, `nhood_enrichment.png` |
| 5 | Differential expression | `differential` | `07_differential/cluster_markers.csv`, `condition_de_binlevel.csv`, `pseudobulk_de.csv` |

```
results/
├── 01_qc/                     adata_qc.h5ad,        qc_violins.png, qc_summary.csv
├── 03_integration/            adata_integrated.h5ad
├── 04_clusters/               umap_overview.png,    spatial_leiden_l1_*.png
├── 05_annotation/             adata_annotated.h5ad, spatial_cell_type_*.png,
│                              markers_leiden_l1.csv, cluster_labels.yaml,
│                              l2_marker_scores_<L1>.csv
├── 06_spatial/                adata_spatial.h5ad,   morans_i.csv, nhood_enrichment.png
├── 07_differential/           adata_final.h5ad,     *_de.csv, de_top5_dotplot.png
└── segmentation/              <sample_id>/adata_cells.h5ad,  adata_cells_merged.h5ad  (only if enabled)
```

---

## Annotation

### Automated marker scoring (default)

`annotation.method: markers` scores ~18 canonical mouse-brain marker
panels (broad cell types + region/layer markers) with
`sc.tl.score_genes` and assigns each Leiden cluster the panel with the
highest mean score. Good first pass; not a substitute for manual review.

### Hierarchical (L1 → L2)

When `annotation.hierarchical: true` (default), each L1 group is
reclustered at lower resolution (default 25 PCs, res 0.1) and rescored
panel-by-panel. Outputs:

- `obs['cell_type_l1']`, `obs['leiden_l1']`
- `obs['cell_type_l2']`, `obs['leiden_l2']` (composed as `<L1>__<best_panel>_c<leiden>`)
- `uns['l2_marker_scores'][<L1>]` per-cluster score tables (also dumped to CSV)

Disable with `annotation.hierarchical: false`.

### Manual (publication-grade)

```bash
# 1. Run through clustering with method=markers (auto first pass)
visium-brain cluster

# 2. Dump markers + cluster_labels.yaml stub
visium-brain export-markers
# -> results/05_annotation/markers_leiden_l1.csv
# -> results/05_annotation/cluster_labels.yaml      (cluster -> "TBD")

# 3. Open notebooks/manual_annotation.ipynb
#    - rank_genes_groups plots, UMAP, per-sample spatial maps
#    - edit cluster_labels.yaml with cell-type names
#    - save next to the notebook

# 4. Re-run with the manual mapping:
#    annotation.method: manual
#    annotation.manual_labels_l1: results/05_annotation/cluster_labels.yaml
visium-brain cluster
```

L2 reclustering then runs on top of the manual L1 labels.

### CellTypist

`annotation.method: celltypist` runs the configured pretrained model
(default `Mouse_Whole_Brain.pkl`) per bin with majority voting.

---

## Nuclear segmentation (2 µm + H&E)

When 8 µm bins under-resolve small cell types (microglia, endothelial,
immune), enable `segmentation` to aggregate 2 µm bins under
StarDist-segmented nuclei via `bin2cell`. Requires the
`segmentation` extra and a path to each sample's full-resolution H&E
TIFF (Space Ranger does not include this in `binned_outputs/`).

```yaml
segmentation:
  enabled: true
  he_mpp: 0.5
  stardist_model: 2D_versatile_he

samples:
  - sample_id: ctrl_m1
    path: data/raw/ctrl_m1
    he_image_path: data/raw/ctrl_m1/he_full.tif
    ...
```

```bash
visium-brain segment
```

Outputs: `results/segmentation/<sample>/adata_cells.h5ad` per sample
plus a merged `adata_cells_merged.h5ad`. Run separately from the
8 µm pipeline; downstream stages 4–7 are not yet wired to consume the
cell-level AnnData.

---

## Notes on Visium HD specifics

- **Bin size.** Start with **8 µm**: 16× more UMIs per bin than 2 µm,
  small enough for cell-type-scale annotation, large enough for
  third-party tools. For sub-bin-resolution work on small cell types,
  enable `segmentation` to fold 2 µm bins into cells.
- **Replication.** With 2 mice × 2 conditions, bin-level DE is heavily
  pseudoreplicated. The `pseudobulk_de` step aggregates raw counts per
  (sample × cluster) and tests at the mouse level. With n=2 per group
  treat p-values as exploratory and lean on effect sizes; for
  publication-grade statistics swap in `pyDESeq2` or `edgeR` (via
  `rpy2`).
- **Integration.** Defaults to Harmony on `sample_id` because the
  default design has 4 mice across 2 conditions; the original *Nat Genet*
  paper analyses one source and skips integration. Set
  `integration.method: scvi` for stronger correction or `none` to
  match the paper exactly.
- **Random seed.** Every stochastic step (sketch, neighbors, UMAP,
  Leiden, scVI) reads `project.random_seed`. Reproducible by default.

---

## Module reference

| Module | Purpose |
|---|---|
| `visium_brain.io` | Read Space Ranger HD per-sample, concatenate with per-sample images |
| `visium_brain.qc` | Mt/Hb/Ribo flagging, QC metrics, per-sample summary, filters |
| `visium_brain.preprocessing` | Log-CPM, HVGs (`seurat_v3` default), scaling, PCA |
| `visium_brain.integration` | Harmony / BBKNN / scVI / none |
| `visium_brain.sketch` | Leverage scores, sketch sampling, kNN label propagation |
| `visium_brain.clustering` | Per-level neighbors / UMAP / Leiden |
| `visium_brain.annotation` | Marker-score panels, CellTypist, manual-YAML labels, marker-table export |
| `visium_brain.hierarchical` | L1 → L2 reclustering and sub-annotation |
| `visium_brain.spatial` | Per-sample spatial graph, Moran's I, neighborhood enrichment |
| `visium_brain.differential` | Cluster markers, condition DE, pseudobulk DE |
| `visium_brain.segmentation` | `bin2cell` 2 µm + H&E nuclear segmentation |
| `visium_brain.plotting` | QC violins, UMAP, per-sample spatial maps, DE dotplots |
| `visium_brain.pipeline` | Stage orchestration; entry points used by the CLI |
| `visium_brain.cli` | `visium-brain` Click CLI |
| `visium_brain.utils` | Config loading, logging, IO paths, `Sample` dataclass |

---

## Testing

```bash
pytest -q
```

The smoke tests build a synthetic Visium-shaped AnnData (4 samples × 2
conditions) and exercise QC, preprocessing, leverage-score sketch +
kNN propagation accuracy, hierarchical L2 reclustering, manual-label
application, and pseudobulk DE — all without needing real Space
Ranger output.

---

## References

- 10x Genomics. *Visium HD Spatial Gene Expression*.
  https://www.10xgenomics.com/support/spatial-gene-expression-hd
- Oh et al., *Nat Genet* 2025. Visium HD characterization workflow this
  pipeline mirrors. https://www.nature.com/articles/s41588-025-02193-3
- Wolf et al., *Genome Biol* 2018. **Scanpy.**
- Virshup et al. **AnnData.**
- Palla et al., *Nat Methods* 2022. **Squidpy.**
- Korsunsky et al., *Nat Methods* 2019. **Harmony.**
- Hao et al., *Nat Biotechnol* 2024. **Seurat v5** sketch / leverage-score
  sampling that the L1 step here mirrors.
- Polanski et al., 2024. **bin2cell** for Visium HD nuclear segmentation.
- Domínguez Conde et al., *Science* 2022. **CellTypist.**
- Lopez et al., *Nat Methods* 2018. **scVI.**
- Schmidt et al., *Nucleic Acids Res* 2018; Stardist:
  Schmidt et al. 2018. **StarDist.**

---

## License

MIT.
