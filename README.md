# Visium-Brain

A Python pipeline for analyzing 10x Genomics **Visium HD** spatial gene
expression data on mouse brain sections. The default configuration assumes
a 2 × 2 design — **2 conditions × 2 mice per condition**, one section per
mouse — but any number of samples can be added via `config/config.yaml`.

The pipeline is built on `scanpy` + `anndata` + `squidpy`, with optional
`spatialdata` / `harmonypy` / `celltypist` integrations.

## Pipeline stages

| # | Stage | Module | Outputs |
|---|---|---|---|
| 1 | Load + QC | `visium_brain.io`, `visium_brain.qc` | `results/01_qc/adata_qc.h5ad`, QC violins, per-sample summary |
| 2 | Normalize / HVG / PCA / integrate | `preprocessing`, `integration` | `results/03_integration/adata_integrated.h5ad` |
| 3 | Cluster + annotate (sketch + hierarchical) | `sketch`, `clustering`, `annotation`, `hierarchical` | UMAP, spatial maps of L1 and L2 cell types |
| 4 | Spatial analyses | `spatial` | Moran's I table, neighborhood-enrichment heatmap |
| 5 | Differential expression | `differential` | Cluster markers, condition DE (bin-level), pseudobulk DE per (sample × cluster) |

Defaults: `square_008um` bins, Harmony integration on `sample_id`,
**leverage-score sketch (15%) + kNN label propagation** for L1 clustering,
**hierarchical L1 → L2 annotation** (recluster each L1 cell type at lower
resolution), marker-score cell-type annotation, Wilcoxon DE plus pseudobulk
DE using mice as the unit of replication.

### Sketch + hierarchical annotation

The L1 step (clustering + annotation) runs on a leverage-score sketch
of the data — leverage scores from PCA oversample rare populations, so
the sketch preserves biological complexity better than uniform
sampling. Labels are then propagated to all bins via a kNN classifier
in the integrated PCA / Harmony space. This mirrors Seurat v5's
`SketchData` / `ProjectData` and the methods in Oh et al., *Nat Genet*
2025 (https://www.nature.com/articles/s41588-025-02193-3).

After L1, each L1 group is reclustered at lower resolution (default
25 PCs, resolution 0.1) and rescored against the marker panels;
results are written to `obs['cell_type_l1']`, `obs['cell_type_l2']`,
and `obs['leiden_l1']` / `obs['leiden_l2']`. Disable either with
`sketch.enabled: false` or `annotation.hierarchical: false` in the
config.

## Layout

```
config/config.yaml          # samples + all pipeline parameters
src/visium_brain/           # importable package (modules per stage)
scripts/0[1-5]_*.py         # thin runners per stage
scripts/run_all.py          # full end-to-end run
tests/test_smoke.py         # synthetic-data smoke tests
results/                    # all outputs (gitignored)
data/raw/                   # Space Ranger outputs (gitignored)
```

## Install

```bash
# conda (recommended for igraph / opencv pinning)
conda env create -f environment.yml
conda activate visium-brain
pip install -e .

# or pure pip
python -m venv .venv && source .venv/bin/activate
pip install -e .
```

## Configure your samples

Edit `config/config.yaml` so each entry under `samples:` points at a
Space Ranger output directory containing `binned_outputs/square_008um/`:

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

## Run

End-to-end:
```bash
python scripts/run_all.py --config config/config.yaml
# or
visium-brain --config config/config.yaml all
```

Stage by stage (each stage reads the previous stage's `.h5ad`):
```bash
visium-brain qc
visium-brain preprocess
visium-brain cluster
visium-brain spatial
visium-brain de
```

## Tests

```bash
pytest -q
```
The smoke tests build a synthetic Visium-shaped AnnData and exercise QC,
preprocessing, and pseudobulk DE without needing real Space Ranger output.

## Notes on Visium HD specifics

* **Bin size.** Space Ranger writes 2 µm, 8 µm, and 16 µm binned outputs.
  Start with **8 µm**; it balances resolution with per-bin counts. For
  cell-resolution work, switch to 2 µm and add nucleus-segmentation
  aggregation (e.g. `bin2cell`) — out of scope here.
* **Replication.** With 2 mice × 2 conditions, bin-level DE is heavily
  pseudoreplicated. The `pseudobulk_de` step aggregates raw counts per
  (sample × cluster) and tests at the mouse level. With n=2 per group
  treat p-values as exploratory and rely on effect sizes; for a
  publication-grade test substitute `pyDESeq2` or `edgeR` via `rpy2`.
* **Integration.** Defaults to Harmony on `sample_id`. For stronger batch
  correction across mice + conditions, set `integration.method: scvi` in
  the config (requires `scvi-tools`).
