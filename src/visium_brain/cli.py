"""Command-line entrypoint: ``visium-brain``."""

from __future__ import annotations

import click

from . import pipeline
from .utils import load_config, set_seed, setup_logging


@click.group()
@click.option("--config", "-c", default="config/config.yaml", show_default=True,
              help="Path to YAML config file.")
@click.option("--log-level", default="INFO", show_default=True)
@click.pass_context
def main(ctx: click.Context, config: str, log_level: str) -> None:
    """Visium HD mouse-brain analysis pipeline."""
    setup_logging(log_level)
    cfg = load_config(config)
    set_seed(cfg["project"].get("random_seed", 0))
    ctx.obj = cfg


@main.command("qc")
@click.pass_obj
def qc_cmd(cfg):
    """Load samples, compute QC, and write a filtered AnnData."""
    pipeline.run_qc(cfg)


@main.command("preprocess")
@click.pass_obj
def preprocess_cmd(cfg):
    """Normalize, HVG, scale, PCA, and integrate samples."""
    pipeline.run_preprocess_integrate(cfg)


@main.command("cluster")
@click.pass_obj
def cluster_cmd(cfg):
    """Compute neighbors / UMAP / Leiden and annotate cell types."""
    pipeline.run_cluster_annotate(cfg)


@main.command("spatial")
@click.pass_obj
def spatial_cmd(cfg):
    """Run spatial-aware analyses (Moran's I, neighborhood enrichment)."""
    pipeline.run_spatial(cfg)


@main.command("de")
@click.pass_obj
def de_cmd(cfg):
    """Run differential expression (cluster markers + condition DE + pseudobulk)."""
    pipeline.run_differential(cfg)


@main.command("all")
@click.pass_context
def all_cmd(ctx):
    """Run the entire pipeline end-to-end."""
    pipeline.run_all(ctx.parent.params["config"])


if __name__ == "__main__":  # pragma: no cover
    main()
