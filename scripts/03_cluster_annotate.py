#!/usr/bin/env python
"""Stage 3: neighbors, UMAP, Leiden, and cell-type annotation."""
from __future__ import annotations

import argparse

from visium_brain import pipeline
from visium_brain.utils import load_config, set_seed, setup_logging


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", "-c", default="config/config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    setup_logging()
    set_seed(cfg["project"].get("random_seed", 0))
    pipeline.run_cluster_annotate(cfg)


if __name__ == "__main__":
    main()
