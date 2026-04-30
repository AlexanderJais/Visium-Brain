"""Utility helpers: config loading, logging, IO paths."""

from __future__ import annotations

import logging
import os
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import yaml


@dataclass
class Sample:
    sample_id: str
    condition: str
    mouse_id: str
    path: Path

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "Sample":
        return cls(
            sample_id=str(d["sample_id"]),
            condition=str(d["condition"]),
            mouse_id=str(d["mouse_id"]),
            path=Path(d["path"]),
        )


def load_config(path: str | Path) -> dict[str, Any]:
    with open(path, "r") as fh:
        cfg = yaml.safe_load(fh)
    cfg["_config_path"] = str(Path(path).resolve())
    return cfg


def get_samples(cfg: dict[str, Any]) -> list[Sample]:
    return [Sample.from_dict(s) for s in cfg["samples"]]


def setup_logging(level: str = "INFO") -> logging.Logger:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s | %(levelname)-7s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    return logging.getLogger("visium_brain")


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def output_paths(cfg: dict[str, Any]) -> dict[str, Path]:
    out = ensure_dir(cfg["project"]["output_dir"])
    figs = ensure_dir(cfg["project"]["figures_dir"])
    return {
        "out": out,
        "figures": figs,
        "qc": ensure_dir(out / "01_qc"),
        "preproc": ensure_dir(out / "02_preprocess"),
        "integration": ensure_dir(out / "03_integration"),
        "clusters": ensure_dir(out / "04_clusters"),
        "annotation": ensure_dir(out / "05_annotation"),
        "spatial": ensure_dir(out / "06_spatial"),
        "de": ensure_dir(out / "07_differential"),
    }
