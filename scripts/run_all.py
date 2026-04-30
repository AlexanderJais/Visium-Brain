#!/usr/bin/env python
"""Run the full Visium HD pipeline end-to-end."""
from __future__ import annotations

import argparse

from visium_brain import pipeline


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", "-c", default="config/config.yaml")
    args = parser.parse_args()
    pipeline.run_all(args.config)


if __name__ == "__main__":
    main()
