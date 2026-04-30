"""Analysis pipeline for 10x Genomics Visium HD mouse brain data."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("visium-brain")
except PackageNotFoundError:  # pragma: no cover - package not installed
    __version__ = "0.1.0"

__all__ = ["__version__"]
