"""Compatibility import and CLI for the Phase 2 saved-data renderer."""

import argparse
import sys
from pathlib import Path

_src = str(Path(__file__).resolve().parent / "src")
if _src not in sys.path:
    sys.path.insert(0, _src)

from starformation.visualization import (AU_M, YEAR_S, ProjectionView, load_dataset,
                                         render, research_figure, validate_dataset)

__all__ = ["AU_M", "YEAR_S", "ProjectionView", "load_dataset", "render",
           "research_figure", "validate_dataset"]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--view", choices=("com", "fixed"), default="com")
    parser.add_argument("--animation", action="store_true")
    args = parser.parse_args()
    render(args.data, args.output, view=args.view, animation=args.animation)
