"""Compatibility import for the validated Phase 1 N-body API."""

import sys
from pathlib import Path

_src = str(Path(__file__).resolve().parent / "src")
if _src not in sys.path:
    sys.path.insert(0, _src)

from starformation.gravity import compute_acceleration, potential_energy
from starformation.integration import leapfrog_step
from starformation.energy import diagnostics

__all__ = ["compute_acceleration", "potential_energy", "leapfrog_step", "diagnostics"]
