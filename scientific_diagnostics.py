"""Compatibility import for the Phase 2 scientific diagnostics API."""

import sys
from pathlib import Path

_src = str(Path(__file__).resolve().parent / "src")
if _src not in sys.path:
    sys.path.insert(0, _src)

from starformation.cloud import (accuracy_history, centre_of_mass, cloud_geometry,
                                 contraction, enclosed_mass, enclosed_mass_radii,
                                 jeans_scales, json_ready, local_particle_density,
                                 radial_flow, save_diagnostics, snapshot,
                                 virial_measurements)

__all__ = ["accuracy_history", "centre_of_mass", "cloud_geometry", "contraction",
           "enclosed_mass", "enclosed_mass_radii", "jeans_scales", "json_ready",
           "local_particle_density", "radial_flow", "save_diagnostics", "snapshot",
           "virial_measurements"]
