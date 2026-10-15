"""Validated fixed-step kick-drift-kick integrator."""

from ._nbody_impl import leapfrog_step

__all__ = ["leapfrog_step"]
