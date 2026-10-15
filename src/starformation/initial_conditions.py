"""Legacy MT19937 state construction, isolated from NumPy's global RNG."""

from dataclasses import dataclass

import numpy as np

from .config import SimulationConfig


@dataclass(frozen=True)
class InitialState:
    pos: np.ndarray
    vel: np.ndarray
    mass: np.ndarray
    method: str
    seed: int | None


def make_initial_state(config: SimulationConfig, *, pos=None, vel=None, mass=None):
    supplied = (pos is not None, vel is not None, mass is not None)
    if any(supplied) and not all(supplied):
        raise ValueError('Explicit initial state requires pos, vel and mass together')
    n = config.particle_count
    if all(supplied):
        p, v, m = (np.asarray(x, dtype=float).copy() for x in (pos, vel, mass))
        method, seed = 'explicit_arrays', None
    else:
        rng = np.random.RandomState(config.seed)
        R = config.radius_m
        r_sphere = R * np.cbrt(rng.rand(n))
        theta = np.arccos(1 - 2*rng.rand(n))
        phi = 2 * np.pi * rng.rand(n)
        p = np.column_stack([r_sphere * np.sin(theta) * np.cos(phi),
                             r_sphere * np.sin(theta) * np.sin(phi),
                             r_sphere * np.cos(theta)])
        virial = np.sqrt(config.gravitational_constant_si * config.cloud_mass_kg / R)
        v = (rng.randn(n, 3) * config.velocity_sigma_over_virial * virial +
             np.column_stack([-p[:, 1] * config.rotation_omega_over_dynamical * virial / R,
                               p[:, 0] * config.rotation_omega_over_dynamical * virial / R,
                               np.zeros(n)]))
        m = np.full(n, config.cloud_mass_kg/n)
        method, seed = 'legacy_uniform_sphere_mt19937', config.seed
    if p.shape != (n, 3) or v.shape != p.shape or m.shape != (n,):
        raise ValueError('Initial state shapes must be (N,3), (N,3), (N,)')
    if not all(np.isfinite(x).all() for x in (p, v, m)) or (m <= 0).any():
        raise ValueError('Initial state must be finite with strictly positive masses')
    if not np.isclose(m.sum(), config.cloud_mass_kg, rtol=1e-12, atol=0):
        raise ValueError('Initial particle masses must sum to configured cloud mass')
    if config.recenter_com:
        p -= np.sum(m[:, None]*p, axis=0)/m.sum()
        method += '+recenter_com'
    if config.recenter_velocity:
        v -= np.sum(m[:, None]*v, axis=0)/m.sum()
        method += '+recenter_velocity'
    return InitialState(p, v, m, method, seed)
