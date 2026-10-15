"""Strict version-1 JSON configuration. Public input units are in field names."""

from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np

from .units import free_fall_time, legacy_si_constants

SCHEMA_VERSION = 1


def legacy_input():
    G, _, _ = legacy_si_constants()
    return {"schema_version": 1, "cloud_mass_msun": 1.0, "radius_pc": 0.1,
            "particle_count": 300, "seed": 42, "gravitational_constant_si": G,
            "softening_over_radius": 0.05, "dt_over_initial_tff": 0.012,
            "steps": 1080, "velocity_sigma_over_virial": 0.05,
            "rotation_omega_over_dynamical": 0.03, "recenter_com": False,
            "recenter_velocity": False, "sound_speed_m_s": 200.0,
            "neighbours": 16, "energy_budget": 0.01,
            "output": {"simulation_path": None, "render_path": None},
            "render": {"view": "com", "animation": False, "fps": 5, "dpi": 120}}


def _number(x, name, *, positive=False, nonnegative=False):
    if isinstance(x, bool) or not isinstance(x, (int, float)) or not np.isfinite(x):
        raise ValueError(f"{name} must be a finite number")
    if positive and x <= 0 or nonnegative and x < 0:
        raise ValueError(f"{name} has invalid sign")
    return float(x)


def _integer(x, name, minimum):
    if isinstance(x, bool) or not isinstance(x, int) or x < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return x


@dataclass(frozen=True)
class SimulationConfig:
    cloud_mass_kg: float
    radius_m: float
    particle_count: int
    seed: int
    gravitational_constant_si: float
    softening_m: float
    dt_s: float
    steps: int
    velocity_sigma_over_virial: float
    rotation_omega_over_dynamical: float
    recenter_com: bool
    recenter_velocity: bool
    sound_speed_m_s: float
    neighbours: int
    energy_budget: float
    simulation_path: str | None
    render_path: str | None
    render_view: str
    render_animation: bool
    render_fps: int
    render_dpi: int
    initial_tff_s: float

    def __post_init__(self):
        for name in ('cloud_mass_kg', 'radius_m', 'gravitational_constant_si',
                     'dt_s', 'sound_speed_m_s', 'energy_budget', 'initial_tff_s'):
            _number(getattr(self, name), name, positive=True)
        for name in ('softening_m', 'velocity_sigma_over_virial',
                     'rotation_omega_over_dynamical'):
            _number(getattr(self, name), name, nonnegative=True)
        _integer(self.particle_count, 'particle_count', 3)
        _integer(self.seed, 'seed', 0)
        _integer(self.steps, 'steps', 1)
        _integer(self.neighbours, 'neighbours', 2)
        _integer(self.render_fps, 'render_fps', 1)
        _integer(self.render_dpi, 'render_dpi', 1)
        if self.seed >= 2**32 or self.neighbours >= self.particle_count:
            raise ValueError('Invalid seed or neighbour count')
        if not isinstance(self.recenter_com, bool) or not isinstance(self.recenter_velocity, bool):
            raise ValueError('COM corrections must be Boolean')
        if not isinstance(self.render_animation, bool) or self.render_view not in ('com', 'fixed'):
            raise ValueError('Invalid render settings')
        for name in ('simulation_path', 'render_path'):
            value = getattr(self, name)
            if value is not None and (not isinstance(value, str) or not value):
                raise ValueError(f'{name} must be a nonempty path or None')
        expected_tff = free_fall_time(self.cloud_mass_kg, self.radius_m,
                                      self.gravitational_constant_si)
        if not np.isclose(self.initial_tff_s, expected_tff, rtol=1e-12, atol=0):
            raise ValueError('initial_tff_s is inconsistent with mass, radius and G')

    def as_si_dict(self):
        return asdict(self)


def from_input(data):
    if not isinstance(data, dict) or set(data) != set(legacy_input()):
        raise ValueError(f"Configuration keys must be exactly {sorted(legacy_input())}")
    if data['schema_version'] != SCHEMA_VERSION or isinstance(data['schema_version'], bool):
        raise ValueError('Unsupported configuration schema_version')
    output, render = data['output'], data['render']
    if not isinstance(output, dict) or set(output) != {'simulation_path', 'render_path'}:
        raise ValueError('output requires simulation_path and render_path')
    if not isinstance(render, dict) or set(render) != {'view', 'animation', 'fps', 'dpi'}:
        raise ValueError('render requires view, animation, fps and dpi')
    for key in output:
        if output[key] is not None and (not isinstance(output[key], str) or not output[key]):
            raise ValueError(f'output.{key} must be a nonempty path or null')
    if render['view'] not in ('com', 'fixed') or not isinstance(render['animation'], bool):
        raise ValueError('render view/animation invalid')
    G0, solar_kg, pc_m = legacy_si_constants()
    del G0
    mass = _number(data['cloud_mass_msun'], 'cloud_mass_msun', positive=True) * solar_kg
    radius = _number(data['radius_pc'], 'radius_pc', positive=True) * pc_m
    G = _number(data['gravitational_constant_si'], 'gravitational_constant_si', positive=True)
    n = _integer(data['particle_count'], 'particle_count', 3)
    k = _integer(data['neighbours'], 'neighbours', 2)
    if k >= n:
        raise ValueError('neighbours must be less than particle_count')
    seed = _integer(data['seed'], 'seed', 0)
    if seed >= 2**32:
        raise ValueError('seed must fit NumPy MT19937 uint32 range')
    steps = _integer(data['steps'], 'steps', 1)
    tff = free_fall_time(mass, radius, G)
    for name in ('recenter_com', 'recenter_velocity'):
        if not isinstance(data[name], bool):
            raise ValueError(f'{name} must be Boolean')
    return SimulationConfig(
        mass, radius, n, seed, G,
        radius * _number(data['softening_over_radius'], 'softening_over_radius', nonnegative=True),
        tff * _number(data['dt_over_initial_tff'], 'dt_over_initial_tff', positive=True),
        steps,
        _number(data['velocity_sigma_over_virial'], 'velocity_sigma_over_virial', nonnegative=True),
        _number(data['rotation_omega_over_dynamical'], 'rotation_omega_over_dynamical', nonnegative=True),
        data['recenter_com'], data['recenter_velocity'],
        _number(data['sound_speed_m_s'], 'sound_speed_m_s', positive=True), k,
        _number(data['energy_budget'], 'energy_budget', positive=True),
        output['simulation_path'], output['render_path'], render['view'], render['animation'],
        _integer(render['fps'], 'render.fps', 1), _integer(render['dpi'], 'render.dpi', 1), tff)


def load_config(path):
    return from_input(json.loads(Path(path).read_text()))
