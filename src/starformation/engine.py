"""Rendering-free, fixed-step N-body experiment orchestration."""

from dataclasses import dataclass

import numpy as np

from .cloud import accuracy_history, cloud_geometry, snapshot
from .config import SimulationConfig
from .gravity import compute_acceleration
from .initial_conditions import InitialState, make_initial_state
from .integration import leapfrog_step


@dataclass(frozen=True)
class SimulationResult:
    config: SimulationConfig
    initial: InitialState
    trajectory: dict
    report: dict


def run(config: SimulationConfig, initial: InitialState | None = None):
    """Save t=0 and all completed steps; no plots, global state or output writes."""
    initial = make_initial_state(config) if initial is None else initial
    pos, vel, mass = initial.pos.copy(), initial.vel.copy(), initial.mass.copy()
    if pos.shape != (config.particle_count, 3) or vel.shape != pos.shape or mass.shape != (len(pos),):
        raise ValueError('Initial state and configuration have incompatible shapes')
    half = cloud_geometry(pos, mass)['r50_m']
    acc = compute_acceleration(pos, mass, config.softening_m, G=config.gravitational_constant_si)
    positions, velocities, summaries, densities, supports, radial = [], [], [], [], [], []
    for step in range(config.steps+1):
        measurement, particles = snapshot(
            pos, vel, mass, G=config.gravitational_constant_si, softening=config.softening_m,
            sound_speed=config.sound_speed_m_s, k=config.neighbours,
            reference_half_radius=half, aperture_radius=config.radius_m/2)
        positions.append(pos.copy())
        velocities.append(vel.copy())
        summaries.append(measurement)
        densities.append(particles['local_density_kg_m3'])
        supports.append(particles['neighbour_radius_m'])
        radial.append(particles['radial_velocity_m_s'])
        if step < config.steps:
            pos, vel, acc = leapfrog_step(pos, vel, acc, config.dt_s, mass,
                                         config.softening_m, G=config.gravitational_constant_si)
    times = np.arange(config.steps+1)*config.dt_s
    accuracy = accuracy_history(summaries, times, mass_total=mass.sum(),
                                reference_radius=config.radius_m, G=config.gravitational_constant_si,
                                energy_budget=config.energy_budget)
    trajectory = {'pos': np.asarray(positions), 'vel': np.asarray(velocities),
                  'mass': mass.copy(), 'time_s': times,
                  'density_kg_m3': np.asarray(densities),
                  'neighbour_radius_m': np.asarray(supports),
                  'radial_velocity_m_s': np.asarray(radial)}
    report = {'schema_version': 1, 'units': 'SI',
              'configuration': {'dt_s': config.dt_s, 'steps': config.steps,
                                'softening_m': config.softening_m,
                                'G_si': config.gravitational_constant_si,
                                'initial_radius_m': config.radius_m,
                                'sound_speed_m_s': config.sound_speed_m_s,
                                'seed': initial.seed, 'particles': len(mass),
                                'k': config.neighbours, 'initial_half_radius_m': half,
                                'aperture_radius_m': config.radius_m/2,
                                'energy_budget': config.energy_budget},
              'definitions': {'dynamics': 'Plummer-softened direct pair gravity; fixed-step KDK',
                              'density': '(k-1)/k neighbour proxy, no edge correction',
                              'energy': 'K+U with pairwise Plummer potential'},
              'accuracy': accuracy, 'snapshots': summaries}
    return SimulationResult(config, initial, trajectory, report)
