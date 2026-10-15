"""Reproduce a short Phase 1 cloud trajectory and save Phase 2 measurements.

An explicit timestep choice is required. End time is fixed at the tested 1.2 t_ff.
No plotting here. No legacy script import. Existing outputs are never overwritten.
"""

import argparse
import csv
import hashlib
import importlib.metadata
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'tests'))
from baseline_support import load_baseline
from nbody import compute_acceleration, leapfrog_step
from scientific_diagnostics import accuracy_history, cloud_geometry, save_diagnostics, snapshot


def generate(output, dt_fraction, k=16):
    if dt_fraction not in (.012, .006, .003):
        raise ValueError('Use a Phase 1 tested timestep: .012, .006 or .003')
    n = load_baseline()
    if isinstance(k, bool) or not isinstance(k, int) or not 2 <= k < n['N_PARTICLES']:
        raise ValueError('k must satisfy 2 <= k < N')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    pos, vel, mass = n['pos'].copy(), n['vel'].copy(), n['mass'].copy()
    # Confirm the corrected seed-42 state against the committed Phase 1 artifact.
    with np.load(ROOT/'validation/phase1/cloud_initial.npz', allow_pickle=False) as initial:
        for actual, key in ((pos, 'pos'), (vel, 'vel'), (mass, 'mass')):
            np.testing.assert_array_equal(actual, initial[key])
    dt = dt_fraction*n['T']
    steps = round(1.2/dt_fraction)
    initial_half = cloud_geometry(pos, mass)['r50_m']
    aperture = n['R_cloud']/2
    summaries, positions, velocities, densities, supports, radial = [], [], [], [], [], []
    acc = compute_acceleration(pos, mass, n['SOFTENING'], G=n['G_si'])
    for step in range(steps+1):
        measurement, particles = snapshot(pos, vel, mass, G=n['G_si'], softening=n['SOFTENING'],
                                         sound_speed=n['c_sound'], k=k,
                                         reference_half_radius=initial_half, aperture_radius=aperture)
        summaries.append(measurement)
        positions.append(pos.copy())
        velocities.append(vel.copy())
        densities.append(particles['local_density_kg_m3'])
        supports.append(particles['neighbour_radius_m'])
        radial.append(particles['radial_velocity_m_s'])
        if step != steps:
            pos, vel, acc = leapfrog_step(pos, vel, acc, dt, mass, n['SOFTENING'], G=n['G_si'])
    times = np.arange(steps+1)*dt
    accuracy = accuracy_history(summaries, times, mass_total=mass.sum(), reference_radius=n['R_cloud'], G=n['G_si'])
    with (ROOT/f'validation/phase1/cloud_dt{dt_fraction:g}.csv').open() as stream:
        reference = list(csv.DictReader(stream))
    reference_energy = np.array([float(row['total_energy']) for row in reference])
    np.testing.assert_allclose([s['energy'] for s in summaries], reference_energy, rtol=1e-12, atol=0)
    with np.load(ROOT/f'validation/phase1/cloud_final_dt{dt_fraction:g}.npz', allow_pickle=False) as final:
        np.testing.assert_allclose(pos, final['pos'], rtol=1e-12, atol=0)
        np.testing.assert_allclose(vel, final['vel'], rtol=1e-12, atol=0)
    np.savez_compressed(output/'trajectory.npz', pos=positions, vel=velocities, mass=mass,
                        time_s=times, density_kg_m3=densities, neighbour_radius_m=supports,
                        radial_velocity_m_s=radial)
    sources = ['scientific_diagnostics.py', 'experiments/phase2_preview.py', 'nbody.py',
               'tests/baseline_support.py', 'baseline/reference.json', 'validation/phase1/results.json',
               f'validation/phase1/cloud_dt{dt_fraction:g}.csv']
    payload = {
        'schema_version': 1, 'units': 'SI',
        'configuration': {'dt_fraction': dt_fraction, 'dt_s': dt, 'steps': steps, 'end_tff': 1.2,
                          'initial_tff_s': n['T'], 'softening_m': n['SOFTENING'], 'G_si': n['G_si'],
                          'initial_radius_m': n['R_cloud'], 'sound_speed_m_s': n['c_sound'], 'seed': 42,
                          'particles': len(mass), 'k': k, 'initial_half_radius_m': initial_half,
                          'aperture_radius_m': aperture, 'energy_budget': .01},
        'definitions': {'density': '(k-1)/k times k-neighbour mass / sphere volume; exclude self; no edge correction',
                        'mean_density': 'actual enclosed mass / sphere volume at COM-centred r50',
                        'jeans_length_comparison': 'L=2*r50 compared with lambda_J using mean_density',
                        'virial_U': '2*K_internal/abs(U)', 'virial_W': '2*K_internal/abs(sum m (r-Rcm).a)',
                        'contraction': '1-r50(t)/r50(0), geometric, unbounded below',
                        'aperture': 'fixed physical radius R_initial/2, centred on measured COM'},
        'accuracy': accuracy, 'snapshots': summaries,
        'trajectory_sha256': hashlib.sha256((output/'trajectory.npz').read_bytes()).hexdigest(),
        'provenance': {'git_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                       'python': sys.version, 'versions': {p: importlib.metadata.version(p) for p in ('numpy', 'astropy', 'matplotlib')},
                       'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
                       'phase1_energy_history_max_relative_difference': float(np.max(abs(
                           (np.array([s['energy'] for s in summaries])-reference_energy)/reference_energy))),
                       'phase1_final_state_regression': 'PASS at rtol=1e-12, atol=0'}
    }
    save_diagnostics(output/'diagnostics.json', payload)
    print(f'{output}: {accuracy["energy_status"]}; peak error {accuracy["peak_energy_error"]:.6%}; '
          f'{steps+1} saved states; Phase 1 history/final-state regression passed.')
    return accuracy


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--dt-fraction', type=float, choices=(.012, .006, .003), required=True)
    parser.add_argument('--neighbours', type=int, default=16)
    args = parser.parse_args()
    result = generate(args.output, args.dt_fraction, args.neighbours)
    raise SystemExit(1 if result['energy_status'] == 'FAIL' else 0)
