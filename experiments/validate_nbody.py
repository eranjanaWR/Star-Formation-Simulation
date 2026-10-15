"""Reproducible Phase 1 experiments; never imports the executable legacy module.

Run: .venv/bin/python -B experiments/validate_nbody.py --output NEW_DIRECTORY
Output directories must not exist. No reference artifact can be overwritten.
"""

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
from baseline_support import load_baseline
from nbody import compute_acceleration, diagnostics, leapfrog_step


def circular_reference(t, mass, separation, eps, G):
    """Exact circular orbit in the Plummer two-body potential, barycentric frame."""
    omega = np.sqrt(G * np.sum(mass) / (separation**2 + eps**2)**1.5)
    angle = omega * t
    r = separation * np.array([np.cos(angle), np.sin(angle), 0.])
    v = separation * omega * np.array([-np.sin(angle), np.cos(angle), 0.])
    weights = np.array([-mass[1], mass[0]]) / np.sum(mass)
    return weights[:, None] * r, weights[:, None] * v


def kepler_reference(t, mass, semimajor, eccentricity, G):
    """Unsoftened ellipse from Kepler's equation, independent of KDK integration."""
    motion = np.sqrt(G * np.sum(mass) / semimajor**3)
    mean_anomaly = (motion * t + np.pi) % (2*np.pi) - np.pi
    anomaly = mean_anomaly
    for _ in range(30):
        delta = (anomaly - eccentricity*np.sin(anomaly) - mean_anomaly) / (
            1 - eccentricity*np.cos(anomaly))
        anomaly -= delta
        if abs(delta) < 1e-14:
            break
    else:
        raise RuntimeError("Kepler reference did not converge")
    r = semimajor * np.array([np.cos(anomaly)-eccentricity,
                              np.sqrt(1-eccentricity**2)*np.sin(anomaly), 0.])
    v = semimajor * motion / (1-eccentricity*np.cos(anomaly)) * np.array([
        -np.sin(anomaly), np.sqrt(1-eccentricity**2)*np.cos(anomaly), 0.])
    weights = np.array([-mass[1], mass[0]]) / np.sum(mass)
    return weights[:, None]*r, weights[:, None]*v


def trajectory(p, v, mass, eps, G, dt, steps, length_scale, exact=None):
    """All diagnostics at t=0 and after every synchronized full KDK step."""
    p, v = p.copy(), v.copy()
    a = compute_acceleration(p, mass, eps, G=G)
    initial = diagnostics(p, v, mass, eps, G=G)
    momentum_scale = np.sum(mass) * np.sqrt(G*np.sum(mass)/length_scale)
    history = []
    for step in range(steps + 1):
        t = step * dt
        d = diagnostics(p, v, mass, eps, G=G)
        row = {"step": step, "time": t, "kinetic": d['kinetic'],
               "potential": d['potential'], "total_energy": d['energy'],
               "energy_error": (d['energy']-initial['energy'])/abs(initial['energy']),
               "momentum_error": float(np.linalg.norm(d['momentum']-initial['momentum'])/momentum_scale),
               "angular_momentum_error": float(np.linalg.norm(d['angular_momentum']-initial['angular_momentum'])/
                                               (momentum_scale*length_scale)),
               "com_residual": float(np.linalg.norm(d['com']-initial['com']-t*initial['com_velocity'])/length_scale)}
        for key in ('momentum', 'angular_momentum', 'com', 'com_velocity'):
            row.update({f'{key}_{axis}': float(d[key][i]) for i, axis in enumerate('xyz')})
        if exact is not None:
            expected, _ = exact(t)
            row['relative_position_error'] = float(np.linalg.norm((p[1]-p[0]) -
                                                       (expected[1]-expected[0]))/length_scale)
        if not all(np.isfinite(value) for value in row.values()):
            raise FloatingPointError(f'Nonfinite diagnostic at step {step}')
        history.append(row)
        if step != steps:
            p, v, a = leapfrog_step(p, v, a, dt, mass, eps, G=G)
    metrics = {f'max_{key}': max(abs(row[key]) for row in history) for key in (
        'energy_error', 'momentum_error', 'angular_momentum_error', 'com_residual')}
    metrics['final_energy_error'] = history[-1]['energy_error']
    if exact is not None:
        metrics['max_relative_position_error'] = max(row['relative_position_error'] for row in history)
        metrics['final_relative_position_error'] = history[-1]['relative_position_error']
    return history, metrics, p, v


def write_csv(path, rows):
    with path.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run_orbits(config, output=None):
    c = config['orbits']
    mass = np.array(c['masses'])
    G, radius = c['G'], c['separation_or_semimajor_axis']
    rows = []
    for kind, eps in [('circular', e) for e in c['circular_softening']] + [('eccentric', 0.)]:
        if kind == 'circular':
            period = 2*np.pi * (radius**2+eps**2)**.75 / np.sqrt(G*mass.sum())
            exact = lambda t: circular_reference(t, mass, radius, eps, G)
        else:
            period = 2*np.pi*np.sqrt(radius**3/(G*mass.sum()))
            exact = lambda t: kepler_reference(t, mass, radius, c['kepler_eccentricity'], G)
        previous_error = None
        for count in c['steps_per_period']:
            p, v = exact(0.)
            history, metrics, _, _ = trajectory(p, v, mass, eps, G, period/count,
                                                c['periods']*count, radius, exact)
            error = metrics['max_relative_position_error']
            order = None if previous_error is None else math.log2(previous_error/error)
            previous_error = error
            rows.append({'kind': kind, 'softening': eps, 'steps_per_period': count,
                         'period': float(period), 'dt': float(period/count),
                         'observed_order': order, **metrics})
            if output is not None and count == max(c['steps_per_period']):
                write_csv(output/f'orbit_{kind}_eps{eps:g}.csv', history)
    return rows


def plot_results(output, cloud_histories, orbit_rows, softening_rows, tff):
    # Only explicit report generation imports plotting; solver/tests stay headless.
    with tempfile.TemporaryDirectory(prefix='nbody-mpl-') as cache:
        old = os.environ.get('MPLCONFIGDIR')
        os.environ['MPLCONFIGDIR'] = cache
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
            for key, ax, title in zip(('energy_error', 'momentum_error', 'angular_momentum_error', 'com_residual'),
                                      axes.flat, ('Signed energy error / |E0|', 'Momentum drift / M V',
                                                  'Angular momentum drift / M R V', 'COM residual / R')):
                for frac, history in cloud_histories:
                    ax.plot([r['time']/tff for r in history], [r[key] for r in history], label=f'dt/tff={frac}')
                ax.set(xlabel='t / initial tff', ylabel=title)
                ax.grid(alpha=.3)
            axes[0, 0].legend(fontsize=8)
            fig.suptitle('Seed 42, N=300, fixed epsilon/R=0.05; uncorrected trajectories')
            fig.savefig(output/'cloud_conservation.png', dpi=150)
            plt.close(fig)
            fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
            for kind, eps in sorted({(r['kind'], r['softening']) for r in orbit_rows}):
                rows = [r for r in orbit_rows if r['kind'] == kind and r['softening'] == eps]
                axes[0].loglog([r['steps_per_period'] for r in rows],
                              [r['max_relative_position_error'] for r in rows], 'o-', label=f'{kind}, eps={eps}')
            axes[0].set(xlabel='Steps per period', ylabel='Max relative-position error / a')
            axes[0].legend(fontsize=7)
            axes[0].grid(alpha=.3)
            axes[1].plot([r['epsilon_over_separation'] for r in softening_rows],
                         [r['force_ratio_to_newtonian'] for r in softening_rows], 'o-', label='Force / Newtonian force')
            axes[1].plot([r['epsilon_over_separation'] for r in softening_rows],
                         [r['period_ratio_to_kepler'] for r in softening_rows], 'o-', label='Circular period / Kepler period')
            axes[1].set(xlabel='epsilon / separation', ylabel='Analytical model ratio')
            axes[1].legend(fontsize=8)
            axes[1].grid(alpha=.3)
            fig.savefig(output/'orbits_and_softening.png', dpi=150)
            plt.close(fig)
        finally:
            if old is None:
                os.environ.pop('MPLCONFIGDIR', None)
            else:
                os.environ['MPLCONFIGDIR'] = old


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    config_path = ROOT/'experiments/phase1_config.json'
    config = json.loads(config_path.read_text())
    n = load_baseline()  # Hash-check legacy source before allocating any output.
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter()
    (output/'config.json').write_text(json.dumps(config, indent=2)+'\n')
    np.savez(output/'cloud_initial.npz', pos=n['pos'], vel=n['vel'], mass=n['mass'])
    orbit_rows = run_orbits(config, output)
    write_csv(output/'orbit_convergence.csv', orbit_rows)
    print('Analytical orbit experiments complete.', flush=True)
    cloud_rows, histories, finals = [], [], []
    for frac in config['cloud']['dt_over_initial_tff']:
        dt = frac * n['T']
        steps = round(config['cloud']['end_time_over_initial_tff']/frac)
        history, metrics, p, v = trajectory(n['pos'], n['vel'], n['mass'], n['SOFTENING'],
                                          n['G_si'], dt, steps, n['R_cloud'])
        histories.append((frac, history))
        finals.append((p, v))
        cloud_rows.append({'dt_over_tff': frac, 'steps': steps, 'dt_s': dt, **metrics})
        write_csv(output/f'cloud_dt{frac:g}.csv', history)
        np.savez(output/f'cloud_final_dt{frac:g}.npz', pos=p, vel=v, mass=n['mass'])
        print(f'Cloud dt/tff={frac:g}: max |dE|/|E0|={metrics["max_energy_error"]:.6g}', flush=True)
    write_csv(output/'cloud_summary.csv', cloud_rows)
    # Self-convergence is empirical; chaotic close encounters need not be asymptotic.
    differences = [float(np.sqrt(np.mean((finals[i][0]-finals[i+1][0])**2))/n['R_cloud'])
                   for i in range(len(finals)-1)]
    cloud_order = math.log2(differences[0]/differences[1])
    softening_rows = []
    p = np.array([[0., 0., 0.], [1., 0., 0.]])
    for eps in config['orbits']['circular_softening']:
        measured = compute_acceleration(p, [0.4, 0.6], eps, G=1.)[0, 0]/.6
        softening_rows.append({'epsilon_over_separation': eps, 'force_ratio_to_newtonian': (1+eps**2)**-1.5,
                              'measured_force_ratio': float(measured),
                              'period_ratio_to_kepler': (1+eps**2)**.75})
    write_csv(output/'softening.csv', softening_rows)
    tol = config['tolerances_declared_before_runs']
    checks = {
        'orbit_second_order': all(tol['convergence_order_min'] <= r['observed_order'] <= tol['convergence_order_max']
                                 for r in orbit_rows if r['observed_order'] is not None),
        'orbit_finest_position_accuracy': all(r['max_relative_position_error'] <= tol['orbit_finest_relative_position']
                                             for r in orbit_rows if r['steps_per_period'] == max(config['orbits']['steps_per_period'])),
        'orbit_energy_budget': all(r['max_energy_error'] <= tol['orbit_max_relative_energy'] for r in orbit_rows),
        'cloud_energy_budget': all(r['max_energy_error'] <= tol['cloud_max_relative_energy'] for r in cloud_rows),
        'momentum_budget': all(r['max_momentum_error'] <= tol['normalized_momentum_drift'] for r in orbit_rows+cloud_rows),
        'angular_momentum_budget': all(r['max_angular_momentum_error'] <= tol['normalized_angular_momentum_drift'] for r in orbit_rows+cloud_rows),
        'com_budget': all(r['max_com_residual'] <= tol['normalized_com_residual'] for r in orbit_rows+cloud_rows)
    }
    plot_results(output, histories, orbit_rows, softening_rows, n['T'])
    tracked_inputs = ['nbody.py', 'experiments/validate_nbody.py', 'experiments/phase1_config.json',
                      'tests/baseline_support.py', 'baseline/manifest.json', 'baseline/reference.json', 'starformation.py']
    metadata = {
        'config': config, 'checks': checks, 'wall_time_s': time.perf_counter()-start,
        'python': sys.version, 'platform': platform.platform(),
        'versions': {p: importlib.metadata.version(p) for p in ('numpy', 'astropy', 'matplotlib')},
        'git_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'git_status': subprocess.check_output(['git', 'status', '--short'], cwd=ROOT, text=True),
        'input_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in tracked_inputs},
        'units': {'cloud': 'SI; normalized errors dimensionless', 'orbits': 'G=M_total=a=1 consistent dimensionless units'},
        'cloud_position_self_difference_rms_over_R': differences,
        'cloud_position_self_convergence_order': cloud_order,
        'warnings': [f'FAILED criterion: {key}' for key, passed in checks.items() if not passed],
        'scope': 'Two orbital periods; cloud 1.2 initial free-fall times. No long-term stability claim.'
    }
    (output/'results.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(json.dumps({'checks': checks, 'cloud_self_order': cloud_order,
                      'wall_time_s': metadata['wall_time_s']}, indent=2))
    # Preserve diagnostics of failed experiments, and signal the failure to automation.
    return 0 if all(checks.values()) else 1


if __name__ == '__main__':
    raise SystemExit(main())
