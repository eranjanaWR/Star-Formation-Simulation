"""Phase 3 boundaries and numerical regressions; no repository outputs."""

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))

from starformation.cloud import snapshot
from starformation.config import from_input, legacy_input
from starformation.engine import run
from starformation.gravity import compute_acceleration
from starformation.initial_conditions import make_initial_state
from starformation.integration import leapfrog_step
from starformation.serialization import load_result, save_result
from scientific_diagnostics import snapshot as legacy_snapshot
from nbody import compute_acceleration as legacy_acceleration


def short_config(n=8, steps=2):
    values = legacy_input()
    values.update(particle_count=n, neighbours=3, steps=steps)
    return from_input(values)


class ArchitectureTests(unittest.TestCase):
    def test_package_and_wrapper_identity(self):
        self.assertIs(legacy_acceleration, compute_acceleration)
        self.assertIs(legacy_snapshot, snapshot)

    def test_seed42_state_bitwise_phase1(self):
        state = make_initial_state(from_input(legacy_input()))
        with np.load(ROOT/'validation/phase1/cloud_initial.npz', allow_pickle=False) as ref:
            for key in ('pos', 'vel', 'mass'):
                np.testing.assert_array_equal(getattr(state, key), ref[key])

    def test_initial_state_does_not_change_global_rng(self):
        before = np.random.get_state()
        make_initial_state(short_config())
        after = np.random.get_state()
        for a, b in zip(before, after):
            np.testing.assert_array_equal(a, b)

    def test_explicit_arrays_and_opt_in_recentering(self):
        config = short_config()
        state = make_initial_state(config)
        explicit = make_initial_state(config, pos=state.pos, vel=state.vel, mass=state.mass)
        np.testing.assert_array_equal(explicit.pos, state.pos)
        self.assertEqual(explicit.method, 'explicit_arrays')
        values = legacy_input()
        values.update(particle_count=8, neighbours=3, steps=2, recenter_com=True,
                      recenter_velocity=True)
        corrected = make_initial_state(from_input(values))
        self.assertLess(np.linalg.norm(np.average(corrected.pos, axis=0, weights=corrected.mass)), 2)
        self.assertLess(np.linalg.norm(np.average(corrected.vel, axis=0, weights=corrected.mass)), 1e-13)
        self.assertFalse(np.array_equal(corrected.pos, state.pos))

    def test_config_rejects_invalid_fields(self):
        for key, bad in [('particle_count', True), ('seed', -1), ('steps', 0),
                         ('dt_over_initial_tff', 0), ('softening_over_radius', -1),
                         ('cloud_mass_msun', float('nan')), ('neighbours', 300)]:
            values = legacy_input()
            values[key] = bad
            with self.subTest(key=key), self.assertRaises(ValueError):
                from_input(values)
        values = legacy_input()
        values['schema_version'] = 2
        with self.assertRaises(ValueError):
            from_input(values)

    def test_trajectory_and_diagnostic_equivalence(self):
        config = from_input({**legacy_input(), 'steps': 1})
        result = run(config)
        with np.load(ROOT/'validation/phase1/cloud_initial.npz', allow_pickle=False) as ref:
            p, v, m = ref['pos'], ref['vel'], ref['mass']
        a = legacy_acceleration(p, m, config.softening_m, G=config.gravitational_constant_si)
        p1, v1, _ = leapfrog_step(p, v, a, config.dt_s, m, config.softening_m,
                                 G=config.gravitational_constant_si)
        np.testing.assert_array_equal(result.trajectory['pos'][1], p1)
        np.testing.assert_array_equal(result.trajectory['vel'][1], v1)
        measured, _ = legacy_snapshot(p, v, m, G=config.gravitational_constant_si,
                                      softening=config.softening_m,
                                      sound_speed=config.sound_speed_m_s, k=config.neighbours,
                                      reference_half_radius=result.report['configuration']['initial_half_radius_m'],
                                      aperture_radius=config.radius_m/2)
        self.assertEqual(result.report['snapshots'][0]['energy'], measured['energy'])

    def test_archived_failed_trajectory_and_status(self):
        config = from_input({**legacy_input(), 'steps': 100})
        result = run(config)
        with np.load(ROOT/'validation/phase2/dt0.012/trajectory.npz', allow_pickle=False) as archived:
            for key in ('pos', 'vel', 'mass', 'time_s'):
                np.testing.assert_array_equal(result.trajectory[key], archived[key])
        recorded = json.loads((ROOT/'validation/phase2/dt0.012/diagnostics.json').read_text())
        self.assertEqual(result.report['accuracy']['energy_status'], 'FAIL')
        for key in ('energy', 'kinetic', 'potential', 'r50_m', 'inflow_mass_fraction'):
            self.assertEqual([s[key] for s in result.report['snapshots']],
                             [s[key] for s in recorded['snapshots']])
        self.assertAlmostEqual(result.report['accuracy']['peak_energy_error'],
                               recorded['accuracy']['peak_energy_error'], places=14)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'failed'
            save_result(result, path)
            saved = load_result(path)
            self.assertEqual(saved['validation']['status'], 'FAIL')
            self.assertIn('FAILED', saved['validation']['warnings'][0])

    def test_determinism_and_round_trip(self):
        config = short_config()
        first, second = run(config), run(config)
        for key in first.trajectory:
            np.testing.assert_array_equal(first.trajectory[key], second.trajectory[key])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'experiment'
            save_result(first, path)
            saved = load_result(path)
            for key in first.trajectory:
                np.testing.assert_array_equal(saved['trajectory'][key], first.trajectory[key])
            self.assertEqual(saved['configuration']['units'], 'SI')
            self.assertEqual(saved['validation']['status'], first.report['accuracy']['energy_status'])
            with self.assertRaises(FileExistsError):
                save_result(first, path)

    def test_unknown_schema_and_incomplete_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'experiment'
            save_result(run(short_config()), path)
            marker = json.loads((path/'COMPLETE.json').read_text())
            marker['schema_version'] = 2
            (path/'COMPLETE.json').write_text(json.dumps(marker))
            with self.assertRaises(ValueError):
                load_result(path)
            (path/'COMPLETE.json').unlink()
            with self.assertRaisesRegex(ValueError, 'Incomplete dataset'):
                load_result(path)

    def test_import_and_cli_isolation(self):
        env = dict(__import__('os').environ, PYTHONPATH=str(ROOT/'src'))
        code = ('import sys; import starformation, starformation.engine, '
                'starformation.cli_simulate; '
                'assert "matplotlib" not in sys.modules; '
                'assert "starformation.cli_render" not in sys.modules; '
                'assert "starformation.py" not in starformation.__file__')
        check = subprocess.run([sys.executable, '-B', '-c', code], cwd=ROOT,
                               env=env, capture_output=True, text=True)
        self.assertEqual(check.returncode, 0, check.stderr)
        help_run = subprocess.run([sys.executable, '-B', '-m',
                                   'starformation.cli_simulate', '--help'], cwd=ROOT,
                                  env=env, capture_output=True, text=True)
        self.assertEqual(help_run.returncode, 0, help_run.stderr)
        self.assertIn('--output', help_run.stdout)
        help_render = subprocess.run([sys.executable, '-B', '-m',
                                      'starformation.cli_render', '--help'], cwd=ROOT,
                                     env=env, capture_output=True, text=True)
        self.assertEqual(help_render.returncode, 0, help_render.stderr)
        self.assertIn('--output', help_render.stdout)
        missing_output = subprocess.run([sys.executable, '-B', '-m',
                                         'starformation.cli_simulate', '--config',
                                         str(ROOT/'configs/phase3_short.json')], cwd=ROOT,
                                        env=env, capture_output=True, text=True)
        self.assertEqual(missing_output.returncode, 2)
        self.assertIn('required: --output', missing_output.stderr)


if __name__ == '__main__':
    unittest.main()
