"""Keep analytic orbital convergence in the normal automated test suite."""

import json
import unittest

from baseline_support import ROOT
from experiments.validate_nbody import run_orbits


class OrbitTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        config = json.loads((ROOT / 'experiments/phase1_config.json').read_text())
        cls.tol = config['tolerances_declared_before_runs']
        # A small subset of the archived experiment; no plots or data files.
        config['orbits']['circular_softening'] = [0., .2]
        config['orbits']['steps_per_period'] = [200, 400, 800]
        cls.rows = run_orbits(config)

    def test_analytic_orbits_show_second_order_position_convergence(self):
        for row in self.rows:
            if row['observed_order'] is not None:
                with self.subTest(kind=row['kind'], eps=row['softening'], steps=row['steps_per_period']):
                    self.assertGreaterEqual(row['observed_order'], self.tol['convergence_order_min'])
                    self.assertLessEqual(row['observed_order'], self.tol['convergence_order_max'])

    def test_finest_orbit_position_error(self):
        for row in self.rows:
            if row['steps_per_period'] == 800:
                self.assertLessEqual(row['max_relative_position_error'], self.tol['orbit_finest_relative_position'])

    def test_orbit_energy_budget(self):
        for row in self.rows:
            self.assertLessEqual(row['max_energy_error'], self.tol['orbit_max_relative_energy'])

    def test_orbit_momentum_angular_momentum_and_com(self):
        for row in self.rows:
            for metric, tolerance in [('max_momentum_error', 'normalized_momentum_drift'),
                                      ('max_angular_momentum_error', 'normalized_angular_momentum_drift'),
                                      ('max_com_residual', 'normalized_com_residual')]:
                self.assertLessEqual(row[metric], self.tol[tolerance])


if __name__ == '__main__':
    unittest.main()
