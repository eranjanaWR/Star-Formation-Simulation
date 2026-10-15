"""Isolated Phase 1 checks. Tolerances precede measured experiments."""

import unittest

import numpy as np

from baseline_support import load_baseline
from nbody import compute_acceleration as acceleration, diagnostics, leapfrog_step, potential_energy


class NBodyTests(unittest.TestCase):
    def setUp(self):
        self.p = np.array([[-0.7, 0.2, 0.1], [0.4, -0.6, 0.3], [0.8, 0.9, -0.5]])
        self.m = np.array([0.2, 0.3, 0.5])
        self.v = np.array([[0.1, 0.3, -0.2], [-0.4, 0.2, 0.1], [0.2, -0.1, 0.3]])

    def test_two_body_direction_mass_weighting_and_symmetry(self):
        p, m = np.array([[0., 0., 0.], [3., 4., 0.]]), np.array([2., 7.])
        for eps in (0., 0.2, 2.):
            with self.subTest(eps=eps):
                actual = acceleration(p, m, eps, G=1.7)
                expected = 1.7 * np.array([7 * p[1], -2 * p[1]]) / (25 + eps**2)**1.5
                np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=0)
                np.testing.assert_allclose(np.sum(m[:, None] * actual, axis=0), 0, atol=1e-14)

    def test_single_particle_has_no_self_force_or_energy(self):
        for eps in (0., 0.5):
            p = np.array([[1., 2., 3.]])
            np.testing.assert_array_equal(acceleration(p, [2.], eps, G=1.), [[0., 0., 0.]])
            self.assertEqual(potential_energy(p, [2.], eps, G=1.), 0.)

    def test_coincident_softened_pair_and_unsoftened_singularity(self):
        p = np.zeros((2, 3))
        np.testing.assert_array_equal(acceleration(p, [2., 3.], 0.5, G=1.), p)
        self.assertEqual(potential_energy(p, [2., 3.], 0.5, G=1.), -12.)
        with self.assertRaises(FloatingPointError):
            acceleration(p, [2., 3.], 0., G=1.)
        with self.assertRaises(FloatingPointError):
            potential_energy(p, [2., 3.], 0., G=1.)

    def test_potential_counts_each_pair_once(self):
        expected = sum(-self.m[i] * self.m[j] /
                       np.sqrt(np.dot(self.p[i]-self.p[j], self.p[i]-self.p[j]) + .2**2)
                       for i in range(3) for j in range(i+1, 3))
        np.testing.assert_allclose(potential_energy(self.p, self.m, .2, G=1.), expected,
                                   rtol=1e-12, atol=0)

    def test_force_is_negative_gradient_of_softened_potential(self):
        for eps in (0., .05, .5):
            forces = self.m[:, None] * acceleration(self.p, self.m, eps, G=1.)
            finite_difference = np.zeros_like(forces)
            for i in range(3):
                for j in range(3):
                    plus, minus = self.p.copy(), self.p.copy()
                    plus[i, j] += 1e-5
                    minus[i, j] -= 1e-5
                    finite_difference[i, j] = -(potential_energy(plus, self.m, eps, G=1.) -
                                                potential_energy(minus, self.m, eps, G=1.)) / 2e-5
            np.testing.assert_allclose(forces, finite_difference, rtol=1e-6, atol=0)

    def test_total_force_and_torque_cancel(self):
        force = self.m[:, None] * acceleration(self.p, self.m, .1, G=1.)
        np.testing.assert_allclose(force.sum(axis=0), 0, atol=1e-14)
        np.testing.assert_allclose(np.cross(self.p, force).sum(axis=0), 0, atol=1e-14)

    def test_translation_rotation_and_permutation_covariance(self):
        original = acceleration(self.p, self.m, .1, G=1.)
        np.testing.assert_allclose(acceleration(self.p + [4., -3., 2.], self.m, .1, G=1.),
                                   original, rtol=1e-12, atol=1e-15)
        q = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        np.testing.assert_allclose(acceleration(self.p @ q.T, self.m, .1, G=1.),
                                   original @ q.T, rtol=1e-12, atol=1e-15)
        order = [2, 0, 1]
        np.testing.assert_allclose(acceleration(self.p[order], self.m[order], .1, G=1.),
                                   original[order], rtol=1e-12, atol=1e-15)

    def test_diagnostics_against_hand_computed_case(self):
        d = diagnostics([[-1., 0., 0.], [2., 0., 0.]], [[0., 2., 0.], [0., -1., 0.]],
                        [2., 1.], 4., G=5.)
        self.assertEqual(d['kinetic'], 4.5)
        self.assertEqual(d['potential'], -2.)
        self.assertEqual(d['energy'], 2.5)
        np.testing.assert_array_equal(d['momentum'], [0., 3., 0.])
        np.testing.assert_array_equal(d['angular_momentum'], [0., 0., -6.])
        np.testing.assert_array_equal(d['com'], [0., 0., 0.])
        np.testing.assert_array_equal(d['com_velocity'], [0., 1., 0.])

    def test_unit_scaling_including_integrator_and_diagnostics(self):
        # L=3e15 m, M=2e30 kg; T chosen so dimensionless G=1.
        L, M, G = 3e15, 2e30, 6.67430e-11
        T = np.sqrt(L**3 / (G * M))
        a = acceleration(self.p, self.m, .1, G=1.)
        a_si = acceleration(self.p*L, self.m*M, .1*L, G=G)
        np.testing.assert_allclose(a_si / (L/T**2), a, rtol=1e-12, atol=0)
        x = leapfrog_step(self.p, self.v, a, .01, self.m, .1, G=1.)
        y = leapfrog_step(self.p*L, self.v*L/T, a_si, .01*T, self.m*M, .1*L, G=G)
        for left, right, scale in zip(x, y, (L, L/T, L/T**2)):
            np.testing.assert_allclose(right/scale, left, rtol=1e-12, atol=0)
        d = diagnostics(self.p, self.v, self.m, .1, G=1.)
        s = diagnostics(self.p*L, self.v*L/T, self.m*M, .1*L, G=G)
        for key, scale in {'energy': M*L**2/T**2, 'kinetic': M*L**2/T**2,
                           'potential': M*L**2/T**2, 'momentum': M*L/T,
                           'angular_momentum': M*L**2/T, 'com': L, 'com_velocity': L/T}.items():
            np.testing.assert_allclose(s[key]/scale, d[key], rtol=1e-12, atol=1e-15)

    def test_no_gravity_exact_linear_motion(self):
        p, v = self.p.copy(), self.v.copy()
        a = acceleration(p, self.m, 0., G=0.)
        for _ in range(100):
            p, v, a = leapfrog_step(p, v, a, .01, self.m, 0., G=0.)
        np.testing.assert_allclose(p, self.p + self.v, rtol=1e-12, atol=1e-14)
        np.testing.assert_array_equal(v, self.v)
        self.assertEqual(potential_energy(p, self.m, 0., G=0.), 0.)

    def test_kdk_reference_and_time_reversal(self):
        a = acceleration(self.p, self.m, .1, G=1.)
        half = self.v + .005 * a
        expected_p = self.p + .01 * half
        expected_a = acceleration(expected_p, self.m, .1, G=1.)
        result = leapfrog_step(self.p, self.v, a, .01, self.m, .1, G=1.)
        for actual, expected in zip(result, (expected_p, half + .005*expected_a, expected_a)):
            np.testing.assert_array_equal(actual, expected)
        p, v, a = result
        p, v, _ = leapfrog_step(p, v, a, -.01, self.m, .1, G=1.)
        np.testing.assert_allclose(p, self.p, rtol=1e-12, atol=0)
        np.testing.assert_allclose(v, self.v, rtol=1e-12, atol=0)

    def test_seed42_legacy_equivalence_for_complete_arrays(self):
        n = load_baseline()
        p, v = n['pos'].copy(), n['vel'].copy()
        old_p, old_v = p.copy(), v.copy()
        a = acceleration(p, n['mass'], n['SOFTENING'], G=n['G_si'])
        old_a = n['compute_acceleration'](old_p, n['mass'], n['SOFTENING'])
        np.testing.assert_array_equal(a, old_a)
        for _ in range(3):
            p, v, a = leapfrog_step(p, v, a, n['dt'], n['mass'], n['SOFTENING'], G=n['G_si'])
            old_p, old_v, old_a = n['leapfrog_step'](old_p, old_v, old_a, n['dt'])
            for actual, expected in zip((p, v, a), (old_p, old_v, old_a)):
                np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=0)
        np.testing.assert_array_equal(n['pos'], load_baseline()['pos'])

    def test_invalid_inputs_are_reported_not_corrected(self):
        for mass, eps, G in (([-1., 1., 1.], .1, 1.), (self.m, -.1, 1.),
                             (self.m, .1, -1.), (self.m, np.nan, 1.)):
            with self.assertRaises(ValueError):
                acceleration(self.p, mass, eps, G=G)
        with self.assertRaises(ValueError):
            acceleration([[np.nan, 0., 0.]], [1.], .1, G=1.)
        with self.assertRaises(ValueError):
            acceleration(self.p, [1.], .1, G=1.)
        with self.assertRaises(ValueError):
            leapfrog_step(self.p, self.v, np.zeros((2, 3)), .1, self.m, .1, G=1.)
        from astropy import units as u
        with self.assertRaises(TypeError):
            acceleration(self.p*u.m, self.m, .1, G=1.)


if __name__ == '__main__':
    unittest.main()
