"""Phase 2 tolerances: analytical 1e-12 relative; uniform interior mean within 10%.

Synthetic random clouds are deterministic. No trajectory correction is permitted.
"""

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scientific_diagnostics import (accuracy_history, centre_of_mass, cloud_geometry, contraction,
    enclosed_mass, enclosed_mass_radii, jeans_scales, local_particle_density, radial_flow,
    save_diagnostics, snapshot, virial_measurements)


def uniform_density_experiment():
    rng = np.random.default_rng(2026)
    pos = rng.uniform(-.5, .5, (2048, 3))
    mass = np.full(len(pos), 1/len(pos))
    density, support = local_particle_density(pos, mass, k=32)
    interior = np.all(abs(pos)+support[:, None] < .5, axis=1)
    return {'interior_count': int(interior.sum()), 'interior_mean': float(density[interior].mean()),
            'interior_std': float(density[interior].std()), 'true_density': 1.,
            'all_particle_mean': float(density.mean()),
            'interior_mean_relative_error': float(abs(density[interior].mean()-1))}


class ScientificDiagnosticsTests(unittest.TestCase):
    def setUp(self):
        self.p = np.array([[-1., 0., 0.], [3., 0., 0.]])
        self.m = np.array([3., 1.])

    def test_mass_weighted_com_and_radii(self):
        np.testing.assert_array_equal(centre_of_mass(self.p, self.m), [0., 0., 0.])
        np.testing.assert_array_equal(enclosed_mass_radii(self.p, self.m, [.1, .5, .75, .9, 1]), [1, 1, 1, 3, 3])
        g = cloud_geometry(self.p, self.m)
        np.testing.assert_allclose(g['rms_radius_m'], np.sqrt(3), rtol=1e-12)
        self.assertEqual(g['mass_within_r50_kg'], 3.)
        self.assertEqual(enclosed_mass(self.p, self.m, 2), 3.)
        self.assertAlmostEqual(g['mean_density_half_kg_m3'], 3/(4*np.pi/3), places=12)

    def test_mass_radius_ties_include_actual_mass_not_nominal_fraction(self):
        p = np.array([[-1., 0, 0], [1., 0, 0]])
        g = cloud_geometry(p, [1., 1.])
        self.assertEqual(g['r50_m'], 1.)
        self.assertEqual(g['mass_within_r50_kg'], 2.)
        self.assertAlmostEqual(g['mean_density_half_kg_m3'], 2/(4*np.pi/3))

    def test_translation_and_uniform_velocity_boost(self):
        vel = np.array([[0., 1., 0.], [0., -3., 0.]])
        g = cloud_geometry(self.p, self.m)
        self.assertEqual(g, cloud_geometry(self.p+[100., 2., -5.], self.m))
        a = virial_measurements(self.p, vel, self.m, .2, G=1.)
        b = virial_measurements(self.p+10, vel+[2, 3, 4], self.m, .2, G=1.)
        np.testing.assert_allclose(a['kinetic_internal'], b['kinetic_internal'], rtol=1e-12)
        self.assertAlmostEqual(b['kinetic'], b['kinetic_internal']+b['kinetic_bulk'])
        np.testing.assert_allclose(radial_flow(self.p, vel, self.m)['radial_velocity_m_s'],
                                   radial_flow(self.p+10, vel+[2, 3, 4], self.m)['radial_velocity_m_s'])

    def test_local_density_hand_case_and_mass_scaling(self):
        p = np.array([[0., 0, 0], [1., 0, 0], [3., 0, 0], [5., 0, 0]])
        mass = np.array([100., 2., 4., 8.])
        density, support = local_particle_density(p, mass, k=2)
        self.assertEqual(support[0], 3.)
        self.assertAlmostEqual(density[0], .5*(2+4)/(4*np.pi*3**3/3))
        scaled, _ = local_particle_density(2*p+7, 3*mass, k=2)
        np.testing.assert_allclose(scaled, density*3/8, rtol=1e-12)

    def test_neighbour_ties_are_permutation_invariant_for_unequal_masses(self):
        p = np.vstack([np.zeros(3), np.eye(3), -np.eye(3)])
        mass = np.arange(1, 8, dtype=float)
        a, _ = local_particle_density(p, mass, k=3)
        order = [5, 1, 3, 0, 6, 2, 4]
        b, _ = local_particle_density(p[order], mass[order], k=3)
        np.testing.assert_allclose(b, a[order], rtol=1e-12)
        self.assertAlmostEqual(a[0], (2/3)*(mass[1:].sum()/2)/(4*np.pi/3))

    def test_uniform_synthetic_density_estimator(self):
        result = uniform_density_experiment()
        self.assertGreater(result['interior_count'], 300)
        self.assertLessEqual(result['interior_mean_relative_error'], .10)

    def test_uniform_volume_sampling_rms_and_half_mass_radius(self):
        rng = np.random.default_rng(12)
        direction = rng.normal(size=(1000, 3))
        direction /= np.linalg.norm(direction, axis=1)[:, None]
        p = direction * (((np.arange(1000)+.5)/1000)**(1/3))[:, None]
        # Antipodal pairs have exactly zero COM up to summation roundoff.
        p = np.vstack([p, -p])
        g = cloud_geometry(p, np.ones(len(p)))
        self.assertLess(abs(g['r50_m']/.5**(1/3)-1), .001)
        self.assertLess(abs(g['rms_radius_m']/np.sqrt(3/5)-1), .001)

    def test_jeans_independent_formula_and_scalings(self):
        d = jeans_scales(2., 3., G=4.)
        length = 3*np.sqrt(np.pi/8)
        np.testing.assert_allclose(d['jeans_length_m'], length, rtol=1e-12)
        np.testing.assert_allclose(d['jeans_mass_kg'], 2*4*np.pi/3*(length/2)**3, rtol=1e-12)
        np.testing.assert_allclose(d['free_fall_s'], np.sqrt(3*np.pi/256), rtol=1e-12)
        denser = jeans_scales(8., 3., G=4.)
        self.assertAlmostEqual(denser['jeans_length_m']/d['jeans_length_m'], .5)
        self.assertAlmostEqual(denser['jeans_mass_kg']/d['jeans_mass_kg'], .5)

    def test_softened_virial_distinguishes_force_and_potential(self):
        p = np.array([[-.5, 0., 0.], [.5, 0., 0.]])
        omega = np.sqrt(2/(1+.5**2)**1.5)
        vel = np.array([[0., -.5*omega, 0.], [0., .5*omega, 0.]])
        d = virial_measurements(p, vel, [1., 1.], .5, G=1.)
        self.assertAlmostEqual(d['virial_ratio_W'], 1.)
        self.assertAlmostEqual(d['virial_ratio_U'], .8)
        self.assertAlmostEqual(d['virial_W'], -1/(1+.5**2)**1.5)
        zero = virial_measurements(p, vel*0, [1., 1.], .5, G=1.)
        self.assertEqual(zero['virial_ratio_U'], 0.)

    def test_expansion_and_contraction_metrics_and_zero_velocity(self):
        for factor in (-2., 0., 3.):
            d = radial_flow(self.p, factor*self.p, self.m)
            self.assertAlmostEqual(d['mean_radial_velocity_m_s'], factor*1.5)
            self.assertEqual(d['inflow_mass_fraction'], float(factor < 0))
        self.assertEqual(contraction(.5, 1.), .5)
        self.assertEqual(contraction(2., 1.), -1.)
        self.assertEqual(contraction(0., 1.), 1.)

    def test_fragmented_cloud_is_not_mistaken_for_low_local_density(self):
        rng = np.random.default_rng(4)
        clump = rng.normal(size=(40, 3))*.01
        p = np.vstack([clump+[-10, 0, 0], clump+[10, 0, 0]])
        m = np.ones(len(p))
        rho, _ = local_particle_density(p, m, k=8)
        g = cloud_geometry(p, m)
        self.assertGreater(np.median(rho)/g['mean_density_half_kg_m3'], 1e6)
        self.assertGreater(g['r50_m'], 9.)

    def test_coincident_positions_report_undefined_estimates(self):
        p = np.zeros((4, 3))
        s, fields = snapshot(p, p, np.ones(4), G=1., softening=.1, sound_speed=1., k=2,
                             reference_half_radius=0., aperture_radius=1.)
        self.assertIsNone(s['mean_density_half_kg_m3'])
        self.assertIsNone(s['jeans_length_m'])
        self.assertIsNone(s['virial_ratio_W'])
        self.assertIsNone(s['contraction_half'])
        self.assertEqual(s['undefined_density_count'], 4)
        self.assertTrue(np.isnan(fields['local_density_kg_m3']).all())
        self.assertIsNone(s['mean_radial_velocity_m_s'])

    def test_snapshot_uses_current_enclosed_diameter_for_jeans_comparison(self):
        p = np.array([[-2., 0, 0], [-1., 0, 0], [1., 0, 0], [2., 0, 0]])
        s, _ = snapshot(p, p*0, np.ones(4), G=1., softening=.1, sound_speed=1., k=2,
                        reference_half_radius=1., aperture_radius=1.)
        self.assertEqual(s['jeans_characteristic_length_m'], 2.)
        self.assertEqual(s['mass_within_r50_kg'], 2.)
        self.assertAlmostEqual(s['diameter_half_over_jeans'], 2/np.sqrt(np.pi/(2/(4*np.pi/3))))
        self.assertEqual(s['aperture_mass_fraction'], .5)

    def test_bulk_translation_com_history(self):
        velocity = np.array([1., 2., -3.])
        samples = [{'energy': -1., 'momentum': velocity*4., 'angular_momentum': np.zeros(3),
                    'com': np.array([3., 4., 5.])+t*velocity, 'com_velocity': velocity} for t in (0., 1., 2.)]
        a = accuracy_history(samples, [0., 1., 2.], mass_total=4., reference_radius=1., G=1.)
        np.testing.assert_array_equal(a['com_residual'], np.zeros(3))
        self.assertEqual(a['energy_status'], 'PASS (sampled interval only)')

    def test_accuracy_keeps_peak_failure_even_after_energy_recovers(self):
        samples = [{'energy': e, 'momentum': np.zeros(3), 'angular_momentum': np.zeros(3),
                    'com': np.zeros(3), 'com_velocity': np.zeros(3)} for e in (-1., -1.02, -1.)]
        a = accuracy_history(samples, [0., 1., 2.], mass_total=1., reference_radius=1., G=1.)
        self.assertEqual(a['energy_status'], 'FAIL')
        self.assertEqual(a['energy_error'][-1], 0.)
        self.assertAlmostEqual(a['running_peak_energy_error'][-1], .02)

    def test_strict_serialization_and_no_overwrite(self):
        payload = {'array': np.array([1., 2., np.nan]), 'undefined': None, 'value': np.float64(.123456789012345)}
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)/'diagnostics.json'
            save_diagnostics(p, payload)
            restored = json.loads(p.read_text(), parse_constant=lambda s: self.fail(f'Invalid JSON {s}'))
            self.assertEqual(restored, {'array': [1., 2., None], 'undefined': None, 'value': float(payload['value'])})
            with self.assertRaises(FileExistsError):
                save_diagnostics(p, payload)

    def test_invalid_inputs(self):
        for p, m in (([], []), ([[0., 0]], [1]), ([[0, 0, np.nan]], [1]), (self.p, [1, -1])):
            with self.assertRaises(ValueError):
                cloud_geometry(p, m)
        for fractions in ([], [0], [1.1], [np.nan]):
            with self.assertRaises(ValueError):
                enclosed_mass_radii(self.p, self.m, fractions)
        for k in (0, 1, 2, 2.5, True):
            with self.assertRaises(ValueError):
                local_particle_density(self.p, self.m, k)
        for density in (0, -1, np.inf):
            with self.assertRaises(ValueError):
                jeans_scales(density, 1., G=1.)
        with self.assertRaises(ValueError):
            contraction(1., 0.)
        with self.assertRaises(ValueError):
            radial_flow(self.p, np.zeros((3, 3)), self.m)


if __name__ == '__main__':
    unittest.main()
