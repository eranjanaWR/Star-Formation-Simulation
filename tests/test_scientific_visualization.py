"""Deterministic rendering checks use synthetic states, not claims of a valid orbit."""

import copy
import hashlib
import json
import os
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scientific_diagnostics import accuracy_history, cloud_geometry, json_ready, save_diagnostics, snapshot
from scientific_visualization import AU_M, ProjectionView, load_dataset, research_figure, validate_dataset


def fixture(coincident=False):
    base = np.array([[-1., 0., 0.], [1., 0., 0.], [0., -2., 1.], [0., 2., -1.]])*AU_M
    if coincident:
        base *= 0
    pos = np.array([base, .5*base+[2*AU_M, 0, 0]])
    mass, time = np.ones(4), np.array([0., 1.])
    half = cloud_geometry(base, mass)['r50_m']
    snapshots, density = [], []
    for p in pos:
        s, f = snapshot(p, np.zeros_like(p), mass, G=1., softening=.1*AU_M, sound_speed=1., k=2,
                        reference_half_radius=half, aperture_radius=AU_M)
        snapshots.append(s)
        density.append(f['local_density_kg_m3'])
    accuracy = accuracy_history(snapshots, time, mass_total=4., reference_radius=AU_M, G=1.)
    data = {'pos': pos, 'vel': pos*0, 'mass': mass, 'time_s': time, 'density_kg_m3': np.array(density)}
    report = json_ready({'units': 'SI', 'snapshots': snapshots, 'accuracy': accuracy,
                         'configuration': {'dt_s': 1., 'softening_m': .1*AU_M}})
    return data, report


class ScientificVisualizationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cache = tempfile.TemporaryDirectory(prefix='phase2-test-mpl-')
        cls.old_cache = os.environ.get('MPLCONFIGDIR')
        os.environ['MPLCONFIGDIR'] = cls.cache.name
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        cls.plt = plt

    @classmethod
    def tearDownClass(cls):
        cls.plt.close('all')
        if cls.old_cache is None:
            os.environ.pop('MPLCONFIGDIR', None)
        else:
            os.environ['MPLCONFIGDIR'] = cls.old_cache
        cls.cache.cleanup()

    def tearDown(self):
        self.plt.close('all')

    def test_projection_values_units_and_following_view(self):
        data, report = fixture()
        for mode in ('com', 'fixed'):
            display = ProjectionView(data, report, view=mode)
            display.update(1)
            expected = data['pos'][1].copy()
            if mode == 'com':
                expected -= report['snapshots'][1]['com']
            for artist, j in zip(display.scatters, (1, 2)):
                np.testing.assert_allclose(artist.get_offsets(), expected[:, [0, j]]/AU_M, rtol=1e-12, atol=0)
            self.assertEqual(display.axes[0].get_aspect(), 1.)
            self.assertIn('[AU]', display.axes[0].get_xlabel())

    def test_fixed_colour_scale_labels_and_persistent_failure(self):
        data, report = fixture()
        display = ProjectionView(data, report)
        limits = (display.norm.vmin, display.norm.vmax)
        self.assertEqual(limits, (data['density_kg_m3'].min(), data['density_kg_m3'].max()))
        for index in (0, 1, 0):
            display.update(index)
            self.assertEqual((display.norm.vmin, display.norm.vmax), limits)
            self.assertIn('FAIL', display.accuracy_text.get_text())
            self.assertIn('run peak=', display.accuracy_text.get_text())
        labels = ' '.join(t.get_text() for t in display.figure.texts)
        self.assertNotIn('Protostar', labels)
        self.assertNotIn('temperature', labels.lower())
        self.assertIn('no pressure forces', labels)

    def test_clipping_notice_and_automatic_full_extent(self):
        data, report = fixture()
        auto = ProjectionView(data, report, view='fixed')
        self.assertGreater(auto.half_width, abs(auto.projected).max())
        display = ProjectionView(data, report, half_width_au=.1)
        self.assertIn('XY 4/4', display.clipping.get_text())
        self.assertIn('XZ 4/4', display.clipping.get_text())

    def test_undefined_density_is_grey_without_fabricated_colourbar(self):
        data, report = fixture(coincident=True)
        display = ProjectionView(data, report)
        self.assertIsNone(display.norm)
        self.assertEqual(len(display.figure.axes), 2)
        self.assertIn('undefined', display.measurements.get_text())
        self.assertIn('density (grey): 4', display.clipping.get_text())
        display.figure.canvas.draw()

    def test_png_render_and_research_diagnostics(self):
        data, report = fixture()
        view = ProjectionView(data, report)
        research = research_figure(data, report)
        self.assertEqual(len(research.axes), 9)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'test.png'
            view.figure.savefig(path, dpi=60)
            self.assertGreater(path.stat().st_size, 1000)
            self.assertEqual(path.read_bytes()[:8], b'\x89PNG\r\n\x1a\n')

    def test_annotation_does_not_overlap_axis_labels(self):
        data, report = fixture()
        display = ProjectionView(data, report)
        display.figure.canvas.draw()
        renderer = display.figure.canvas.get_renderer()
        annotation = display.measurements.get_window_extent(renderer)
        for axis in display.axes:
            self.assertFalse(annotation.overlaps(axis.xaxis.label.get_window_extent(renderer)))

    def test_mixed_undefined_density_keeps_particle_visible_in_grey(self):
        data, report = fixture()
        data['density_kg_m3'][0, 0] = np.nan
        display = ProjectionView(data, report)
        display.figure.canvas.draw()
        np.testing.assert_allclose(display.scatters[0].get_facecolors()[0], [.6, .6, .6, 1.])
        self.assertEqual(len(display.scatters[0].get_offsets()), 4)

    def test_serialized_dataset_and_hash_validation(self):
        data, report = fixture()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            np.savez(path/'trajectory.npz', **data)
            report['trajectory_sha256'] = hashlib.sha256((path/'trajectory.npz').read_bytes()).hexdigest()
            save_diagnostics(path/'diagnostics.json', report)
            restored, restored_report = load_dataset(path)
            for key in data:
                np.testing.assert_array_equal(data[key], restored[key])
            self.assertEqual(report, restored_report)
            # A mismatched reference is rejected before plotting.
            report['trajectory_sha256'] = 'incorrect'
            with (path/'diagnostics.json').open('w') as stream:
                json.dump(report, stream)
            with self.assertRaises(ValueError):
                load_dataset(path)

    def test_false_accuracy_label_and_invalid_shapes_are_rejected(self):
        data, report = fixture()
        bad_report = copy.deepcopy(report)
        bad_report['accuracy']['energy_status'] = 'PASS (sampled interval only)'
        with self.assertRaises(ValueError):
            validate_dataset(data, bad_report)
        bad = dict(data, pos=data['pos'][:, :, :2])
        with self.assertRaises(ValueError):
            validate_dataset(bad, report)
        with self.assertRaises(ValueError):
            ProjectionView(data, report, half_width_au=-1)
        bad_report = copy.deepcopy(report)
        bad_report['snapshots'][1]['com'][0] += AU_M
        with self.assertRaises(ValueError):
            validate_dataset(data, bad_report)
        bad_report = copy.deepcopy(report)
        bad_report['accuracy']['sample_end_s'] = 100.
        with self.assertRaises(ValueError):
            validate_dataset(data, bad_report)


if __name__ == '__main__':
    unittest.main()
