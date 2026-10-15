"""Phase 0 characterization only; no full simulation, rendering, or installation."""

import ast
import hashlib
import json
import math
import subprocess
import sys
import unittest

import numpy as np

from baseline_support import MANIFEST, ROOT, load_baseline, summary


class BaselineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = json.loads(MANIFEST.read_text())
        cls.n = load_baseline()

    def test_preserved_source_and_media(self):
        for name, entry in self.manifest["preserved_files"].items():
            with self.subTest(file=name):
                path = ROOT / name
                with path.open("rb") as stream:
                    actual = hashlib.file_digest(stream, "sha256").hexdigest()
                self.assertEqual(actual, entry["sha256"])
                if "bytes" in entry:
                    self.assertEqual(path.stat().st_size, entry["bytes"])

    def test_source_syntax_without_execution(self):
        compile(ast.parse((ROOT / "starformation.py").read_text()), "starformation.py", "exec")

    def test_recorded_parameters(self):
        n, p = self.n, self.manifest["parameters"]
        for key, expected in {"N_PARTICLES": p["particles"], "N_FRAMES": p["frames"],
                              "c_sound": p["sound_speed_m_s"],
                              "DT_FRAC": p["dt_over_initial_free_fall_time"]}.items():
            self.assertEqual(n[key], expected)
        self.assertEqual(n["SOFTENING"], p["softening_over_initial_radius"] * n["R_cloud"])

    def test_independent_si_reference(self):
        # Explicit CODATA 2018 G; IAU nominal GM_sun, exact AU and Julian year.
        # pc = AU / (one arcsecond in radians). No Astropy-derived references.
        g = 6.67430e-11  # m^3 kg^-1 s^-2
        mass = 1.3271244e20 / g  # kg
        radius = 0.1 * 149597870700.0 * 648000.0 / math.pi  # m
        year = 31557600.0  # s
        density = mass / (4 * math.pi * radius**3 / 3)
        free_fall = math.pi / math.sqrt(8) * math.sqrt(radius**3 / (g * mass))
        expected = [g, mass, radius, free_fall / year, 0.012 * free_fall]
        actual = [self.n[k] for k in ("G_si", "M_cloud", "R_cloud", "t_ff_yr", "dt")]
        np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=0)
        self.assertAlmostEqual(density / 1.615718379859146e-17, 1.0, places=12)
        # Units: G*rho has dimension s^-2; c_s/sqrt(G*rho) has dimension m.
        jeans = math.sqrt(math.pi * 200.0**2 / (g * density))
        jeans_mass = math.pi / 6 * density * jeans**3  # radius lambda_J/2 convention
        np.testing.assert_allclose([jeans / 149597870700.0, jeans_mass / mass],
                                   [72159.55153478573, 5.351999815649333], rtol=1e-10, atol=0)

    def test_initial_state_and_total_mass(self):
        n = self.n
        for key, shape in (("pos", (300, 3)), ("vel", (300, 3)), ("mass", (300,))):
            self.assertEqual(n[key].shape, shape)
            self.assertTrue(np.isfinite(n[key]).all())
        self.assertTrue((np.linalg.norm(n["pos"], axis=1) <= n["R_cloud"]).all())
        self.assertTrue((n["mass"] > 0).all())
        np.testing.assert_allclose(n["mass"].sum(), n["M_cloud"], rtol=1e-12, atol=0)
        np.testing.assert_array_equal(n["mass"], np.full(300, n["M_cloud"] / 300))

    def test_seeded_repeatability_and_rng_isolation(self):
        before = np.random.get_state()
        repeated = load_baseline()
        after = np.random.get_state()
        self.assertEqual(before[0], after[0])
        np.testing.assert_array_equal(before[1], after[1])
        self.assertEqual(before[2:], after[2:])
        for key in ("pos", "vel", "mass"):
            np.testing.assert_array_equal(self.n[key], repeated[key])

    def test_initial_and_one_step_regression(self):
        reference = json.loads((ROOT / "baseline" / "reference.json").read_text())
        for key, actual in summary(self.n).items():
            with self.subTest(quantity=key):
                np.testing.assert_allclose(actual, reference[key], rtol=1e-12, atol=0)

    def test_one_step_smoke_without_mutating_initial_state(self):
        n = self.n
        p, v = n["pos"].copy(), n["vel"].copy()
        a = n["compute_acceleration"](p, n["mass"], n["SOFTENING"])
        for result in n["leapfrog_step"](p, v, a, n["dt"]):
            self.assertEqual(result.shape, (300, 3))
            self.assertTrue(np.isfinite(result).all())
        np.testing.assert_array_equal(p, n["pos"])
        np.testing.assert_array_equal(v, n["vel"])

    def test_frame_time_convention(self):
        n = self.n
        # The last stored frame is before the last of 1080 integration steps.
        elapsed = (n["N_FRAMES"] - 1) * n["dt"] / 31557600.0
        np.testing.assert_allclose(elapsed, 6780685.533686639, rtol=1e-10, atol=0)
        np.testing.assert_allclose(elapsed / n["t_ff_yr"], 12.948, rtol=1e-12, atol=0)
        self.assertEqual(n["N_FRAMES"] / 18, 60.0)

    def test_adapter_does_not_initialize_staging_or_plotting(self):
        code = (
            "import sys; from baseline_support import load_baseline; n=load_baseline(); "
            "assert 'matplotlib' not in sys.modules; "
            "assert not ({'frames_pos', 'frames_vel', 'frames_jeans', 'ani', 'fig', 'OUT'} & n.keys())"
        )
        result = subprocess.run([sys.executable, "-B", "-c", code], cwd=ROOT / "tests",
                                capture_output=True, text=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
