"""Q4 phase normalization, direction conventions, and independent references."""

import shutil
import unittest

import numpy as np

from HW5.hw5_question4 import run_sweep
from models.mie import mie_phase_function, mie_sweep

try:
    from HW5.reference_mie import reference_phase_function
except ImportError:
    reference_phase_function = None


@unittest.skipUnless(shutil.which("gfortran"), "gfortran is required")
class Question4Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = run_sweep()

    def test_requested_angles_and_positive_finite_output(self):
        for values in self.results.values():
            np.testing.assert_array_equal(values["angle_deg"], np.arange(181))
            self.assertTrue(np.all(np.isfinite(values["phase_4pi"])))
            self.assertTrue(np.all(values["phase_4pi"] > 0))
            self.assertGreater(values["phase_4pi"][0], values["phase_4pi"][180])

    def test_normalization_and_first_moment(self):
        for values in self.results.values():
            self.assertAlmostEqual(values["normalization"], 1., delta=5e-10)
            self.assertAlmostEqual(values["angular_g"], values["g"], delta=5e-10)
            np.testing.assert_allclose(values["phase_per_sr"] * (4*np.pi),
                                       values["phase_4pi"], rtol=1e-14)

    def test_metadata_matches_efficiency_interface(self):
        for (n, k), values in self.results.items():
            expected = mie_sweep([values["x"]], n, k)
            for name in ("qext", "qscat", "g"):
                self.assertAlmostEqual(values[name], expected[name][0], delta=1e-13)

    def test_mirrored_angles_unsorted_and_repeated(self):
        angles = np.array([180, 0, 90, 2, 178, 2])
        for (n, k), values in self.results.items():
            sample = mie_phase_function(4., n, k, angles_degrees=angles)
            np.testing.assert_allclose(sample["phase_4pi"], values["phase_4pi"][angles],
                                       rtol=1e-12, atol=1e-12)

    @unittest.skipIf(reference_phase_function is None, "SciPy reference requires scipy")
    def test_independent_phase_and_angle_precision_regression(self):
        for (n, k), values in self.results.items():
            expected = reference_phase_function(values["x"], n, k, values["angle_deg"])
            np.testing.assert_allclose(values["phase_4pi"], expected, rtol=1e-9, atol=1e-10)

    def test_rayleigh_phase_limit(self):
        angles = np.arange(0, 181, 10)
        radius = 0.001 * 0.55 / (2*np.pi)
        values = mie_phase_function(radius, 1.5, 0., angles_degrees=angles)
        expected = 0.75 * (1 + np.cos(np.deg2rad(angles))**2)
        np.testing.assert_allclose(values["phase_4pi"], expected, rtol=5e-6, atol=1e-8)

    def test_invalid_inputs_and_workspace_limit(self):
        for options in ({"radius_um": 0}, {"real_index": 0},
                        {"absorption_index": -1}, {"wavelength_um": 0},
                        {"angles_degrees": []}, {"angles_degrees": [181]},
                        {"angles_degrees": [-1]}, {"angles_degrees": [np.nan]},
                        {"angles_degrees": np.linspace(0, 89, 184)}):
            arguments = {"radius_um": 4., "real_index": 1.5, "absorption_index": 0.}
            arguments.update(options)
            with self.assertRaises(ValueError):
                mie_phase_function(**arguments)


if __name__ == "__main__":
    unittest.main()
