"""Mass normalization, limits, and independent validation for HW5 Q3."""

import contextlib
import io
import shutil
import unittest

import numpy as np

from HW5.hw5_question3 import RADII_UM, SPECIFIC_COLUMNS, run_sweep
from models.mie import ROOT, mass_specific_cross_section, mie_sweep

try:
    from HW5.reference_mie import reference_efficiencies
except ImportError:
    reference_efficiencies = None


class MassSpecificTests(unittest.TestCase):
    def test_conversion_against_explicit_area_and_mass(self):
        radius_um, density, q = 0.2, 1.7, 2.3
        radius_m = radius_um * 1e-6
        area_m2 = np.pi * radius_m**2
        mass_g = (4 / 3) * np.pi * radius_m**3 * density * 1e6
        self.assertAlmostEqual(
            mass_specific_cross_section(q, radius_um, density), q * area_m2 / mass_g
        )

    def test_density_scaling_and_signed_residual(self):
        first = mass_specific_cross_section([1., -1e-15], [0.1, 1.], 1.)
        second = mass_specific_cross_section([1., -1e-15], [0.1, 1.], 2.)
        np.testing.assert_allclose(second, first / 2)
        self.assertLess(first[1], 0)

    def test_invalid_inputs(self):
        for q, r, rho in ((1, 0, 1), (1, -1, 1), (1, 1, 0),
                          (1, 1, np.nan), (np.nan, 1, 1)):
            with self.assertRaises(ValueError):
                mass_specific_cross_section(q, r, rho)


@unittest.skipUnless(shutil.which("gfortran"), "gfortran is required")
class Question3Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with contextlib.redirect_stdout(io.StringIO()):
            cls.results = run_sweep()

    def test_full_grid_and_specific_partition(self):
        self.assertEqual(len(self.results), 2)
        self.assertEqual(RADII_UM[0], 0.001)
        self.assertEqual(RADII_UM[-1], 10.)
        for values in self.results.values():
            self.assertTrue(all(np.all(np.isfinite(v)) for v in values.values()))
            np.testing.assert_allclose(values[SPECIFIC_COLUMNS[0]],
                                       values[SPECIFIC_COLUMNS[1]] + values[SPECIFIC_COLUMNS[2]],
                                       rtol=0, atol=1e-12)

    def test_nonabsorbing_case_is_zero_within_roundoff(self):
        values = self.results[(1.5, 0.)]
        self.assertLess(np.max(np.abs(values["qabs"])), 1e-13)
        self.assertLess(np.max(np.abs(values[SPECIFIC_COLUMNS[2]])), 1e-12)

    def test_rayleigh_absorber_specific_absorption_plateau(self):
        m = complex(2., -0.66)
        polarizability = (m**2 - 1) / (m**2 + 2)
        expected = -6 * np.pi * polarizability.imag / 0.55
        values = self.results[(2., 0.66)]
        self.assertAlmostEqual(values[SPECIFIC_COLUMNS[2]][0] / expected, 1., delta=2e-4)

    def test_rayleigh_nonabsorber_specific_extinction_scales_as_radius_cubed(self):
        radii = np.array([0.001, 0.002])
        values = mie_sweep(2 * np.pi * radii / 0.55, 1.5, 0.)
        specific = mass_specific_cross_section(values["qext"], radii, 1.)
        self.assertAlmostEqual(specific[1] / specific[0], 8., delta=0.01)
        self.assertTrue(np.all(np.abs(values["g"]) < 0.001))

    @unittest.skipIf(reference_efficiencies is None, "SciPy reference requires scipy")
    def test_independent_reference_and_early_termination_regression(self):
        radii = np.array([0.001, 0.018706821403658005, 0.1, 0.55, 1., 8.554, 10.])
        x = 2 * np.pi * radii / 0.55
        for n, k in ((2., 0.66), (1.5, 0.)):
            actual = mie_sweep(x, n, k)
            expected = np.array([reference_efficiencies(value, n, k) for value in x])
            np.testing.assert_allclose(
                np.column_stack([actual["qext"], actual["qscat"], actual["g"]]),
                expected, rtol=0, atol=1e-8,
            )

    def test_question2_saved_curve_compatibility(self):
        saved = np.genfromtxt(ROOT / "HW5/results/question2_mie_results.csv",
                              delimiter=",", names=True)
        for n, k in ((1.3873, 0.0066), (1.2993, 0.037)):
            rows = saved[(saved["real_index"] == n) & (saved["absorption_index"] == k)]
            actual = mie_sweep(rows["x"], n, k)
            for name in ("qext", "qscat", "g"):
                np.testing.assert_allclose(actual[name], rows[name], rtol=0, atol=5e-7)


if __name__ == "__main__":
    unittest.main()
