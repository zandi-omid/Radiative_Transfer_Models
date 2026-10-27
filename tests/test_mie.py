"""Focused checks of the course Mie code and its batch interface."""

from pathlib import Path
import shutil
import unittest

import numpy as np

from models.mie import mie_sweep

try:
    from HW5.reference_mie import reference_efficiencies
except ImportError:
    reference_efficiencies = None


@unittest.skipUnless(shutil.which("gfortran"), "gfortran is required")
class MieTests(unittest.TestCase):
    def test_supplied_subroutine_is_preserved(self):
        root = Path(__file__).resolve().parents[1] / "HW5"
        original = (root / "reference/2026_HW5_F90.f.txt").read_text()
        supplied = original[original.index("      SUBROUTINE DMIESS("):].rstrip()
        self.assertEqual((root / "fortran/dmiess.f").read_text().rstrip(), supplied)

    def test_complete_question2_grid(self):
        x = np.arange(1, 1001) * 0.05
        for n, k in ((1.2993, 0.037), (1.3873, 0.0066)):
            result = mie_sweep(x, n, k)
            self.assertEqual(len(result["x"]), 1000)
            self.assertTrue(all(np.all(np.isfinite(v)) for v in result.values()))
            self.assertTrue(np.all(result["qabs"] >= 0))
            self.assertTrue(np.all((result["g"] >= -1) & (result["g"] <= 1)))

    @unittest.skipIf(reference_efficiencies is None, "SciPy reference requires scipy")
    def test_independent_bessel_reference(self):
        x = np.array([0.05, 0.1, 0.25, 1., 3., 5., 10., 25., 50.])
        for n, k in ((1.2993, 0.037), (1.3873, 0.0066)):
            actual = mie_sweep(x, n, k)
            expected = np.array([reference_efficiencies(v, n, k) for v in x])
            np.testing.assert_allclose(
                np.column_stack([actual["qext"], actual["qscat"], actual["g"]]),
                expected, rtol=0, atol=5e-7,
            )

    def test_nonabsorbing_control(self):
        result = mie_sweep([0.05, 1., 10., 50.], 1.5, 0.)
        np.testing.assert_allclose(result["qext"], result["qscat"], rtol=0, atol=1e-12)

    def test_fixed_size_parameter_is_wavelength_invariant(self):
        first = mie_sweep([1., 10.], 1.2993, 0.037, wavelength_um=0.55)
        second = mie_sweep([1., 10.], 1.2993, 0.037, wavelength_um=1.1)
        for key in ("qext", "qscat", "g"):
            np.testing.assert_allclose(first[key], second[key], rtol=0, atol=1e-12)
        np.testing.assert_allclose(second["radius_um"], 2 * first["radius_um"])

    def test_invalid_inputs(self):
        for x, n, k in (([0], 1.3, 0.1), ([-1], 1.3, 0),
                        ([np.nan], 1.3, 0), ([], 1.3, 0),
                        ([1], 0, 0), ([1], 1.3, -0.1)):
            with self.assertRaises(ValueError):
                mie_sweep(x, n, k)


if __name__ == "__main__":
    unittest.main()
