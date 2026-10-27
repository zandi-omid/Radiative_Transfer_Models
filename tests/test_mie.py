"""Focused checks of the course Mie code and its batch interface."""

from pathlib import Path
import shutil
import unittest

import numpy as np

from models.mie import mass_specific_cross_section, mie_sweep

try:
    from HW5.reference_mie import reference_efficiencies
except ImportError:
    reference_efficiencies = None


@unittest.skipUnless(shutil.which("gfortran"), "gfortran is required")
class MieTests(unittest.TestCase):
    def test_only_documented_numerical_safeguards(self):
        root = Path(__file__).resolve().parents[1] / "HW5"
        original = (root / "reference/2026_HW5_F90.f.txt").read_text()
        supplied = original[original.index("      SUBROUTINE DMIESS("):].rstrip()
        expected = supplied.replace(
            "      IFLAG = 1\n",
            "      IFLAG = 1\n"
            "C     Equal core/shell material with equal radii is exactly homogeneous.\n"
            "C     Avoid removable singularities in the unnecessary coated-sphere algebra.\n"
            "      IF ( R .EQ. RO .AND. RE2 .EQ. RFR .AND.\n"
            "     1     TMAG2 .EQ. RFI ) IFLAG = 2\n",
            1,
        ).replace(
            "      IF ( N .LT. M )   GO TO 1002\n"
            "      IF ( IFLAG .EQ. 2 )   GO TO 1001\n",
            "      IF ( IFLAG .EQ. 2 )   GO TO 1001\n"
            "      IF ( N .LT. M )   GO TO 1002\n",
            1,
        )
        expected = expected.replace(
            "      X  =  RO * WVNO\n",
            "      X  =  RO * WVNO\n"
            "C     Do not terminate at an isolated small coefficient before the tail.\n"
            "      NSTOP = CEILING( X + 4.0D0 * X**(1.0D0/3.0D0) + 2.0D0 )\n",
            1,
        ).replace(
            "      IF ( T(4) .LT. 1.0E-10 )   GO TO 100\n",
            "      IF ( T(4) .LT. 1.0E-10 .AND. N .GE. NSTOP ) GO TO 100\n",
            1,
        )
        expected = expected.replace(
            "      T(1)      =  ( 3.14159265359 * THETD(J) ) / 180.0\n",
            "      T(1)      =  ( DACOS(-1.0D0) * THETD(J) ) / 180.0D0\n",
            1,
        )
        self.assertEqual((root / "fortran/dmiess.f").read_text().rstrip(), expected)

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

    @unittest.skipIf(reference_efficiencies is None, "SciPy reference requires scipy")
    def test_equal_material_core_shell_removable_singularity(self):
        x = np.array([2 * np.pi, 4 * np.pi, 20 * np.pi])
        actual = mie_sweep(x, 1.5, 0.)
        expected = np.array([reference_efficiencies(value, 1.5, 0.) for value in x])
        np.testing.assert_allclose(
            np.column_stack([actual["qext"], actual["qscat"], actual["g"]]),
            expected, rtol=0, atol=5e-7,
        )

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
