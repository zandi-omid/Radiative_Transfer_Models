"""Tests for the Liou (1974) two-stream discrete-ordinate model."""

import unittest

import numpy as np

from models.liou_two_stream import liou_coefficients, liou_reflect_transmit
from models.two_stream import reflect_transmit
from HW3.validate_question1 import reference_reflect_transmit


class LiouTwoStreamTests(unittest.TestCase):
    def setUp(self):
        self.depths = np.linspace(0.0, 10.0, 201)

    def test_zero_optical_depth(self):
        for omega0 in (1.0, 0.9):
            for mu0 in (0.9, 0.5, 0.1):
                reflected, transmitted = liou_reflect_transmit(
                    0.0, omega0, mu0, 0.75
                )
                self.assertAlmostEqual(reflected, 0.0, places=14)
                self.assertAlmostEqual(transmitted, 1.0, places=14)

    def test_requested_sweep_is_finite(self):
        for omega0 in (1.0, 0.9):
            for mu0 in (0.9, 0.5, 0.1):
                reflected, transmitted = liou_reflect_transmit(
                    self.depths, omega0, mu0, 0.75
                )
                self.assertTrue(np.all(np.isfinite(reflected)))
                self.assertTrue(np.all(np.isfinite(transmitted)))

    def test_conservative_energy_balance(self):
        for mu0 in (0.9, 0.5, 0.1):
            reflected, transmitted = liou_reflect_transmit(
                self.depths, 1.0, mu0, 0.75
            )
            np.testing.assert_allclose(reflected + transmitted, 1.0, atol=1e-14)

    def test_reference_values_at_tau_four(self):
        expected = {
            (1.0, 0.5): (0.4999879574334767, 0.5000120425665233),
            (0.9, 0.5): (0.3003520895659419, 0.2429683466722991),
            (1.0, 0.1): (0.6856406460551019, 0.3143593539448981),
            (0.9, 0.1): (0.5092556730499518, 0.1418884771959186),
        }
        for (omega0, mu0), reference in expected.items():
            actual = liou_reflect_transmit(4.0, omega0, mu0, 0.75)
            np.testing.assert_allclose(actual, reference, atol=5e-15)

    def test_matches_meador_weaver_quadrature(self):
        for omega0 in (1.0, 0.9):
            for mu0 in (0.9, 0.5, 0.1):
                liou = liou_reflect_transmit(self.depths, omega0, mu0, 0.75)
                quadrature = reflect_transmit(
                    self.depths, omega0, mu0, 0.75, method="quadrature"
                )
                np.testing.assert_allclose(liou[0], quadrature[0], atol=1e-14)
                np.testing.assert_allclose(liou[1], quadrature[1], atol=1e-14)

    def test_matches_independent_ode_reference_cases(self):
        depths = np.array([0.0, 0.2, 1.0, 4.0, 10.0])
        for omega0, mu0 in ((1.0, 0.9), (1.0, 0.1), (0.9, 0.5)):
            liou = liou_reflect_transmit(depths, omega0, mu0, 0.75)
            numerical = reference_reflect_transmit(
                depths, omega0, mu0, 0.75, method="quadrature"
            )
            np.testing.assert_allclose(liou[0], numerical[0], atol=2e-9)
            np.testing.assert_allclose(liou[1], numerical[1], atol=2e-9)

    def test_phase_moments_and_eigenvalue(self):
        omega_bar_0, omega_bar_1, eigenvalue = liou_coefficients(0.9, 0.75)
        self.assertAlmostEqual(omega_bar_0, 0.9)
        self.assertAlmostEqual(omega_bar_1, 2.025)
        self.assertAlmostEqual(eigenvalue, np.sqrt(0.3 * 0.325))

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            liou_reflect_transmit(-1.0)
        with self.assertRaises(ValueError):
            liou_reflect_transmit(1.0, omega0=1.1)
        with self.assertRaises(ValueError):
            liou_reflect_transmit(1.0, mu0=0.0)
        with self.assertRaises(ValueError):
            liou_reflect_transmit(1.0, g=1.1)


if __name__ == "__main__":
    unittest.main()
