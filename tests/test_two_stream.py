"""Validation tests for the Meador-Weaver two-stream implementation."""

import unittest

import numpy as np

from models.two_stream import METHODS, coefficients, reflect_transmit
from HW3.validate_question1 import reference_coefficients, reference_reflect_transmit


class TwoStreamTests(unittest.TestCase):
    def setUp(self):
        self.tau = np.linspace(0.0, 10.0, 201)

    def test_zero_optical_depth(self):
        for omega0 in (1.0, 0.9):
            for mu0 in (0.9, 0.5, 0.1):
                for method in METHODS:
                    reflected, transmitted = reflect_transmit(
                        0.0,
                        omega0,
                        mu0,
                        0.75,
                        method,
                    )
                    self.assertAlmostEqual(reflected, 0.0, places=12)
                    self.assertAlmostEqual(transmitted, 1.0, places=12)

    def test_conservative_energy_balance(self):
        for mu0 in (0.9, 0.5, 0.1):
            for method in METHODS:
                reflected, transmitted = reflect_transmit(
                    self.tau,
                    1.0,
                    mu0,
                    0.75,
                    method,
                )
                np.testing.assert_allclose(reflected + transmitted, 1.0, atol=1e-12)

    def test_absorbing_energy_partition(self):
        for mu0 in (0.9, 0.5, 0.1):
            for method in METHODS:
                reflected, transmitted = reflect_transmit(
                    self.tau,
                    0.9,
                    mu0,
                    0.75,
                    method,
                )
                absorbed = 1.0 - reflected - transmitted
                np.testing.assert_allclose(
                    reflected + transmitted + absorbed,
                    1.0,
                    atol=1e-12,
                )
                self.assertGreater(float(absorbed[-1]), 0.0)

    def test_stable_form_matches_direct_equations(self):
        omega0 = 0.9
        mu0 = 0.5
        tau = 1.8
        for method in METHODS:
            gamma1, _, gamma3, gamma4, alpha1, alpha2, k = coefficients(
                omega0,
                mu0,
                0.75,
                method,
            )
            d_value = (k + gamma1) * np.exp(k * tau) + (
                k - gamma1
            ) * np.exp(-k * tau)
            denominator = (1.0 - (k * mu0) ** 2) * d_value
            reflected_direct = omega0 / denominator * (
                (1.0 - k * mu0) * (alpha2 + k * gamma3) * np.exp(k * tau)
                - (1.0 + k * mu0)
                * (alpha2 - k * gamma3)
                * np.exp(-k * tau)
                - 2.0 * k * (gamma3 - alpha2 * mu0) * np.exp(-tau / mu0)
            )
            transmitted_direct = np.exp(-tau / mu0) * (
                1.0
                - omega0
                / denominator
                * (
                    (1.0 + k * mu0)
                    * (alpha1 + k * gamma4)
                    * np.exp(k * tau)
                    - (1.0 - k * mu0)
                    * (alpha1 - k * gamma4)
                    * np.exp(-k * tau)
                    - 2.0
                    * k
                    * (gamma4 + alpha1 * mu0)
                    * np.exp(tau / mu0)
                )
            )
            reflected, transmitted = reflect_transmit(
                tau,
                omega0,
                mu0,
                0.75,
                method,
            )
            self.assertAlmostEqual(reflected, reflected_direct, places=12)
            self.assertAlmostEqual(transmitted, transmitted_direct, places=12)

    def test_requested_sweep_is_finite(self):
        combinations = 0
        for omega0 in (1.0, 0.9):
            for mu0 in (0.9, 0.5, 0.1):
                for method in METHODS:
                    reflected, transmitted = reflect_transmit(
                        self.tau,
                        omega0,
                        mu0,
                        0.75,
                        method,
                    )
                    self.assertTrue(np.all(np.isfinite(reflected)))
                    self.assertTrue(np.all(np.isfinite(transmitted)))
                    combinations += 1
        self.assertEqual(combinations, 12)

    def test_table_1_coefficients_against_independent_values(self):
        expected = {
            "eddington": (0.34375, 0.14375, 0.21875, 0.78125),
            "quadrature": (
                0.3680607966083864,
                0.1948557158514987,
                0.1752404735808355,
                0.8247595264191645,
            ),
        }
        for method in METHODS:
            actual = coefficients(0.9, 0.5, 0.75, method)[:4]
            np.testing.assert_allclose(actual, expected[method], atol=1.0e-15)
            np.testing.assert_allclose(
                reference_coefficients(0.9, 0.5, 0.75, method),
                expected[method],
                atol=1.0e-15,
            )

    def test_analytical_solution_matches_independent_ode_solver(self):
        depths = np.linspace(0.0, 10.0, 41)
        for omega0 in (1.0, 0.9):
            for mu0 in (0.9, 0.5, 0.1):
                for method in METHODS:
                    analytical = reflect_transmit(depths, omega0, mu0, 0.75, method)
                    numerical = reference_reflect_transmit(
                        depths, omega0, mu0, 0.75, method
                    )
                    np.testing.assert_allclose(analytical[0], numerical[0], atol=2e-9)
                    np.testing.assert_allclose(analytical[1], numerical[1], atol=2e-9)

    def test_thin_layer_negative_reflectivity_is_reproduced_by_ode_solver(self):
        depth = 0.05
        for method in METHODS:
            analytical = reflect_transmit(depth, 1.0, 0.9, 0.75, method)
            numerical = reference_reflect_transmit(depth, 1.0, 0.9, 0.75, method)
            self.assertLess(analytical[0], 0.0)
            self.assertGreater(analytical[1], 1.0)
            np.testing.assert_allclose(analytical, numerical, atol=2e-10)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            reflect_transmit(-1.0)
        with self.assertRaises(ValueError):
            reflect_transmit(1.0, omega0=1.1)
        with self.assertRaises(ValueError):
            reflect_transmit(1.0, mu0=-0.5)
        with self.assertRaises(ValueError):
            reflect_transmit(1.0, method="unknown")


if __name__ == "__main__":
    unittest.main()
