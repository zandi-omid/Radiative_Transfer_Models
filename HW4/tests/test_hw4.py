"""Validate infinite reflections against finite sums and an independent ODE."""

import unittest

import numpy as np

from models.liou_two_stream import liou_reflect_transmit
from models.surface_two_stream import diffuse_reflect_transmit, surface_transport


def reflecting_boundary_reference(tau, albedo, omega0, mu0=0.9, g=0.75):
    """RK4 solution of the discrete-ordinate equations with a reflecting bottom.

    Integrate homogeneous and beam-forced intensities, plus their volume
    absorption integrals. Shoot directly for U(tau)=As*(D(tau)+beam).
    No layer R/T or geometric-series coefficients are used here.
    """
    mu = 1.0 / np.sqrt(3.0)
    phase_same = omega0 * (1.0 + 3.0 * g * mu**2)
    phase_opposite = omega0 * (1.0 - 3.0 * g * mu**2)
    matrix = np.array([
        [1.0 - phase_same / 2.0, -phase_opposite / 2.0],
        [phase_opposite / 2.0, -1.0 + phase_same / 2.0],
    ]) / mu
    # Upward scattering from the downward solar beam has the opposite sign.
    source = omega0 * np.array([
        -(1.0 - 3.0 * g * mu * mu0),
        1.0 + 3.0 * g * mu * mu0,
    ]) / (2.0 * mu0)

    def derivative(depth, state):
        beam = np.exp(-depth / mu0)
        h = matrix @ state[:2]
        p = matrix @ state[2:4] + source * beam
        h_abs = (1.0 - omega0) * np.sum(state[:2]) / mu
        p_abs = (1.0 - omega0) * (np.sum(state[2:4]) / mu + beam / mu0)
        return np.array([*h, *p, h_abs, p_abs])

    state = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    steps = max(1, int(np.ceil(tau / 0.002)))
    step = tau / steps
    for index in range(steps):
        depth = index * step
        k1 = derivative(depth, state)
        k2 = derivative(depth + step / 2.0, state + step * k1 / 2.0)
        k3 = derivative(depth + step / 2.0, state + step * k2 / 2.0)
        k4 = derivative(depth + step, state + step * k3)
        state += step * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    h_up, h_down, p_up, p_down, h_abs, p_abs = state
    beam = np.exp(-tau / mu0)
    reflected = (albedo * (p_down + beam) - p_up) / (h_up - albedo * h_down)
    downward = p_down + reflected * h_down + beam
    absorption = p_abs + reflected * h_abs
    return reflected, downward, (1.0 - albedo) * downward, absorption


class HW4Tests(unittest.TestCase):
    def test_black_surface_recovers_hw3(self):
        for tau in (1.0, 4.0):
            for omega0 in (0.99999, 0.8):
                result = surface_transport(tau, 0.0, omega0)
                r, t = liou_reflect_transmit(tau, omega0, 0.9, 0.75)
                self.assertAlmostEqual(result.toa_reflectivity, r, places=13)
                self.assertAlmostEqual(result.surface_net_flux, t, places=13)

    def test_perfect_surface_and_exact_conservative_limit(self):
        for tau in (0.0, 1.0, 4.0):
            result = surface_transport(tau, 1.0, 1.0)
            self.assertEqual(result.surface_net_flux, 0.0)
            self.assertAlmostEqual(result.toa_reflectivity, 1.0, places=13)
        for tau in (1.0, 4.0):
            result = surface_transport(tau, 1.0, 0.8)
            self.assertLess(result.toa_reflectivity, 1.0)
            self.assertGreater(result.atmospheric_absorption, 0.0)

    def test_energy_balance_and_flux_normalization(self):
        for tau in (1.0, 4.0):
            for omega0 in (1.0, 0.99999, 0.8):
                result = surface_transport(tau, np.linspace(0, 1, 11), omega0)
                np.testing.assert_allclose(
                    result.toa_reflectivity + result.surface_net_flux
                    + result.atmospheric_absorption, 1.0, rtol=0.0, atol=2e-13,
                )
                self.assertTrue(np.all(result.reflection_denominator > 0.0))
                self.assertTrue(np.all(np.isfinite(result.surface_downward_flux)))
                scaled = surface_transport(tau, 0.5, omega0, incident_flux=0.9)
                unit = surface_transport(tau, 0.5, omega0)
                self.assertAlmostEqual(scaled.toa_reflectivity, unit.toa_reflectivity)
                self.assertAlmostEqual(scaled.surface_net_flux, 0.9 * unit.surface_net_flux)

    def test_finite_reflection_series(self):
        for tau in (1.0, 4.0):
            for omega0 in (0.99999, 0.8):
                r, t = liou_reflect_transmit(tau, omega0, 0.9, 0.75)
                rd, td = diffuse_reflect_transmit(tau, omega0, 0.75)
                for albedo in (0.0, 0.5, 1.0):
                    downward, reflected, atmospheric = 0.0, r, 1.0 - r - t
                    arriving = t
                    for _ in range(200):
                        downward += arriving
                        upward = albedo * arriving
                        reflected += td * upward
                        atmospheric += (1.0 - rd - td) * upward
                        arriving = rd * upward
                    result = surface_transport(tau, albedo, omega0)
                    np.testing.assert_allclose(
                        [result.toa_reflectivity, result.surface_downward_flux,
                         result.atmospheric_absorption],
                        [reflected, downward, atmospheric], rtol=0.0, atol=2e-13,
                    )

    def test_independent_reflecting_boundary_and_volume_absorption(self):
        for tau in (1.0, 4.0):
            for omega0 in (1.0, 0.99999, 0.8):
                for albedo in (0.0, 0.5, 1.0):
                    with self.subTest(tau=tau, omega0=omega0, albedo=albedo):
                        result = surface_transport(tau, albedo, omega0)
                        reference = reflecting_boundary_reference(tau, albedo, omega0)
                        np.testing.assert_allclose(
                            [result.toa_reflectivity, result.surface_downward_flux,
                             result.surface_net_flux, result.atmospheric_absorption],
                            reference, rtol=0.0, atol=2e-10,
                        )

    def test_diffuse_limits_and_no_atmosphere(self):
        for omega0 in (0.0, 0.8, 0.99999, 1.0):
            np.testing.assert_allclose(diffuse_reflect_transmit(0.0, omega0), [0, 1])
        r, t = diffuse_reflect_transmit(4.0, 1.0)
        self.assertAlmostEqual(r + t, 1.0, places=14)
        albedos = np.linspace(0, 1, 11)
        result = surface_transport(0.0, albedos, 1.0)
        np.testing.assert_allclose(result.toa_reflectivity, albedos)
        np.testing.assert_allclose(result.surface_net_flux, 1.0 - albedos)
        # Pure absorption: no layer reflection, including a reflected surface beam.
        result = surface_transport(1.0, 0.5, 0.0)
        self.assertAlmostEqual(
            result.toa_reflectivity,
            0.5 * np.exp(-1.0 / 0.9) * np.exp(-np.sqrt(3.0)),
        )

    def test_invalid_inputs(self):
        for kwargs in (
            {"tau": -1}, {"surface_albedo": 1.1}, {"omega0": 1.1},
            {"mu0": 0}, {"g": np.nan}, {"incident_flux": 0},
        ):
            options = {"tau": 1.0, "surface_albedo": 0.5}
            options.update(kwargs)
            with self.assertRaises(ValueError):
                surface_transport(**options)


if __name__ == "__main__":
    unittest.main()
