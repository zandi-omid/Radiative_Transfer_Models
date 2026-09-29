"""Checks for the shared photon transport model."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from models.monte_carlo import (
    MonteCarloTransport,
    sample_hg_cos_theta,
    scatter_direction,
)


class MonteCarloTransportTests(unittest.TestCase):
    def test_conservation_and_scattering(self):
        result = MonteCarloTransport(workers=2).run(1000, seed=101)
        self.assertAlmostEqual(
            result.reflectivity + result.transmissivity + result.absorptivity,
            1.0,
        )
        self.assertEqual(result.absorptivity, 0.0)
        self.assertGreater(result.mean_scatterings, 0)
        self.assertGreaterEqual(result.elapsed_seconds, 0)

    def test_reproducible_with_same_worker_count(self):
        model = MonteCarloTransport(workers=2)
        first = model.run(1000, seed=42)
        second = model.run(1000, seed=42)
        self.assertEqual(first.reflectivity, second.reflectivity)
        self.assertEqual(first.transmissivity, second.transmissivity)
        self.assertEqual(first.mean_scatterings, second.mean_scatterings)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            MonteCarloTransport(tau=0)
        with self.assertRaises(ValueError):
            MonteCarloTransport(mu0=0)
        with self.assertRaises(ValueError):
            MonteCarloTransport().run(0)
        with self.assertRaises(ValueError):
            MonteCarloTransport(g=1.0)
        with self.assertRaises(ValueError):
            MonteCarloTransport(omega0=1.1)

    def test_hg_samples_have_expected_mean_cosine(self):
        rng = np.random.default_rng(123)
        samples = np.array([sample_hg_cos_theta(rng, 0.75) for _ in range(50_000)])
        self.assertAlmostEqual(float(samples.mean()), 0.75, delta=0.01)
        self.assertTrue(np.all((-1.0 <= samples) & (samples <= 1.0)))

    def test_hg_rotation_preserves_direction_on_average(self):
        rng = np.random.default_rng(456)
        mu_old = -0.7
        samples = np.array(
            [scatter_direction(mu_old, rng, 0.75) for _ in range(50_000)]
        )
        self.assertAlmostEqual(float(samples.mean()), mu_old * 0.75, delta=0.01)

    def test_hg_transport_conserves_photons(self):
        result = MonteCarloTransport(g=0.75, workers=2).run(2_000, seed=42)
        self.assertAlmostEqual(result.reflectivity + result.transmissivity, 1.0)

    def test_absorbing_transport_conserves_photons(self):
        result = MonteCarloTransport(g=0.75, omega0=0.85, workers=2).run(
            2_000,
            seed=42,
        )
        self.assertAlmostEqual(
            result.reflectivity + result.transmissivity + result.absorptivity,
            1.0,
        )
        self.assertGreater(result.absorptivity, 0.0)


if __name__ == "__main__":
    unittest.main()
