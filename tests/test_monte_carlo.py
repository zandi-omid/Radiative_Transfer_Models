"""Checks for the shared photon transport model."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from models.monte_carlo import MonteCarloTransport


class MonteCarloTransportTests(unittest.TestCase):
    def test_conservation_and_scattering(self):
        result = MonteCarloTransport(workers=2).run(1000, seed=101)
        self.assertAlmostEqual(result.reflectivity + result.transmissivity, 1.0)
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


if __name__ == "__main__":
    unittest.main()
