"""Photon transport in a nonabsorbing, isotropically scattering layer."""

from dataclasses import dataclass
from multiprocessing import Pool, cpu_count
from time import perf_counter

import numpy as np


@dataclass(frozen=True)
class SimulationResult:
    reflectivity: float
    transmissivity: float
    mean_scatterings: float
    elapsed_seconds: float


def _simulate_chunk(args):
    n_photons, tau, mu0, l0, seed_sequence = args
    rng = np.random.default_rng(seed_sequence)
    z_top = tau * l0
    reflected = transmitted = scatterings = 0

    for _ in range(n_photons):
        z = z_top
        mu = mu0
        while True:
            z += mu * (-l0 * np.log(1.0 - rng.random()))
            if z > z_top:
                reflected += 1
                break
            if z < 0.0:
                transmitted += 1
                break
            scatterings += 1
            mu = rng.uniform(-1.0, 1.0)

    return reflected, transmitted, scatterings


class MonteCarloTransport:
    """Run independent photon histories with reproducible worker streams."""

    def __init__(self, tau=4.0, mu0=-0.7, l0=1.0, workers=24):
        if tau <= 0 or l0 <= 0:
            raise ValueError("tau and l0 must be positive")
        if not -1.0 <= mu0 < 0.0:
            raise ValueError("mu0 must be in [-1, 0)")
        if workers < 1:
            raise ValueError("workers must be positive")
        self.tau = tau
        self.mu0 = mu0
        self.l0 = l0
        self.workers = min(workers, cpu_count())

    def run(self, n_photons, seed=42):
        if n_photons < 1:
            raise ValueError("n_photons must be positive")
        n_workers = min(self.workers, n_photons)
        chunk_sizes = np.full(n_workers, n_photons // n_workers, dtype=int)
        chunk_sizes[: n_photons % n_workers] += 1
        seeds = np.random.SeedSequence(seed).spawn(n_workers)
        args = [
            (int(size), self.tau, self.mu0, self.l0, stream)
            for size, stream in zip(chunk_sizes, seeds)
        ]

        start = perf_counter()
        if n_workers == 1:
            counts = [_simulate_chunk(args[0])]
        else:
            with Pool(processes=n_workers) as pool:
                counts = pool.map(_simulate_chunk, args)
        elapsed = perf_counter() - start

        reflected = sum(count[0] for count in counts)
        transmitted = sum(count[1] for count in counts)
        scatterings = sum(count[2] for count in counts)
        return SimulationResult(
            reflected / n_photons,
            transmitted / n_photons,
            scatterings / n_photons,
            elapsed,
        )
