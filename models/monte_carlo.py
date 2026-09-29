"""Monte Carlo photon transport in a nonabsorbing plane-parallel layer."""

from dataclasses import dataclass
from multiprocessing import Pool, cpu_count
from time import perf_counter

import numpy as np


@dataclass(frozen=True)
class SimulationResult:
    reflectivity: float
    transmissivity: float
    absorptivity: float
    mean_scatterings: float
    elapsed_seconds: float


def sample_hg_cos_theta(rng, g):
    """Sample the cosine of an HG scattering angle by inverse CDF."""
    xi = rng.random()
    if abs(g) < 1.0e-12:
        return 2.0 * xi - 1.0

    ratio = (1.0 - g * g) / (1.0 - g + 2.0 * g * xi)
    cos_theta = (1.0 + g * g - ratio * ratio) / (2.0 * g)
    return float(np.clip(cos_theta, -1.0, 1.0))


def scatter_direction(mu_old, rng, g):
    """Rotate an HG scattering direction into the vertical coordinate system.

    ``cos_theta`` is relative to the photon's current direction, whereas
    ``mu_old`` and the returned value are direction cosines relative to the
    vertical. A random azimuth completes the three-dimensional rotation.
    """
    if abs(g) < 1.0e-12:
        return rng.uniform(-1.0, 1.0)

    cos_theta = sample_hg_cos_theta(rng, g)
    sin_theta = np.sqrt(max(0.0, 1.0 - cos_theta * cos_theta))
    sin_old = np.sqrt(max(0.0, 1.0 - mu_old * mu_old))
    phi = rng.uniform(0.0, 2.0 * np.pi)
    mu_new = mu_old * cos_theta + sin_old * sin_theta * np.cos(phi)
    return float(np.clip(mu_new, -1.0, 1.0))


def _simulate_chunk(args):
    n_photons, tau, mu0, l0, g, omega0, seed_sequence = args
    rng = np.random.default_rng(seed_sequence)
    z_top = tau * l0
    reflected = transmitted = absorbed = scatterings = 0

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

            # omega0 is the probability of scattering at an interaction.
            # Avoid drawing an extra random number when omega0=1 so the
            # validated nonabsorbing simulations retain their exact sequence.
            if omega0 < 1.0 and rng.random() >= omega0:
                absorbed += 1
                break
            scatterings += 1
            mu = scatter_direction(mu, rng, g)

    return reflected, transmitted, absorbed, scatterings


class MonteCarloTransport:
    """Run independent photon histories with reproducible worker streams."""

    def __init__(
        self,
        tau=4.0,
        mu0=-0.7,
        l0=1.0,
        g=0.0,
        omega0=1.0,
        workers=24,
    ):
        if tau <= 0 or l0 <= 0:
            raise ValueError("tau and l0 must be positive")
        if not -1.0 <= mu0 < 0.0:
            raise ValueError("mu0 must be in [-1, 0)")
        if not -1.0 < g < 1.0:
            raise ValueError("g must be between -1 and 1")
        if not 0.0 <= omega0 <= 1.0:
            raise ValueError("omega0 must be between 0 and 1")
        if workers < 1:
            raise ValueError("workers must be positive")
        self.tau = tau
        self.mu0 = mu0
        self.l0 = l0
        self.g = g
        self.omega0 = omega0
        self.workers = min(workers, cpu_count())

    def run(self, n_photons, seed=42):
        if n_photons < 1:
            raise ValueError("n_photons must be positive")
        n_workers = min(self.workers, n_photons)
        chunk_sizes = np.full(n_workers, n_photons // n_workers, dtype=int)
        chunk_sizes[: n_photons % n_workers] += 1
        seeds = np.random.SeedSequence(seed).spawn(n_workers)
        args = [
            (int(size), self.tau, self.mu0, self.l0, self.g, self.omega0, stream)
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
        absorbed = sum(count[2] for count in counts)
        scatterings = sum(count[3] for count in counts)
        return SimulationResult(
            reflectivity=reflected / n_photons,
            transmissivity=transmitted / n_photons,
            absorptivity=absorbed / n_photons,
            mean_scatterings=scatterings / n_photons,
            elapsed_seconds=elapsed,
        )
