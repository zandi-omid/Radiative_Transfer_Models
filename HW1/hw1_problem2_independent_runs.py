#!/usr/bin/env python3
"""HW1 Problem 2: independent Monte Carlo runs."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from models.monte_carlo import MonteCarloTransport


def main():

    # --------------------------------------------------------
    # Problem 2 parameters
    # --------------------------------------------------------

    tau = 4.0
    mu0 = -0.7
    l0 = 1.0

    # Adequate photon number determined from Problem 1
    n_photons = 100_000

    # Homework asks for 5-10 iterations
    n_runs = 10

    # Use up to 24 CPUs
    model = MonteCarloTransport(tau=tau, mu0=mu0, l0=l0)
    n_workers = model.workers

    print(f"Using {n_workers} worker processes")

    # --------------------------------------------------------
    # Different seed for every independent run
    # --------------------------------------------------------

    seeds = [
        101,
        202,
        303,
        404,
        505,
        606,
        707,
        808,
        909,
        1010,
    ]

    # Storage
    R_values = []
    T_values = []

    print()
    print("HW1 - Problem 2")
    print("----------------")
    print(f"tau       = {tau}")
    print(f"mu0       = {mu0}")
    print(f"N photons = {n_photons:,}")
    print(f"N runs    = {n_runs}")
    print()

    print(
        f"{'Run':>6} "
        f"{'Seed':>8} "
        f"{'R':>12} "
        f"{'T':>12} "
        f"{'R+T':>12}"
    )

    print("-" * 55)

    # --------------------------------------------------------
    # Perform independent Monte Carlo runs
    # --------------------------------------------------------

    for i in range(n_runs):

        result = model.run(n_photons, seed=seeds[i])
        R, T = result.reflectivity, result.transmissivity

        R_values.append(R)
        T_values.append(T)

        print(
            f"{i+1:6d} "
            f"{seeds[i]:8d} "
            f"{R:12.6f} "
            f"{T:12.6f} "
            f"{R+T:12.6f}"
        )

    # Convert to arrays
    R_values = np.array(R_values)
    T_values = np.array(T_values)

    # --------------------------------------------------------
    # Calculate statistics
    # --------------------------------------------------------

    mean_R = np.mean(R_values)
    mean_T = np.mean(T_values)

    # ddof=1 -> sample standard deviation
    std_R = np.std(R_values, ddof=1)
    std_T = np.std(T_values, ddof=1)

    # Theoretical Monte Carlo standard error for comparison
    theoretical_sigma_R = np.sqrt(
        mean_R * (1.0 - mean_R) / n_photons
    )

    theoretical_sigma_T = np.sqrt(
        mean_T * (1.0 - mean_T) / n_photons
    )

    # --------------------------------------------------------
    # Print summary
    # --------------------------------------------------------

    print()
    print("=" * 55)
    print("Summary")
    print("=" * 55)

    print(f"Mean Reflectivity       = {mean_R:.6f}")
    print(f"Std. Dev. Reflectivity  = {std_R:.6f}")

    print()

    print(f"Mean Transmissivity      = {mean_T:.6f}")
    print(f"Std. Dev. Transmissivity = {std_T:.6f}")

    print()
    print("Theoretical standard error:")
    print(f"sigma_R = {theoretical_sigma_R:.6f}")
    print(f"sigma_T = {theoretical_sigma_T:.6f}")

    print()
    print(f"Mean R + Mean T = {mean_R + mean_T:.6f}")

    # ========================================================
    # Plot
    # ========================================================

    output_dir = Path(__file__).resolve().parent

    runs = np.arange(1, n_runs + 1)

    plt.figure(figsize=(9, 6))

    plt.plot(
        runs,
        R_values,
        marker="o",
        linewidth=2,
        label="Reflectivity (R)",
    )

    plt.plot(
        runs,
        T_values,
        marker="s",
        linewidth=2,
        label="Transmissivity (T)",
    )

    # Mean values
    plt.axhline(
        mean_R,
        linestyle="--",
        label=f"Mean R = {mean_R:.4f}",
    )

    plt.axhline(
        mean_T,
        linestyle="--",
        label=f"Mean T = {mean_T:.4f}",
    )

    plt.xlabel("Independent Monte Carlo Run")
    plt.ylabel("Reflectivity / Transmissivity")

    plt.title(
        r"Monte Carlo Run-to-Run Variability: "
        r"$\tau=4$, $\mu_0=-0.7$, $N=100{,}000$"
    )

    plt.xticks(runs)

    plt.grid(True, alpha=0.3)
    plt.legend()

    plt.tight_layout()

    plt.savefig(
        output_dir / "problem2_independent_runs.png",
        dpi=300,
    )

    plt.show()


if __name__ == "__main__":
    main()
