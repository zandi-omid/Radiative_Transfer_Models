#!/usr/bin/env python3
"""HW1 Problem 3: incident direction sweep."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from monte_carlo import MonteCarloTransport


def main():

    # ========================================================
    # Problem 3 parameters
    # ========================================================

    tau = 4.0
    l0 = 1.0

    # Number selected from Problem 1
    n_photons = 100_000

    # Incident direction cosines
    mu0_values = np.array([
        -0.1,
        -0.2,
        -0.3,
        -0.4,
        -0.5,
        -0.6,
        -0.7,
        -0.8,
        -0.9,
    ])

    # Storage
    R_values = []
    T_values = []
    sigma_R_values = []
    sigma_T_values = []

    print("Using up to 24 worker processes")

    print()
    print("HW1 - Problem 3")
    print("----------------")
    print(f"tau       = {tau}")
    print(f"N photons = {n_photons:,}")
    print()

    print(
        f"{'mu0':>8} "
        f"{'R':>12} "
        f"{'T':>12} "
        f"{'R+T':>12} "
        f"{'sigma_R':>12}"
    )

    print("-" * 62)

    # ========================================================
    # Run simulation for each mu0
    # ========================================================

    for i, mu0 in enumerate(mu0_values):

        # Different seed for each mu0 experiment
        seed = 1000 + i

        model = MonteCarloTransport(tau=tau, mu0=float(mu0), l0=l0)
        result = model.run(n_photons, seed=seed)
        R, T = result.reflectivity, result.transmissivity

        # Estimated Monte Carlo standard errors
        sigma_R = np.sqrt(R * (1.0 - R) / n_photons)
        sigma_T = np.sqrt(T * (1.0 - T) / n_photons)

        R_values.append(R)
        T_values.append(T)

        sigma_R_values.append(sigma_R)
        sigma_T_values.append(sigma_T)

        print(
            f"{mu0:8.1f} "
            f"{R:12.6f} "
            f"{T:12.6f} "
            f"{R+T:12.6f} "
            f"{sigma_R:12.6f}"
        )

    # Convert to NumPy arrays
    R_values = np.array(R_values)
    T_values = np.array(T_values)

    sigma_R_values = np.array(sigma_R_values)
    sigma_T_values = np.array(sigma_T_values)

    # ========================================================
    # Save numerical results
    # ========================================================

    output_dir = Path(__file__).resolve().parent

    output_table = np.column_stack(
        (
            mu0_values,
            R_values,
            T_values,
            sigma_R_values,
            sigma_T_values,
        )
    )

    np.savetxt(
        output_dir / "problem3_results.txt",
        output_table,
        header="mu0 R T sigma_R sigma_T",
        fmt="%.6f",
    )

    # ========================================================
    # Plot: R and T versus mu0
    # ========================================================

    plt.figure(figsize=(9, 6))

    plt.plot(
        mu0_values,
        R_values,
        marker="o",
        linewidth=2,
        label="Reflectivity (R)",
    )

    plt.plot(
        mu0_values,
        T_values,
        marker="s",
        linewidth=2,
        label="Transmissivity (T)",
    )

    plt.xlabel(
        r"Incident Direction Cosine, $\mu_0$"
    )

    plt.ylabel(
        "Reflectivity / Transmissivity"
    )

    plt.title(
        r"Effect of Incident Direction: "
        r"$\tau=4$, $N=100{,}000$"
    )

    plt.grid(True, alpha=0.3)
    plt.legend()

    plt.tight_layout()

    plt.savefig(
        output_dir / "problem3_mu0_reflectivity_transmissivity.png",
        dpi=300,
    )

    plt.show()


if __name__ == "__main__":
    main()
