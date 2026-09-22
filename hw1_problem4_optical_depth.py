#!/usr/bin/env python3
"""HW1 Problem 4: optical depth sweep."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from monte_carlo import MonteCarloTransport


def main():

    # --------------------------------------------------------
    # Problem 4 parameters
    # --------------------------------------------------------

    N = 100_000
    mu0 = -0.7
    l0 = 1.0

    optical_depths = np.arange(1, 9)

    seed = 42

    print("Using up to 24 worker processes")

    # --------------------------------------------------------
    # Storage
    # --------------------------------------------------------

    R_values = []
    T_values = []

    sigma_R_values = []

    computation_times = []

    mean_scattering_values = []

    # --------------------------------------------------------
    # Run simulations
    # --------------------------------------------------------

    for tau in optical_depths:

        print(f"Running tau = {tau} ...")

        model = MonteCarloTransport(tau=float(tau), mu0=mu0, l0=l0)
        result = model.run(N, seed=seed + int(tau) * 100)
        R, T = result.reflectivity, result.transmissivity
        mean_scatterings = result.mean_scatterings
        elapsed_time = result.elapsed_seconds

        sigma_R = np.sqrt(
            R * (1.0 - R) / N
        )

        R_values.append(R)
        T_values.append(T)

        sigma_R_values.append(sigma_R)

        computation_times.append(elapsed_time)

        mean_scattering_values.append(mean_scatterings)

    # Convert to arrays
    R_values = np.array(R_values)
    T_values = np.array(T_values)

    sigma_R_values = np.array(sigma_R_values)

    computation_times = np.array(computation_times)

    mean_scattering_values = np.array(mean_scattering_values)

    # ========================================================
    # Print results
    # ========================================================

    print()
    print("HW1 - Problem 4")
    print("----------------")
    print(f"N photons = {N:,}")
    print(f"mu0       = {mu0}")
    print()

    print(
        f"{'tau':>6} "
        f"{'R':>12} "
        f"{'T':>12} "
        f"{'R+T':>12} "
        f"{'sigma_R':>12} "
        f"{'Mean Scat.':>14} "
        f"{'Time (s)':>12}"
    )

    print("-" * 86)

    for i, tau in enumerate(optical_depths):

        print(
            f"{tau:6.1f} "
            f"{R_values[i]:12.6f} "
            f"{T_values[i]:12.6f} "
            f"{R_values[i] + T_values[i]:12.6f} "
            f"{sigma_R_values[i]:12.6f} "
            f"{mean_scattering_values[i]:14.3f} "
            f"{computation_times[i]:12.4f}"
        )

    # ========================================================
    # Save numerical results
    # ========================================================

    output_dir = Path(__file__).resolve().parent

    output_data = np.column_stack(
        (
            optical_depths,
            R_values,
            T_values,
            sigma_R_values,
            mean_scattering_values,
            computation_times,
        )
    )

    np.savetxt(
        output_dir / "problem4_results.txt",
        output_data,
        header=(
            "tau R T sigma_R "
            "mean_scatterings_per_photon computation_time_seconds"
        ),
        fmt="%.8f",
    )

    # ========================================================
    # Figure 1: Reflectivity and transmissivity
    # ========================================================

    plt.figure(figsize=(9, 6))

    plt.plot(
        optical_depths,
        R_values,
        marker="o",
        linewidth=2,
        label="Reflectivity (R)",
    )

    plt.plot(
        optical_depths,
        T_values,
        marker="s",
        linewidth=2,
        label="Transmissivity (T)",
    )

    plt.xlabel("Optical Depth, τ")
    plt.ylabel("Reflectivity / Transmissivity")

    plt.title(
        rf"Effect of Optical Depth: "
        rf"$\mu_0={mu0}$, $N={N:,}$"
    )

    plt.xticks(optical_depths)

    plt.grid(True, alpha=0.3)
    plt.legend()

    plt.tight_layout()

    plt.savefig(
        output_dir /
        "problem4_reflectivity_transmissivity_vs_tau.png",
        dpi=300,
    )

    plt.close()

    # ========================================================
    # Figure 2: Computation time
    # ========================================================

    plt.figure(figsize=(9, 6))

    plt.plot(
        optical_depths,
        computation_times,
        marker="o",
        linewidth=2,
    )

    plt.xlabel("Optical Depth, τ")
    plt.ylabel("Computation Time (s)")

    plt.title(
        rf"Monte Carlo Computation Time: "
        rf"$N={N:,}$"
    )

    plt.xticks(optical_depths)

    plt.grid(True, alpha=0.3)

    plt.tight_layout()

    plt.savefig(
        output_dir /
        "problem4_computation_time_vs_tau.png",
        dpi=300,
    )

    plt.close()

    # ========================================================
    # Figure 3: Mean number of scatterings
    # ========================================================

    plt.figure(figsize=(9, 6))

    plt.plot(
        optical_depths,
        mean_scattering_values,
        marker="o",
        linewidth=2,
    )

    plt.xlabel("Optical Depth, τ")
    plt.ylabel("Mean Number of Scatterings per Photon")

    plt.title(
        rf"Photon Scattering vs Optical Depth: "
        rf"$\mu_0={mu0}$"
    )

    plt.xticks(optical_depths)

    plt.grid(True, alpha=0.3)

    plt.tight_layout()

    plt.savefig(
        output_dir /
        "problem4_mean_scatterings_vs_tau.png",
        dpi=300,
    )

    plt.close()

    print()
    print("Results saved to:")
    print(output_dir / "problem4_results.txt")

    print()
    print("Figures saved:")
    print(
        output_dir /
        "problem4_reflectivity_transmissivity_vs_tau.png"
    )
    print(
        output_dir /
        "problem4_computation_time_vs_tau.png"
    )
    print(
        output_dir /
        "problem4_mean_scatterings_vs_tau.png"
    )


if __name__ == "__main__":
    main()
