#!/usr/bin/env python3
"""HW1 Problem 1 Monte Carlo analysis."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from models.monte_carlo import MonteCarloTransport


def main():
    # ============================================================
    # Problem 1 parameters
    # ============================================================
    tau = 4.0
    mu0 = -0.7
    l0 = 1.0

    photon_numbers = np.array([
        100,
        500,
        1_000,
        5_000,
        10_000,
        50_000,
        100_000,
        500_000,
        1_000_000,
    ])

    # ============================================================
    # Multiprocessing settings
    # ============================================================
    seed = 42
    model = MonteCarloTransport(tau=tau, mu0=mu0, l0=l0)
    n_workers = model.workers

    print(f"Using {n_workers} worker processes")

    # ============================================================
    # Storage
    # ============================================================
    R_values = []
    T_values = []

    sigma_R_values = []
    sigma_T_values = []

    relative_change_R = []
    relative_change_T = []

    # ============================================================
    # Run simulations
    # ============================================================
    previous_R = None
    previous_T = None

    for N in photon_numbers:
        result = model.run(int(N), seed=seed + int(N))
        R, T = result.reflectivity, result.transmissivity

        # Binomial Monte Carlo standard errors
        sigma_R = np.sqrt(R * (1.0 - R) / N)
        sigma_T = np.sqrt(T * (1.0 - T) / N)

        # Relative change from previous photon number
        if previous_R is None:
            change_R = np.nan
            change_T = np.nan
        else:
            change_R = abs(R - previous_R) / previous_R * 100.0
            change_T = abs(T - previous_T) / previous_T * 100.0

        R_values.append(R)
        T_values.append(T)

        sigma_R_values.append(sigma_R)
        sigma_T_values.append(sigma_T)

        relative_change_R.append(change_R)
        relative_change_T.append(change_T)

        previous_R = R
        previous_T = T

    # Convert to NumPy arrays
    R_values = np.array(R_values)
    T_values = np.array(T_values)

    sigma_R_values = np.array(sigma_R_values)
    sigma_T_values = np.array(sigma_T_values)

    relative_change_R = np.array(relative_change_R)
    relative_change_T = np.array(relative_change_T)

    # ============================================================
    # Print results
    # ============================================================
    print()
    print("Problem 1: Monte Carlo Convergence")
    print(f"tau = {tau}, mu0 = {mu0}")
    print()

    print(
        f"{'N':>10} "
        f"{'R':>10} "
        f"{'T':>10} "
        f"{'sigma':>12} "
        f"{'dR (%)':>12} "
        f"{'dT (%)':>12}"
    )

    print("-" * 72)

    for i, N in enumerate(photon_numbers):
        if i == 0:
            print(
                f"{N:10d} "
                f"{R_values[i]:10.5f} "
                f"{T_values[i]:10.5f} "
                f"{sigma_R_values[i]:12.6f} "
                f"{'---':>12} "
                f"{'---':>12}"
            )
        else:
            print(
                f"{N:10d} "
                f"{R_values[i]:10.5f} "
                f"{T_values[i]:10.5f} "
                f"{sigma_R_values[i]:12.6f} "
                f"{relative_change_R[i]:12.3f} "
                f"{relative_change_T[i]:12.3f}"
            )

    # ============================================================
    # Find first point where both R and T change by < 1%
    # ============================================================
    threshold = 1.0
    converged_N = None

    for i in range(1, len(photon_numbers)):
        if (
            relative_change_R[i] < threshold
            and relative_change_T[i] < threshold
        ):
            converged_N = photon_numbers[i]
            break

    print()
    if converged_N is not None:
        print(
            f"First photon number where both R and T change "
            f"by less than {threshold:.1f}%: N = {converged_N:,}"
        )
    else:
        print(
            f"No photon number satisfied the "
            f"{threshold:.1f}% criterion."
        )

    output_dir = Path(__file__).resolve().parent

    # ============================================================
    # Plot 1: R and T convergence
    # ============================================================
    plt.figure(figsize=(9, 6))

    plt.plot(
        photon_numbers,
        R_values,
        marker="o",
        linewidth=2,
        label="Reflectivity (R)",
    )

    plt.plot(
        photon_numbers,
        T_values,
        marker="s",
        linewidth=2,
        label="Transmissivity (T)",
    )

    plt.xscale("log")

    plt.xlabel("Number of Photons")
    plt.ylabel("Reflectivity / Transmissivity")

    plt.title(
        r"Monte Carlo Convergence: "
        r"$\tau=4$, $\mu_0=-0.7$"
    )

    plt.grid(True, alpha=0.3)
    plt.legend()

    plt.tight_layout()
    plt.savefig(output_dir / "reflectivity_transmissivity_convergence.png", dpi=300)
    plt.show()

    # ============================================================
    # Plot 2: Standard error
    # ============================================================
    plt.figure(figsize=(9, 6))

    plt.plot(
        photon_numbers,
        sigma_R_values,
        marker="o",
        linewidth=2,
    )

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("Number of Photons")
    plt.ylabel("Standard Error")

    plt.title(
        r"Monte Carlo Standard Error: "
        r"$\sigma_R=\sqrt{R(1-R)/N}$"
    )

    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / "standard_error_convergence.png", dpi=300)
    plt.show()

    # ============================================================
    # Plot 3: Relative change between successive simulations
    # ============================================================
    plt.figure(figsize=(9, 6))

    plt.plot(
        photon_numbers[1:],
        relative_change_R[1:],
        marker="o",
        linewidth=2,
        label="Relative change in R",
    )

    plt.plot(
        photon_numbers[1:],
        relative_change_T[1:],
        marker="s",
        linewidth=2,
        label="Relative change in T",
    )

    # 1% convergence reference
    plt.axhline(
        y=1.0,
        linestyle="--",
        label="1% criterion",
    )

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("Number of Photons")
    plt.ylabel("Relative Change (%)")

    plt.title("Change Between Successive Monte Carlo Estimates")

    plt.grid(True, alpha=0.3)
    plt.legend()

    plt.tight_layout()
    plt.savefig(output_dir / "relative_change_convergence.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    main()
