"""HW2 Part 1: repeat the HW1 convergence study with HG scattering."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from models.monte_carlo import MonteCarloTransport


def save_line_plot(x, series, ylabel, title, output_path, log_y=False):
    """Save one consistently formatted convergence figure."""
    fig, ax = plt.subplots(figsize=(9, 6))
    for values, marker, label in series:
        ax.plot(x, values, marker=marker, linewidth=2, label=label)
    ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
    ax.set_xlabel("Number of Photons")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    if any(label for _, _, label in series):
        ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main():
    photon_numbers = np.array([100, 500, 1_000, 5_000, 10_000])
    tau = 4.0
    mu0 = -0.7
    g = 0.75
    omega0 = 1.0
    seed = 42

    model = MonteCarloTransport(tau=tau, mu0=mu0, g=g, workers=24)
    results = [model.run(int(n), seed=seed + int(n)) for n in photon_numbers]

    reflectivity = np.array([result.reflectivity for result in results])
    transmissivity = np.array([result.transmissivity for result in results])
    sigma_r = np.sqrt(reflectivity * (1.0 - reflectivity) / photon_numbers)
    relative_r = np.full(photon_numbers.size, np.nan)
    relative_t = np.full(photon_numbers.size, np.nan)
    relative_r[1:] = np.abs(np.diff(reflectivity) / reflectivity[:-1]) * 100.0
    relative_t[1:] = np.abs(np.diff(transmissivity) / transmissivity[:-1]) * 100.0

    print("HW2 - Part 1: HG Convergence Study")
    print("----------------------------------")
    print(f"tau={tau}, mu0={mu0}, g={g}, omega0={omega0}")
    print()
    print(
        f"{'N':>10} {'R':>10} {'T':>10} {'R+T':>10} "
        f"{'sigma_R':>12} {'dR (%)':>10} {'dT (%)':>10}"
    )
    print("-" * 80)
    for i, n in enumerate(photon_numbers):
        dr = "---" if i == 0 else f"{relative_r[i]:.3f}"
        dt = "---" if i == 0 else f"{relative_t[i]:.3f}"
        print(
            f"{n:10d} {reflectivity[i]:10.5f} {transmissivity[i]:10.5f} "
            f"{reflectivity[i] + transmissivity[i]:10.5f} "
            f"{sigma_r[i]:12.6f} {dr:>10} {dt:>10}"
        )

    final = results[-1]
    print()
    print("Final N=10,000 result")
    print(f"Reflectivity        : {final.reflectivity:.6f}")
    print(f"Transmissivity      : {final.transmissivity:.6f}")
    print(f"R + T               : {final.reflectivity + final.transmissivity:.6f}")
    print(f"Mean scatterings    : {final.mean_scatterings:.3f}")
    print(f"Elapsed time        : {final.elapsed_seconds:.3f} s")

    output_dir = Path(__file__).resolve().parent
    np.savetxt(
        output_dir / "part1_hg_convergence_results.txt",
        np.column_stack(
            (photon_numbers, reflectivity, transmissivity, sigma_r, relative_r, relative_t)
        ),
        header="N R T sigma_R relative_change_R_percent relative_change_T_percent",
        fmt=["%d", "%.8f", "%.8f", "%.8f", "%.8f", "%.8f"],
    )

    parameters = r"$\tau=4$, $\mu_0=-0.7$, $g=0.75$, $\omega_0=1$"
    save_line_plot(
        photon_numbers,
        [
            (reflectivity, "o", "Reflectivity (R)"),
            (transmissivity, "s", "Transmissivity (T)"),
        ],
        "Reflectivity / Transmissivity",
        f"HG Monte Carlo Convergence: {parameters}",
        output_dir / "part1_hg_reflectivity_transmissivity.png",
    )
    save_line_plot(
        photon_numbers,
        [(sigma_r, "o", None)],
        "Standard Error of Reflectivity",
        f"HG Monte Carlo Standard Error: {parameters}",
        output_dir / "part1_hg_standard_error.png",
        log_y=True,
    )

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(photon_numbers[1:], relative_r[1:], marker="o", linewidth=2, label="R")
    ax.plot(photon_numbers[1:], relative_t[1:], marker="s", linewidth=2, label="T")
    ax.axhline(1.0, linestyle="--", color="black", label="1% criterion")
    ax.set(xlabel="Number of Photons", ylabel="Relative Change (%)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(f"Change Between Successive HG Estimates: {parameters}")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "part1_hg_relative_change.png", dpi=300)
    plt.close(fig)

    print()
    print(f"Results and three figures saved in: {output_dir}")


if __name__ == "__main__":
    main()
