"""Run the HW1-style experiments requested for HW2."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from models.monte_carlo import MonteCarloTransport


N_PHOTONS = 10_000
G = 0.75
ALBEDOS = np.array([1.0, 0.95, 0.85, 0.75])
SEEDS = np.array([101, 202, 303, 404, 505, 606, 707, 808, 909, 1010])
OUTPUT_DIR = Path(__file__).resolve().parent / "results"


def run_model(tau, mu0, omega0, seed):
    model = MonteCarloTransport(
        tau=float(tau),
        mu0=float(mu0),
        g=G,
        omega0=float(omega0),
        workers=24,
    )
    return model.run(N_PHOTONS, seed=int(seed))


def save_table(name, header, rows):
    OUTPUT_DIR.mkdir(exist_ok=True)
    np.savetxt(OUTPUT_DIR / name, np.asarray(rows), header=header, fmt="%.8f")


def independent_runs():
    rows = []
    for omega0 in ALBEDOS:
        for run, seed in enumerate(SEEDS, start=1):
            result = run_model(4.0, -0.7, omega0, seed)
            rows.append(
                [
                    omega0,
                    run,
                    seed,
                    result.reflectivity,
                    result.transmissivity,
                    result.absorptivity,
                    result.mean_scatterings,
                    result.elapsed_seconds,
                ]
            )
    rows = np.asarray(rows)
    save_table(
        "independent_runs.txt",
        "omega0 run seed R T A mean_scatterings elapsed_seconds",
        rows,
    )

    nonabsorbing = rows[rows[:, 0] == 1.0]
    runs = nonabsorbing[:, 1]
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(runs, nonabsorbing[:, 3], "o-", label="Reflectivity (R)")
    ax.plot(runs, nonabsorbing[:, 4], "s-", label="Transmissivity (T)")
    ax.axhline(nonabsorbing[:, 3].mean(), linestyle="--", color="C0")
    ax.axhline(nonabsorbing[:, 4].mean(), linestyle="--", color="C1")
    ax.set(
        xlabel="Independent Monte Carlo Run",
        ylabel="Fraction of Incident Photons",
        title=r"Independent HG Runs: $\tau=4$, $\mu_0=-0.7$, $\omega_0=1$",
    )
    ax.set_xticks(runs)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "omega1_independent_runs.png", dpi=300)
    plt.close(fig)

    summary = []
    print("Independent-run statistics (10 runs each)")
    print(f"{'omega0':>8} {'mean R':>10} {'std R':>10} {'mean T':>10} "
          f"{'std T':>10} {'mean A':>10} {'std A':>10}")
    for omega0 in ALBEDOS:
        subset = rows[rows[:, 0] == omega0]
        means = subset[:, 3:6].mean(axis=0)
        stds = subset[:, 3:6].std(axis=0, ddof=1)
        summary.append([omega0, *means, *stds])
        print(
            f"{omega0:8.2f} {means[0]:10.6f} {stds[0]:10.6f} "
            f"{means[1]:10.6f} {stds[1]:10.6f} "
            f"{means[2]:10.6f} {stds[2]:10.6f}"
        )
    save_table(
        "independent_run_summary.txt",
        "omega0 mean_R mean_T mean_A std_R std_T std_A",
        summary,
    )
    return np.asarray(summary)


def incident_angle_sweep():
    mu_values = -np.arange(1, 10) / 10.0
    rows = []
    for j, omega0 in enumerate(ALBEDOS):
        for i, mu0 in enumerate(mu_values):
            result = run_model(4.0, mu0, omega0, 2_000 + 100 * j + i)
            rows.append(
                [omega0, mu0, result.reflectivity, result.transmissivity,
                 result.absorptivity, result.mean_scatterings, result.elapsed_seconds]
            )
    rows = np.asarray(rows)
    save_table(
        "incident_angle_sweep.txt",
        "omega0 mu0 R T A mean_scatterings elapsed_seconds",
        rows,
    )
    plot_rta_comparison(
        rows,
        x_column=1,
        xlabel=r"Incident Direction Cosine, $\mu_0$",
        title=r"HG Transport versus Incident Direction ($\tau=4$)",
        filename="albedo_comparison_vs_mu0.png",
    )
    plot_nonabsorbing_rt(
        rows,
        xlabel=r"Incident Direction Cosine, $\mu_0$",
        title=r"HG Transport versus Incident Direction ($\tau=4$, $\omega_0=1$)",
        filename="omega1_reflectivity_transmissivity_vs_mu0.png",
    )
    return rows


def optical_depth_sweep():
    tau_values = np.arange(1, 9)
    rows = []
    for j, omega0 in enumerate(ALBEDOS):
        for tau in tau_values:
            result = run_model(tau, -0.7, omega0, 3_000 + 100 * j + tau)
            rows.append(
                [omega0, tau, result.reflectivity, result.transmissivity,
                 result.absorptivity, result.mean_scatterings, result.elapsed_seconds]
            )
    rows = np.asarray(rows)
    save_table(
        "optical_depth_sweep.txt",
        "omega0 tau R T A mean_scatterings elapsed_seconds",
        rows,
    )
    plot_rta_comparison(
        rows,
        x_column=1,
        xlabel=r"Optical Depth, $\tau$",
        title=r"HG Transport versus Optical Depth ($\mu_0=-0.7$)",
        filename="albedo_comparison_vs_tau.png",
    )
    plot_nonabsorbing_rt(
        rows,
        xlabel=r"Optical Depth, $\tau$",
        title=r"HG Transport versus Optical Depth ($\mu_0=-0.7$, $\omega_0=1$)",
        filename="omega1_reflectivity_transmissivity_vs_tau.png",
    )

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for omega0 in ALBEDOS:
        subset = rows[rows[:, 0] == omega0]
        label = rf"$\omega_0={omega0:.2f}$"
        axes[0].plot(subset[:, 1], subset[:, 5], "o-", label=label)
        axes[1].plot(subset[:, 1], subset[:, 6], "o-", label=label)
    axes[0].set(xlabel=r"Optical Depth, $\tau$", ylabel="Mean Scatterings per Photon")
    axes[1].set(xlabel=r"Optical Depth, $\tau$", ylabel="Computation Time (s)")
    for ax in axes:
        ax.grid(True, alpha=0.3)
        ax.legend()
    fig.suptitle(r"Computational Behavior versus Optical Depth ($g=0.75$)")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "performance_vs_tau.png", dpi=300)
    plt.close(fig)
    return rows


def plot_rta_comparison(rows, x_column, xlabel, title, filename):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True)
    active_axes = axes.flat[:3]
    labels = ["Reflectivity (R)", "Transmissivity (T)", "Absorptivity (A)"]
    for omega0 in ALBEDOS:
        subset = rows[rows[:, 0] == omega0]
        for ax, column in zip(active_axes, (2, 3, 4)):
            ax.plot(
                subset[:, x_column],
                subset[:, column],
                "o-",
                linewidth=2.5,
                markersize=7,
                label=rf"$\omega_0={omega0:.2f}$",
            )
    for ax, ylabel in zip(active_axes, labels):
        ax.set_xlabel(xlabel, fontsize=16)
        ax.set_ylabel(ylabel, fontsize=16)
        ax.tick_params(axis="both", labelsize=14)
        ax.grid(True, alpha=0.3)
    active_axes[-1].legend(fontsize=13)
    axes.flat[-1].axis("off")
    fig.suptitle(title, fontsize=21)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(OUTPUT_DIR / filename, dpi=600)
    plt.close(fig)


def plot_nonabsorbing_rt(rows, xlabel, title, filename):
    subset = rows[rows[:, 0] == 1.0]
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(subset[:, 1], subset[:, 2], "o-", label="Reflectivity (R)")
    ax.plot(subset[:, 1], subset[:, 3], "s-", label="Transmissivity (T)")
    ax.set(xlabel=xlabel, ylabel="Fraction of Incident Photons", title=title)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / filename, dpi=300)
    plt.close(fig)


def print_conservation(name, rows):
    error = np.max(np.abs(rows[:, 2] + rows[:, 3] + rows[:, 4] - 1.0))
    print(f"Maximum |R + T + A - 1| for {name}: {error:.3e}")


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    print(f"Running HW2 experiments with N={N_PHOTONS:,} and g={G}")
    summary = independent_runs()
    angle_rows = incident_angle_sweep()
    tau_rows = optical_depth_sweep()
    print()
    print_conservation("incident-angle sweep", angle_rows)
    print_conservation("optical-depth sweep", tau_rows)
    print(f"Numerical results and figures saved in: {OUTPUT_DIR}")

    # Suppress a linter warning while retaining the saved summary for callers.
    return summary, angle_rows, tau_rows


if __name__ == "__main__":
    main()
