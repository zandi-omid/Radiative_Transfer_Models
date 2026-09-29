"""Compare the computational cost of isotropic and HG scattering."""

from multiprocessing import cpu_count
from pathlib import Path
from time import perf_counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from models.monte_carlo import MonteCarloTransport


N_PHOTONS = 10_000
MU0 = -0.7
OMEGA0 = 1.0
OPTICAL_DEPTHS = np.arange(1, 9)
N_REPEATS = 5
N_WORKERS = min(cpu_count(), 24)
MODELS = ((0.0, "Isotropic ($g=0$)"), (0.75, "HG ($g=0.75$)"))
OUTPUT_DIR = Path(__file__).resolve().parent / "results"


def run_benchmark():
    """Return one row per timed repetition.

    The same seed is used for both scattering models in each tau/repeat pair.
    Their execution order alternates to reduce systematic first-run bias.
    """
    rows = []
    for tau in OPTICAL_DEPTHS:
        models = {
            g: MonteCarloTransport(
                tau=float(tau),
                mu0=MU0,
                g=g,
                omega0=OMEGA0,
                workers=N_WORKERS,
            )
            for g, _ in MODELS
        }
        for repeat in range(1, N_REPEATS + 1):
            seed = 5_000 + 100 * int(tau) + repeat
            order = (0.0, 0.75) if repeat % 2 else (0.75, 0.0)
            for g in order:
                start = perf_counter()
                result = models[g].run(N_PHOTONS, seed=seed)
                wall_time = perf_counter() - start
                rows.append(
                    [
                        tau,
                        g,
                        repeat,
                        seed,
                        wall_time,
                        result.mean_scatterings,
                        result.reflectivity,
                        result.transmissivity,
                    ]
                )
    return np.asarray(rows)


def summarize(raw_rows):
    rows = []
    for tau in OPTICAL_DEPTHS:
        for g, _ in MODELS:
            subset = raw_rows[(raw_rows[:, 0] == tau) & (raw_rows[:, 1] == g)]
            rows.append(
                [
                    tau,
                    g,
                    N_WORKERS,
                    subset[:, 4].mean(),
                    subset[:, 4].std(ddof=1),
                    subset[:, 5].mean(),
                    subset[:, 5].std(ddof=1),
                ]
            )
    return np.asarray(rows)


def save_results(raw_rows, summary):
    OUTPUT_DIR.mkdir(exist_ok=True)
    np.savetxt(
        OUTPUT_DIR / "isotropic_hg_timing_raw.txt",
        raw_rows,
        header="tau g repeat seed wall_time_seconds mean_scatterings R T",
        fmt=["%.0f", "%.2f", "%.0f", "%.0f", "%.8f", "%.8f", "%.8f", "%.8f"],
    )
    np.savetxt(
        OUTPUT_DIR / "isotropic_hg_timing_summary.txt",
        summary,
        header=(
            "tau g workers mean_time_seconds std_time_seconds "
            "mean_scatterings std_scatterings"
        ),
        fmt=["%.0f", "%.2f", "%.0f", "%.8f", "%.8f", "%.8f", "%.8f"],
    )


def make_figure(summary):
    fig, ax = plt.subplots(figsize=(9, 6))
    for g, label in MODELS:
        subset = summary[summary[:, 1] == g]
        ax.errorbar(
            subset[:, 0],
            subset[:, 3],
            yerr=subset[:, 4],
            marker="o",
            linewidth=2,
            capsize=4,
            label=label,
        )

    ax.set(
        xlabel=r"Optical Depth, $\tau$",
        ylabel="Mean Wall-Clock Time (s)",
        title=f"Computation Time ({N_REPEATS} Repeats)",
    )
    ax.set_xticks(OPTICAL_DEPTHS)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.suptitle(
        rf"Isotropic and HG Cost Comparison: $N={N_PHOTONS:,}$, "
        rf"$\mu_0={MU0}$, $\omega_0={OMEGA0}$, {N_WORKERS} workers"
    )
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "isotropic_hg_computational_comparison.png", dpi=300)
    plt.close(fig)


def print_summary(summary):
    print(
        f"Computational comparison: N={N_PHOTONS:,}, mu0={MU0}, "
        f"omega0={OMEGA0}, workers={N_WORKERS}, repeats={N_REPEATS}"
    )
    print(
        f"{'tau':>4} {'model':>10} {'mean scat.':>12} "
        f"{'mean time (s)':>14} {'std time (s)':>13}"
    )
    print("-" * 60)
    names = {0.0: "isotropic", 0.75: "HG"}
    for tau in OPTICAL_DEPTHS:
        for g, _ in MODELS:
            row = summary[(summary[:, 0] == tau) & (summary[:, 1] == g)][0]
            print(
                f"{int(tau):4d} {names[g]:>10} {row[5]:12.3f} "
                f"{row[3]:14.4f} {row[4]:13.4f}"
            )


def main():
    raw_rows = run_benchmark()
    summary = summarize(raw_rows)
    save_results(raw_rows, summary)
    make_figure(summary)
    print_summary(summary)
    print(f"\nResults saved in: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
