"""ATMO 656A HW3 Question 1: Eddington and quadrature two-stream R/T."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from models.two_stream import METHODS, reflect_transmit


OPTICAL_DEPTHS = np.linspace(0.0, 10.0, 201)
MU_VALUES = (0.9, 0.5, 0.1)
ALBEDOS = (1.0, 0.9)
G = 0.75
OUTPUT_DIR = Path(__file__).resolve().parent / "results"


def run_sweep():
    """Evaluate every requested parameter combination."""
    results = {}
    warnings = []
    for omega0 in ALBEDOS:
        for mu0 in MU_VALUES:
            for method in METHODS:
                reflected, transmitted = reflect_transmit(
                    OPTICAL_DEPTHS,
                    omega0=omega0,
                    mu0=mu0,
                    g=G,
                    method=method,
                )
                absorbed = 1.0 - reflected - transmitted
                if not all(
                    np.all(np.isfinite(values))
                    for values in (reflected, transmitted, absorbed)
                ):
                    raise RuntimeError("two-stream sweep produced a nonfinite value")
                if omega0 == 1.0 and not np.allclose(
                    reflected + transmitted,
                    1.0,
                    atol=1.0e-12,
                ):
                    raise RuntimeError("conservative two-stream solution lost energy")

                minimum_r = float(np.min(reflected))
                minimum_t = float(np.min(transmitted))
                if minimum_r < -1.0e-8 or minimum_t < -1.0e-8:
                    warnings.append((omega0, mu0, method, minimum_r, minimum_t))
                results[(omega0, mu0, method)] = (
                    reflected,
                    transmitted,
                    absorbed,
                )
    return results, warnings


def save_csv(results):
    OUTPUT_DIR.mkdir(exist_ok=True)
    with (OUTPUT_DIR / "hw3_question1_results.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["omega0", "mu0", "g", "method", "tau", "R", "T", "A"])
        for omega0 in ALBEDOS:
            for mu0 in MU_VALUES:
                for method in METHODS:
                    reflected, transmitted, absorbed = results[(omega0, mu0, method)]
                    writer.writerows(
                        (
                            omega0,
                            mu0,
                            G,
                            method,
                            float(tau),
                            float(r_value),
                            float(t_value),
                            float(a_value),
                        )
                        for tau, r_value, t_value, a_value in zip(
                            OPTICAL_DEPTHS,
                            reflected,
                            transmitted,
                            absorbed,
                        )
                    )


def save_figures(results):
    OUTPUT_DIR.mkdir(exist_ok=True)
    colors = {0.9: "C0", 0.5: "C1", 0.1: "C2"}
    styles = {"eddington": "-", "quadrature": "--"}

    for omega0 in ALBEDOS:
        for quantity, index, ylabel in (
            ("reflectivity", 0, "Reflectivity, R"),
            ("transmissivity", 1, "Transmissivity, T"),
        ):
            fig, ax = plt.subplots(figsize=(9, 6))
            for mu0 in MU_VALUES:
                for method in METHODS:
                    values = results[(omega0, mu0, method)][index]
                    ax.plot(
                        OPTICAL_DEPTHS,
                        values,
                        color=colors[mu0],
                        linestyle=styles[method],
                        linewidth=2.2,
                        label=rf"$\mu_0={mu0}$, {method.capitalize()}",
                    )
            ax.axhline(0.0, color="black", linewidth=0.8, alpha=0.45)
            ax.set(
                xlabel=r"Optical Depth, $\tau$",
                ylabel=ylabel,
                title=(
                    rf"Meador-Weaver Two-Stream {quantity.capitalize()}: "
                    rf"$g={G}$, $\omega_0={omega0}$"
                ),
                xlim=(0.0, 10.0),
                ylim=(-0.04, 1.04),
            )
            ax.tick_params(labelsize=12)
            ax.xaxis.label.set_size(14)
            ax.yaxis.label.set_size(14)
            ax.title.set_size(16)
            ax.grid(True, alpha=0.25)
            ax.legend(fontsize=10, ncol=2)
            fig.tight_layout()
            fig.savefig(
                OUTPUT_DIR / f"hw3_q1_{quantity}_omega_{omega0:.1f}.png",
                dpi=400,
            )
            plt.close(fig)


def main():
    results, warnings = run_sweep()
    save_csv(results)
    save_figures(results)

    if warnings:
        print("Two-stream closure warnings (values retained, not clipped):")
        for omega0, mu0, method, minimum_r, minimum_t in warnings:
            print(
                f"  omega0={omega0}, mu0={mu0}, {method}: "
                f"min R={minimum_r:.8f}, min T={minimum_t:.8f}"
            )
    print(f"Saved four figures and 2,412 CSV rows in: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
