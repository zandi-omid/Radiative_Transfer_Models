"""HW3 Question 2: compare Liou and Meador-Weaver quadrature solutions."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from models.liou_two_stream import liou_reflect_transmit
from models.two_stream import reflect_transmit


OPTICAL_DEPTHS = np.linspace(0.0, 10.0, 201)
MU_VALUES = (0.9, 0.5, 0.1)
ALBEDOS = (1.0, 0.9)
G = 0.75
OUTPUT_DIR = Path(__file__).resolve().parent / "results"


def run_sweep():
    """Evaluate both formulations and collect unphysical-value warnings."""
    results = {}
    warnings = []
    for omega0 in ALBEDOS:
        for mu0 in MU_VALUES:
            liou = liou_reflect_transmit(OPTICAL_DEPTHS, omega0, mu0, G)
            quadrature = reflect_transmit(
                OPTICAL_DEPTHS, omega0, mu0, G, method="quadrature"
            )
            for method, values in (("liou", liou), ("quadrature", quadrature)):
                reflected, transmitted = values
                if not np.all(np.isfinite(reflected)) or not np.all(
                    np.isfinite(transmitted)
                ):
                    raise RuntimeError(f"{method} produced a nonfinite value")
                tolerance = 1.0e-12
                questionable = (reflected < -tolerance) | (
                    reflected > 1.0 + tolerance
                )
                questionable |= (transmitted < -tolerance) | (
                    transmitted > 1.0 + tolerance
                )
                if np.any(questionable):
                    warnings.append(
                        (
                            omega0,
                            mu0,
                            method,
                            float(np.min(reflected)),
                            float(np.max(reflected)),
                            float(np.min(transmitted)),
                            float(np.max(transmitted)),
                        )
                    )
                results[(omega0, mu0, method)] = values
    return results, warnings


def save_csv(results):
    """Write both solutions and their pointwise differences."""
    OUTPUT_DIR.mkdir(exist_ok=True)
    output = OUTPUT_DIR / "hw3_question2_results.csv"
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "omega0",
                "mu0",
                "g",
                "tau",
                "liou_R",
                "liou_T",
                "quadrature_R",
                "quadrature_T",
                "delta_R",
                "delta_T",
            ]
        )
        for omega0 in ALBEDOS:
            for mu0 in MU_VALUES:
                liou_r, liou_t = results[(omega0, mu0, "liou")]
                quad_r, quad_t = results[(omega0, mu0, "quadrature")]
                writer.writerows(
                    (
                        omega0,
                        mu0,
                        G,
                        float(depth),
                        float(lr),
                        float(lt),
                        float(qr),
                        float(qt),
                        float(lr - qr),
                        float(lt - qt),
                    )
                    for depth, lr, lt, qr, qt in zip(
                        OPTICAL_DEPTHS, liou_r, liou_t, quad_r, quad_t
                    )
                )
    return output


def save_figures(results):
    """Save R and T comparisons for each single-scattering albedo."""
    OUTPUT_DIR.mkdir(exist_ok=True)
    colors = {0.9: "C0", 0.5: "C1", 0.1: "C2"}
    styles = {"quadrature": "-", "liou": "--"}
    outputs = []
    for omega0 in ALBEDOS:
        for quantity, index, ylabel in (
            ("reflectivity", 0, "Reflectivity, R"),
            ("transmissivity", 1, "Transmissivity, T"),
        ):
            fig, ax = plt.subplots(figsize=(9, 6))
            for mu0 in MU_VALUES:
                for method in ("quadrature", "liou"):
                    ax.plot(
                        OPTICAL_DEPTHS,
                        results[(omega0, mu0, method)][index],
                        color=colors[mu0],
                        linestyle=styles[method],
                        linewidth=2.8 if method == "quadrature" else 2.0,
                        alpha=0.4 if method == "quadrature" else 1.0,
                        label=rf"$\mu_0={mu0}$, {method.capitalize()}",
                    )
            ax.axhline(0.0, color="black", linewidth=0.8, alpha=0.45)
            ax.set(
                xlabel=r"Optical Depth, $\tau$",
                ylabel=ylabel,
                title=(
                    rf"Liou vs. Meador-Weaver Quadrature: {ylabel}, "
                    rf"$g={G}$, $\omega_0={omega0}$"
                ),
                xlim=(0.0, 10.0),
                ylim=(-0.04, 1.04),
            )
            ax.tick_params(labelsize=12)
            ax.xaxis.label.set_size(14)
            ax.yaxis.label.set_size(14)
            ax.title.set_size(15)
            ax.grid(True, alpha=0.25)
            ax.legend(fontsize=10, ncol=2)
            fig.tight_layout()
            output = OUTPUT_DIR / f"hw3_q2_{quantity}_omega_{omega0:.1f}.png"
            fig.savefig(output, dpi=400)
            plt.close(fig)
            outputs.append(output)
    return outputs


def main():
    results, warnings = run_sweep()
    csv_path = save_csv(results)
    figures = save_figures(results)
    maximum_r_difference = max(
        float(
            np.max(
                np.abs(
                    results[(omega0, mu0, "liou")][0]
                    - results[(omega0, mu0, "quadrature")][0]
                )
            )
        )
        for omega0 in ALBEDOS
        for mu0 in MU_VALUES
    )
    maximum_t_difference = max(
        float(
            np.max(
                np.abs(
                    results[(omega0, mu0, "liou")][1]
                    - results[(omega0, mu0, "quadrature")][1]
                )
            )
        )
        for omega0 in ALBEDOS
        for mu0 in MU_VALUES
    )

    if warnings:
        print("Physically questionable values (retained, not clipped):")
        for omega0, mu0, method, min_r, max_r, min_t, max_t in warnings:
            print(
                f"  omega0={omega0}, mu0={mu0}, {method}: "
                f"R=[{min_r:.8f}, {max_r:.8f}], T=[{min_t:.8f}, {max_t:.8f}]"
            )
    print(f"Maximum |Liou R - quadrature R|: {maximum_r_difference:.3e}")
    print(f"Maximum |Liou T - quadrature T|: {maximum_t_difference:.3e}")
    print(f"Saved CSV: {csv_path}")
    print(f"Saved {len(figures)} figures in: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
