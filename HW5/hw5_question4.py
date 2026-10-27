"""HW5 Question 4: angular Mie phase functions for 4-micrometer spheres."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from models.mie import mie_phase_function


OUTPUT_DIR = Path(__file__).resolve().parent
RADIUS_UM = 4.0
WAVELENGTH_UM = 0.55
CASES = ((1.5, 0.0, "C0"), (2.0, 0.66, "C3"))
ANGLES_DEG = np.arange(181, dtype=float)


def phase_moments(real_index, absorption_index):
    """Integrate the Fortran angular distribution with 128-point quadrature."""
    mu, weights = np.polynomial.legendre.leggauss(128)
    values = mie_phase_function(
        RADIUS_UM, real_index, absorption_index, WAVELENGTH_UM,
        angles_degrees=np.rad2deg(np.arccos(mu)),
    )
    normalization = 0.5 * np.dot(weights, values["phase_4pi"])
    angular_g = 0.5 * np.dot(weights * mu, values["phase_4pi"])
    return float(normalization), float(angular_g), float(values["g"])


def run_sweep():
    results = {}
    for n, k, _ in CASES:
        values = mie_phase_function(RADIUS_UM, n, k, WAVELENGTH_UM, ANGLES_DEG)
        if np.any(values["phase_4pi"] <= 0):
            raise RuntimeError("Nonpositive phase values cannot be plotted on a log axis")
        norm, angular_g, series_g = phase_moments(n, k)
        if abs(norm - 1) > 5e-10 or abs(angular_g - series_g) > 5e-10:
            raise RuntimeError("Phase normalization or angular first moment check failed")
        theta = np.deg2rad(values["angle_deg"])
        # Audit the requested grid without renormalizing away discretization error.
        trapezoid = getattr(np, "trapezoid", None) or np.trapz
        coarse_norm = 0.5 * trapezoid(values["phase_4pi"] * np.sin(theta), theta)
        values.update(normalization=norm, angular_g=angular_g,
                      one_degree_normalization=float(coarse_norm))
        results[(n, k)] = values
    return results


def save_results(results):
    directory = OUTPUT_DIR / "results"
    directory.mkdir(exist_ok=True)
    output = directory / "question4_phase_function.csv"
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["real_index", "absorption_index", "radius_um", "wavelength_um",
                         "x", "angle_deg", "phase_4pi", "phase_per_sr"])
        for (n, k), values in results.items():
            writer.writerows(
                (n, k, RADIUS_UM, WAVELENGTH_UM, values["x"], *row)
                for row in zip(values["angle_deg"], values["phase_4pi"], values["phase_per_sr"])
            )
    audit = directory / "question4_moments.csv"
    with audit.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["real_index", "absorption_index", "qext", "qscat", "g_series",
                         "g_angular", "phase_normalization", "one_degree_normalization"])
        for (n, k), values in results.items():
            writer.writerow([n, k, values["qext"], values["qscat"], values["g"],
                             values["angular_g"], values["normalization"],
                             values["one_degree_normalization"]])
    return output


def save_figure(results):
    directory = OUTPUT_DIR / "figures"
    directory.mkdir(exist_ok=True)
    with matplotlib.rc_context({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Liberation Serif", "STIXGeneral"],
        "mathtext.fontset": "stix", "pdf.fonttype": 42,
    }):
        fig, ax = plt.subplots(figsize=(9, 5.6))
        for n, k, color in CASES:
            values = results[(n, k)]
            ax.semilogy(values["angle_deg"], values["phase_4pi"],
                        color=color, linewidth=1.8, label=rf"$m={n:g}-{k:g}i$")
        ax.set(xlim=(0, 180), xticks=np.arange(0, 181, 30),
               xlabel=r"Scattering angle, $\Theta$ (degrees)",
               ylabel=r"Phase function, $P(\Theta)$")
        ax.set_title(r"Mie phase functions: $r=4\,\mu\mathrm{m}$, "
                     r"$\lambda=0.55\,\mu\mathrm{m}$", fontsize=16, pad=12)
        ax.tick_params(labelsize=12)
        ax.xaxis.label.set_size(14)
        ax.yaxis.label.set_size(14)
        ax.grid(which="major", alpha=0.25)
        ax.legend(fontsize=12)
        fig.tight_layout()
        for extension in ("png", "pdf"):
            fig.savefig(directory / f"question4_phase_function.{extension}", dpi=300)
        plt.close(fig)


def main():
    results = run_sweep()
    output = save_results(results)
    save_figure(results)
    for (n, k), values in results.items():
        print(f"m={n}-{k}i: x={values['x']:.8f}, g={values['g']:.9f}")
        for angle in (0, 90, 180):
            print(f"  P({angle:3d} deg)={values['phase_4pi'][angle]:.9g}")
        print(f"  Quadrature normalization={values['normalization']:.12f}; "
              f"angular g={values['angular_g']:.12f}")
        print(f"  1-degree trapezoidal normalization={values['one_degree_normalization']:.9f}")
    print(f"Saved 362 angular samples: {output}")
    print(f"Saved PNG/PDF figure: {OUTPUT_DIR / 'figures'}")


if __name__ == "__main__":
    main()
