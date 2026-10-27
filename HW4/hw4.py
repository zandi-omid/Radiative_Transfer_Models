"""ATMO 656A HW4: discrete ordinates with all surface-cloud reflections."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from models.surface_two_stream import surface_transport


OUTPUT_DIR = Path(__file__).resolve().parent
MU0 = 0.9
G = 0.75
SURFACE_ALBEDOS = np.linspace(0.0, 1.0, 11)
CASES = ((1.0, 0.99999), (4.0, 0.99999), (1.0, 0.8), (4.0, 0.8))


def energy_budget_residual(result, incident_flux=1.0):
    """Validate normalized energy fractions without hiding numerical artifacts."""
    fractions = (
        result.toa_reflectivity,
        result.surface_net_flux / incident_flux,
        result.atmospheric_absorption / incident_flux,
    )
    if not all(np.all(np.isfinite(value)) for value in fractions):
        raise RuntimeError("HW4 energy fractions contain nonfinite values")
    residual = float(np.max(np.abs(sum(fractions) - 1.0)))
    if residual > 1.0e-11:
        raise RuntimeError(f"HW4 normalized energy residual is {residual:.3e}")
    return residual


def run_sweep():
    """Run four atmospheric cases with pi*F0=1 and report questionable values."""
    results = {}
    for tau, omega0 in CASES:
        result = surface_transport(
            tau, SURFACE_ALBEDOS, omega0, MU0, G, incident_flux=MU0
        )
        values = (
            result.toa_reflectivity,
            result.surface_downward_flux,
            result.surface_upward_flux,
            result.surface_net_flux,
            result.atmospheric_absorption,
        )
        if not all(np.all(np.isfinite(value)) for value in values):
            raise RuntimeError("HW4 sweep produced nonfinite output")
        balance = MU0 * result.toa_reflectivity
        balance += result.surface_net_flux + result.atmospheric_absorption
        if not np.allclose(balance, MU0, rtol=0.0, atol=1e-11):
            raise RuntimeError("HW4 sweep failed energy balance")
        if (
            np.any(result.toa_reflectivity < -1e-12)
            or np.any(result.toa_reflectivity > 1.0 + 1e-12)
            or np.any(result.surface_net_flux < -1e-12)
            or np.any(result.atmospheric_absorption < -1e-12)
        ):
            print(f"WARNING: unphysical fractions retained for tau={tau}, omega0={omega0}")
        results[(tau, omega0)] = result
    return results


def save_csv(results):
    """Export actual fluxes for pi*F0=1 (horizontal incident flux is 0.9)."""
    output = OUTPUT_DIR / "hw4_results.csv"
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow([
            "tau", "omega0", "g", "mu0", "surface_albedo", "toa_reflectivity",
            "surface_downward_flux", "surface_upward_flux", "surface_net_flux",
            "atmospheric_absorption", "normalized_surface_net_flux",
            "normalized_atmospheric_absorption",
        ])
        for (tau, omega0), result in results.items():
            for index, albedo in enumerate(SURFACE_ALBEDOS):
                writer.writerow([
                    tau, omega0, G, MU0, albedo, result.toa_reflectivity[index],
                    result.surface_downward_flux[index], result.surface_upward_flux[index],
                    result.surface_net_flux[index], result.atmospheric_absorption[index],
                    result.surface_net_flux[index] / MU0,
                    result.atmospheric_absorption[index] / MU0,
                ])
    return output


def save_figures():
    """Sample the same solution finely for smooth publication figures."""
    figure_dir = OUTPUT_DIR / "figures"
    figure_dir.mkdir(exist_ok=True)
    albedos = np.linspace(0.0, 1.0, 501)
    plot_results = {
        (tau, omega0): surface_transport(tau, albedos, omega0, MU0, G)
        for tau, omega0 in CASES
    }
    maximum_residual = max(
        energy_budget_residual(result) for result in plot_results.values()
    )
    for quantity, filename, ylabel in (
        ("toa_reflectivity", "reflectivity_vs_albedo", "Top-of-atmosphere reflectivity"),
        ("surface_net_flux", "surface_net_flux_vs_albedo", "Net surface flux / incident horizontal flux"),
    ):
        fig, ax = plt.subplots(figsize=(8, 5.5))
        for tau, omega0 in CASES:
            result = plot_results[(tau, omega0)]
            ax.plot(
                albedos, getattr(result, quantity),
                color="C0" if tau == 1.0 else "C1",
                linestyle="-" if omega0 == 0.99999 else "--",
                linewidth=2.3,
                label=rf"$\tau={tau:g},\ \omega_0={omega0:g}$",
            )
        ax.set(xlabel=r"Surface albedo, $A_s$", ylabel=ylabel, xlim=(0, 1))
        ax.set_title(r"Liou two-stream model: $\mu_0=0.9$, $g=0.75$", fontsize=15)
        ax.tick_params(labelsize=12)
        ax.xaxis.label.set_size(14)
        ax.yaxis.label.set_size(13)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=11)
        fig.tight_layout()
        for extension in ("png", "pdf"):
            fig.savefig(figure_dir / f"{filename}.{extension}", dpi=300)
        plt.close(fig)

    save_energy_partition(plot_results, albedos)
    return maximum_residual


@matplotlib.rc_context({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Liberation Serif", "STIXGeneral"],
    "mathtext.fontset": "stix",
    "pdf.fonttype": 42,
})
def save_energy_partition(plot_results, albedos):
    """Plot normalized fractions; infer panel (c) from the energy budget.

    Inputs have unit horizontal incident flux. The residual closure is an
    arithmetic identity, not independent validation of the transport solver.
    """
    figure_dir = OUTPUT_DIR / "figures"
    figure_dir.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.4), sharex=True, sharey=True)
    panels = (
        ("toa_reflectivity", "(a) Top-of-atmosphere reflectivity", r"$R_{\mathrm{TOA}}$"),
        ("surface_net_flux", "(b) Normalized net surface flux",
         r"$F_{\mathrm{net},s}/F_{\mathrm{incident}}$"),
        ("atmospheric_absorption", "(c) Atmospheric absorption",
         r"$A_{\mathrm{atmosphere}}$"),
    )
    for ax, (quantity, title, ylabel) in zip(axes, panels):
        for tau, omega0 in CASES:
            result = plot_results[(tau, omega0)]
            values = getattr(result, quantity)
            if quantity == "atmospheric_absorption":
                values = 1.0 - result.toa_reflectivity - result.surface_net_flux
            ax.plot(
                albedos, values,
                color="C0" if tau == 1.0 else "C1",
                linestyle="-" if omega0 == 0.99999 else "--",
                linewidth=2.4,
                label=rf"$\tau={tau:g},\ \omega_0={omega0:g}$",
            )
        ax.set(xlabel=r"Surface albedo, $A_s$", ylabel=ylabel,
               xlim=(0.0, 1.0), ylim=(-0.015, 1.015))
        ax.set_title(title, fontsize=14, pad=12)
        ax.xaxis.label.set_size(14)
        ax.yaxis.label.set_size(15)
        ax.tick_params(labelsize=12, labelleft=True)
        ax.set_xticks(np.linspace(0.0, 1.0, 6))
        ax.set_yticks(np.linspace(0.0, 1.0, 6))
        ax.grid(alpha=0.25)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.015),
               ncol=4, frameon=False, fontsize=12)
    fig.subplots_adjust(left=0.065, right=0.985, bottom=0.23, top=0.87, wspace=0.29)
    for extension in ("png", "pdf"):
        fig.savefig(figure_dir / f"energy_partition_three_panels.{extension}", dpi=300)
    plt.close(fig)


def main():
    results = run_sweep()
    csv_path = save_csv(results)
    plot_residual = save_figures()
    csv_residual = max(
        energy_budget_residual(result, incident_flux=MU0)
        for result in results.values()
    )
    print(
        "Maximum normalized energy-budget residual "
        f"(CSV and plotted cases): {max(csv_residual, plot_residual):.3e}"
    )
    print("pi*F0=1; incident horizontal flux=0.9; flux columns below use these units")
    print(" tau  omega0   As       R_TOA     F_net      A_atm")
    for (tau, omega0), result in results.items():
        for index in (0, 5, 10):
            print(
                f" {tau:3g}  {omega0:7g}  {SURFACE_ALBEDOS[index]:.1f}"
                f"  {result.toa_reflectivity[index]:.6f}"
                f"  {result.surface_net_flux[index]:.6f}"
                f"  {result.atmospheric_absorption[index]:.6f}"
            )
        print(f"  minimum reflection denominator: {np.min(result.reflection_denominator):.6f}")
    print(f"Saved {csv_path} and PNG/PDF figures in {OUTPUT_DIR / 'figures'}")


if __name__ == "__main__":
    main()
