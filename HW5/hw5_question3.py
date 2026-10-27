"""HW5 Question 3: mass-specific extinction/absorption and g versus radius."""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from models.mie import COLUMNS, mass_specific_cross_section, mie_sweep


OUTPUT_DIR = Path(__file__).resolve().parent
WAVELENGTH_UM = 0.55
DENSITY_G_CM3 = 1.0
CASES = ((2.0, 0.66, "C3"), (1.5, 0.0, "C0"))
# Log spacing resolves the Rayleigh regime; a fine linear tail resolves
# interference features without smoothing in the large-size regime.
RADII_UM = np.unique(np.concatenate([
    np.geomspace(0.001, 0.1, 501), np.linspace(0.1, 10.0, 4951),
]))
SPECIFIC_COLUMNS = (
    "specific_extinction_m2_g", "specific_scattering_m2_g",
    "specific_absorption_m2_g",
)


def run_sweep(wavelength_um=WAVELENGTH_UM, density_g_cm3=DENSITY_G_CM3):
    """Run both required materials with identical core and outer radii."""
    if not np.isfinite(wavelength_um) or wavelength_um <= 0:
        raise ValueError("wavelength_um must be positive and finite")
    # Validate density before invoking the Fortran calculation.
    mass_specific_cross_section(1.0, 1.0, density_g_cm3)
    results = {}
    x = 2 * np.pi * RADII_UM / wavelength_um
    for real_index, absorption_index, _ in CASES:
        values = mie_sweep(x, real_index, absorption_index, wavelength_um)
        for quantity, column in zip(("qext", "qscat", "qabs"), SPECIFIC_COLUMNS):
            values[column] = mass_specific_cross_section(
                values[quantity], values["radius_um"], density_g_cm3
            )
        np.testing.assert_allclose(
            values[SPECIFIC_COLUMNS[0]],
            values[SPECIFIC_COLUMNS[1]] + values[SPECIFIC_COLUMNS[2]],
            rtol=1e-13, atol=1e-12,
        )
        if absorption_index == 0.0:
            residual = float(np.max(np.abs(values["qabs"])))
            print(f"Nonabsorbing control: max |Qabs|={residual:.3e}")
            print(f"  Signed specific absorption range: "
                  f"[{np.min(values[SPECIFIC_COLUMNS[2]]):.3e}, "
                  f"{np.max(values[SPECIFIC_COLUMNS[2]]):.3e}] m^2/g")
            if residual > 1e-10:
                raise RuntimeError("Nonabsorbing sphere has significant absorption")
        results[(real_index, absorption_index)] = values
    return results


def save_results(results, wavelength_um, density_g_cm3):
    directory = OUTPUT_DIR / "results"
    directory.mkdir(exist_ok=True)
    output = directory / "question3_mie_results.csv"
    columns = (*COLUMNS, *SPECIFIC_COLUMNS)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["real_index", "absorption_index", "wavelength_um",
                         "density_g_cm3", *columns])
        for (real_index, absorption_index), values in results.items():
            writer.writerows(
                (real_index, absorption_index, wavelength_um, density_g_cm3, *row)
                for row in zip(*(values[column] for column in columns))
            )
    return output


def save_figures(results, wavelength_um, density_g_cm3):
    directory = OUTPUT_DIR / "figures"
    directory.mkdir(exist_ok=True)
    with matplotlib.rc_context({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Liberation Serif", "STIXGeneral"],
        "mathtext.fontset": "stix", "pdf.fonttype": 42,
    }):
        for panels, filename, figsize in (
            (("specific", "g"), "question3_mie_radius", (9, 8.8)),
            (("specific",), "question3_specific_cross_sections", (9, 5.8)),
            (("g",), "question3_asymmetry", (9, 5.4)),
        ):
            fig, axes = plt.subplots(len(panels), 1, sharex=True,
                                     figsize=figsize, squeeze=False)
            for index, panel in enumerate(panels):
                ax = axes[index, 0]
                for real_index, absorption_index, color in CASES:
                    values = results[(real_index, absorption_index)]
                    radius = values["radius_um"]
                    label = rf"$m={real_index:g}-{absorption_index:g}i$"
                    if panel == "specific":
                        ax.loglog(radius, values[SPECIFIC_COLUMNS[0]],
                                  color=color, linewidth=1.8, label=label + ", extinction")
                        if absorption_index > 0.0:
                            ax.loglog(radius, values[SPECIFIC_COLUMNS[2]],
                                      color=color, linestyle="--", linewidth=1.8,
                                      label=label + ", absorption")
                    else:
                        ax.semilogx(radius, values["g"], color=color,
                                    linewidth=1.8, label=label)
                if panel == "specific":
                    # Zero has no logarithm. Do not substitute an artificial
                    # positive floor for the nonabsorbing material's residual.
                    handles, labels = ax.get_legend_handles_labels()
                    handles.append(Line2D([], [], color="C0", linestyle="--"))
                    labels.append(r"$m=1.5$: absorption $\simeq0$ (not on log axis)")
                    ax.legend(handles, labels, fontsize=10, ncol=2, loc="upper right")
                    ax.set(ylim=(0.1, 100),
                           ylabel=r"Specific cross-section, $\beta$ (m$^2$ g$^{-1}$)")
                    title = "Specific extinction and absorption"
                else:
                    ax.legend(fontsize=11, loc="lower right")
                    ax.set(ylim=(-0.015, 1.015), ylabel=r"Asymmetry factor, $g$")
                    title = "Scattering directionality"
                ax.set_xlim(0.001, 10.0)
                ax.set_title((f"({chr(97 + index)}) " if len(panels) > 1 else "")
                             + title, fontsize=15)
                ax.tick_params(labelsize=12)
                ax.xaxis.label.set_size(14)
                ax.yaxis.label.set_size(14)
                ax.grid(which="major", alpha=0.25)
            axes[-1, 0].set_xlabel(r"Particle radius, $r$ ($\mu$m)")
            fig.suptitle(
                rf"$\lambda={wavelength_um:g}\,\mu\mathrm{{m}}$; "
                rf"$\rho_p={density_g_cm3:g}\,\mathrm{{g\,cm^{{-3}}}}$ "
                "(assumed density)", fontsize=12, y=0.995,
            )
            fig.tight_layout(rect=(0, 0, 1, 0.96 if len(panels) > 1 else 0.94))
            for extension in ("png", "pdf"):
                fig.savefig(directory / f"{filename}.{extension}", dpi=300)
            plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wavelength-um", type=float, default=WAVELENGTH_UM)
    parser.add_argument("--density-g-cm3", type=float, default=DENSITY_G_CM3)
    args = parser.parse_args()
    print(f"Wavelength={args.wavelength_um:g} um; "
          f"assumed density={args.density_g_cm3:g} g/cm^3")
    results = run_sweep(args.wavelength_um, args.density_g_cm3)
    output = save_results(results, args.wavelength_um, args.density_g_cm3)
    save_figures(results, args.wavelength_um, args.density_g_cm3)
    for case, values in results.items():
        peak = int(np.argmax(values[SPECIFIC_COLUMNS[0]]))
        print(f"m={case[0]}-{case[1]}i: peak specific extinction="
              f"{values[SPECIFIC_COLUMNS[0]][peak]:.8f} m^2/g "
              f"at r={values['radius_um'][peak]:.8f} um")
        for radius in (0.001, 0.01, 0.1, 1, 4, 10):
            sample = int(np.argmin(abs(values["radius_um"] - radius)))
            print(f"  r={values['radius_um'][sample]:.6f}: "
                  f"beta_ext={values[SPECIFIC_COLUMNS[0]][sample]:.8g}, "
                  f"beta_abs={values[SPECIFIC_COLUMNS[2]][sample]:.8g}, "
                  f"g={values['g'][sample]:.8g}")
    print(f"Saved {len(RADII_UM) * len(CASES)} cases: {output}")
    print(f"Saved combined and separate PNG/PDF figures: {OUTPUT_DIR / 'figures'}")


if __name__ == "__main__":
    main()
