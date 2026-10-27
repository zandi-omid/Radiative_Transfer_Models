"""Reproduce HW5 Question 2: extinction efficiency and g versus size."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from models.mie import COLUMNS, mie_sweep


OUTPUT_DIR = Path(__file__).resolve().parent
SIZE_PARAMETERS = np.arange(1, 1001, dtype=float) * 0.05
CASES = ((1.3873, 0.0066, "C0"), (1.2993, 0.037, "C1"))
WAVELENGTH_UM = 0.55


def run_sweep():
    return {
        (real_index, absorption_index): mie_sweep(
            SIZE_PARAMETERS, real_index, absorption_index, WAVELENGTH_UM
        )
        for real_index, absorption_index, _ in CASES
    }


def save_results(results):
    result_dir = OUTPUT_DIR / "results"
    result_dir.mkdir(exist_ok=True)
    output = result_dir / "question2_mie_results.csv"
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["real_index", "absorption_index", "wavelength_um", *COLUMNS])
        for (real_index, absorption_index), values in results.items():
            writer.writerows(
                (real_index, absorption_index, WAVELENGTH_UM, *row)
                for row in zip(*(values[name] for name in COLUMNS))
            )
    return output


def save_figures(results):
    figure_dir = OUTPUT_DIR / "figures"
    figure_dir.mkdir(exist_ok=True)
    with matplotlib.rc_context({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Liberation Serif", "STIXGeneral"],
        "mathtext.fontset": "stix", "pdf.fonttype": 42,
    }):
        # Match the assignment's stacked-panel layout and add separate panels
        # for reuse on presentation slides.
        for selected, filename, figsize in (
            (("qext", "g"), "question2_mie_size_parameter", (9, 8)),
            (("qext",), "question2_extinction", (9, 5)),
            (("g",), "question2_asymmetry", (9, 5)),
        ):
            fig, axes = plt.subplots(len(selected), 1, figsize=figsize,
                                     sharex=True, squeeze=False)
            for index, quantity in enumerate(selected):
                ax = axes[index, 0]
                for real_index, absorption_index, color in CASES:
                    values = results[(real_index, absorption_index)]
                    ax.plot(values["x"], values[quantity], color=color, linewidth=1.9,
                            label=rf"$m={real_index:g}-{absorption_index:g}i$")
                ax.set(xlim=(0, 50), ylabel=(
                    r"Extinction efficiency, $Q_{\mathrm{ext}}$" if quantity == "qext"
                    else r"Asymmetry factor, $g$"
                ))
                ax.set_ylim((0, 5) if quantity == "qext" else (0, 1))
                if quantity == "qext":
                    ax.axhline(2.0, color="0.45", linestyle=":", linewidth=1)
                ax.set_title(
                    (f"({chr(97 + index)}) " if len(selected) > 1 else "")
                    + ("Mie extinction" if quantity == "qext" else "Scattering directionality"),
                    fontsize=16,
                )
                ax.grid(alpha=0.23)
                ax.tick_params(labelsize=12)
                ax.xaxis.label.set_size(14)
                ax.yaxis.label.set_size(14)
            axes[-1, 0].set_xlabel(r"Size parameter, $x=2\pi r/\lambda$")
            handles, labels = axes[0, 0].get_legend_handles_labels()
            fig.legend(handles, labels, loc="upper center", ncol=2,
                       bbox_to_anchor=(0.5, 0.995), fontsize=12, frameon=False)
            fig.tight_layout(rect=(0, 0, 1, 0.94 if len(selected) > 1 else 0.90))
            for extension in ("png", "pdf"):
                fig.savefig(figure_dir / f"{filename}.{extension}", dpi=300)
            plt.close(fig)


def main():
    results = run_sweep()
    csv_path = save_results(results)
    save_figures(results)
    for (real_index, absorption_index), values in results.items():
        peak = int(np.argmax(values["qext"]))
        print(f"m={real_index}-{absorption_index}i:")
        print(f"  Qext maximum={values['qext'][peak]:.8f} at x={values['x'][peak]:.2f}")
        for x in (0.05, 1, 5, 10, 25, 50):
            index = int(round(x / 0.05)) - 1
            print(f"  x={x:5.2f}: Qext={values['qext'][index]:.8f}, "
                  f"Qscat={values['qscat'][index]:.8f}, g={values['g'][index]:.8f}")
    print(f"Saved 2,000 cases: {csv_path}")
    print(f"Saved combined and individual PNG/PDF figures: {OUTPUT_DIR / 'figures'}")


if __name__ == "__main__":
    main()
