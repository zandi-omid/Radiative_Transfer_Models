"""Compare the full saved radius sweep with the independent Mie series."""

import csv

import numpy as np

from HW5.hw5_question3 import OUTPUT_DIR
from HW5.reference_mie import reference_efficiencies


def main():
    data = np.genfromtxt(OUTPUT_DIR / "results/question3_mie_results.csv",
                         delimiter=",", names=True)
    errors = np.empty((len(data), 5))
    for index, row in enumerate(data):
        qe, qs, g = reference_efficiencies(
            row["x"], row["real_index"], row["absorption_index"]
        )
        # Independent dimensional conversion: m^2 cross-section / gram mass.
        radius_m = row["radius_um"] * 1e-6
        density_g_m3 = row["density_g_cm3"] * 1e6
        area_m2 = np.pi * radius_m**2
        mass_g = (4 / 3) * np.pi * radius_m**3 * density_g_m3
        expected = np.array([qe, qs, g, qe * area_m2 / mass_g,
                             (qe - qs) * area_m2 / mass_g])
        actual = np.array([row["qext"], row["qscat"], row["g"],
                           row["specific_extinction_m2_g"],
                           row["specific_absorption_m2_g"]])
        errors[index] = abs(actual - expected)
    quantities = ("qext", "qscat", "g", "specific_extinction_m2_g",
                  "specific_absorption_m2_g")
    maxima = np.max(errors, axis=0)
    if not np.all(np.isfinite(errors)):
        raise RuntimeError("Independent reference produced nonfinite differences")
    if np.any(maxima[:3] > 1e-8) or np.any(maxima[3:] > 1e-7):
        raise RuntimeError(f"Radius-sweep reference comparison failed: {maxima}")
    output = OUTPUT_DIR / "results/question3_validation.csv"
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["quantity", "max_absolute_error", "radius_um",
                         "real_index", "absorption_index"])
        for column, quantity in enumerate(quantities):
            row = data[np.argmax(errors[:, column])]
            writer.writerow([quantity, maxima[column], row["radius_um"],
                             row["real_index"], row["absorption_index"]])
            print(f"Maximum absolute error {quantity}: {maxima[column]:.9e}")
    print(f"Independently checked all {len(data)} cases; saved {output}")


if __name__ == "__main__":
    main()
