"""Independent Bessel/Legendre validation of the saved Q4 phase functions."""

import csv

import numpy as np

from HW5.hw5_question4 import ANGLES_DEG, CASES, OUTPUT_DIR, phase_moments
from HW5.reference_mie import reference_phase_function


def main():
    saved = np.genfromtxt(OUTPUT_DIR / "results/question4_phase_function.csv",
                          delimiter=",", names=True)
    output = OUTPUT_DIR / "results/question4_validation.csv"
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["real_index", "absorption_index", "max_absolute_phase_error",
                         "max_relative_phase_error", "normalization_residual",
                         "g_moment_error"])
        for n, k, _ in CASES:
            rows = saved[(saved["real_index"] == n) & (saved["absorption_index"] == k)]
            np.testing.assert_array_equal(rows["angle_deg"], ANGLES_DEG)
            expected = reference_phase_function(rows["x"][0], n, k, rows["angle_deg"])
            np.testing.assert_allclose(rows["phase_4pi"], expected, rtol=1e-9, atol=1e-10)
            np.testing.assert_allclose(rows["phase_per_sr"], rows["phase_4pi"] / (4*np.pi))
            difference = abs(rows["phase_4pi"] - expected)
            maximum = float(np.max(difference))
            relative = float(np.max(difference / expected))
            normalization, angular_g, series_g = phase_moments(n, k)
            if abs(normalization - 1) > 5e-10 or abs(angular_g - series_g) > 5e-10:
                raise RuntimeError("Angular normalization/moment check failed")
            writer.writerow([n, k, maximum, relative,
                             abs(normalization - 1), abs(angular_g - series_g)])
            print(f"m={n}-{k}i: max absolute phase error={maximum:.9e}, "
                  f"max relative error={relative:.9e}")
            print(f"  normalization residual={abs(normalization-1):.9e}, "
                  f"g moment error={abs(angular_g-series_g):.9e}")
    print(f"Independently checked all 362 angular values; saved {output}")


if __name__ == "__main__":
    main()
