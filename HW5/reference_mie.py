"""Independent homogeneous-sphere series using SciPy spherical Bessel functions.

This validation calculation uses direct special functions instead of the
supplied DMIESS logarithmic-derivative recurrence or coated-sphere algebra.
SciPy is needed only for this reference calculation and its tests.
"""

import numpy as np
from scipy.special import spherical_jn, spherical_yn


def reference_efficiencies(x, real_index, absorption_index):
    """Return Qext, Qscat and g for a positive scalar size parameter.

    Use m=n-ik and the outgoing Hankel function of the second kind, matching
    the time/sign convention of the supplied Fortran routine.
    """
    if not np.isfinite(x) or x <= 0:
        raise ValueError("x must be positive and finite")
    m = complex(real_index, -absorption_index)
    orders = np.arange(1, int(np.ceil(x + 4 * np.cbrt(x) + 2)) + 1)
    psi_x = x * spherical_jn(orders, x)
    psi_x_prime = spherical_jn(orders, x) + x * spherical_jn(orders, x, derivative=True)
    psi_mx = m * x * spherical_jn(orders, m * x)
    psi_mx_prime = (spherical_jn(orders, m * x)
                    + m * x * spherical_jn(orders, m * x, derivative=True))
    xi_x = psi_x - 1j * x * spherical_yn(orders, x)
    xi_x_prime = psi_x_prime - 1j * (
        spherical_yn(orders, x) + x * spherical_yn(orders, x, derivative=True)
    )
    a = (m * psi_mx * psi_x_prime - psi_x * psi_mx_prime) / (
        m * psi_mx * xi_x_prime - xi_x * psi_mx_prime
    )
    b = (psi_mx * psi_x_prime - m * psi_x * psi_mx_prime) / (
        psi_mx * xi_x_prime - m * xi_x * psi_mx_prime
    )
    qext = 2 / x**2 * np.sum((2 * orders + 1) * (a + b).real)
    qscat = 2 / x**2 * np.sum((2 * orders + 1) * (abs(a)**2 + abs(b)**2))
    n = orders[:-1]
    adjacent = np.sum(n * (n + 2) / (n + 1) * (
        a[:-1] * a[1:].conj() + b[:-1] * b[1:].conj()
    ).real)
    cross = np.sum((2 * orders + 1) / (orders * (orders + 1)) * (a * b.conj()).real)
    g = 4 * (adjacent + cross) / (x**2 * qscat)
    return float(qext), float(qscat), float(g)


def main():
    """Compare all saved Question 2 rows with this independent series."""
    from HW5.hw5_question2 import OUTPUT_DIR

    data = np.genfromtxt(OUTPUT_DIR / "results/question2_mie_results.csv",
                         delimiter=",", names=True)
    errors = []
    for row in data:
        expected = reference_efficiencies(row["x"], row["real_index"], row["absorption_index"])
        actual = np.array([row["qext"], row["qscat"], row["g"]])
        errors.append(np.abs(actual - expected))
    maxima = np.max(errors, axis=0)
    print("Independent series comparison over", len(data), "cases")
    for name, error in zip(("Qext", "Qscat", "g"), maxima):
        print(f"Maximum absolute error {name}: {error:.9e}")
    if np.max(maxima) > 5e-7:
        raise RuntimeError("DMIESS disagrees with independent series")


if __name__ == "__main__":
    main()
