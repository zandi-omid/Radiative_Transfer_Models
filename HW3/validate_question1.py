"""Independent numerical validation for HW3 Question 1.

This module solves the coupled Meador-Weaver two-stream differential
equations with an RK4 shooting method. It intentionally does not call the
closed-form ``reflect_transmit`` implementation.
"""

import csv
from pathlib import Path

import numpy as np

from models.two_stream import METHODS, reflect_transmit


MU_VALUES = (0.9, 0.5, 0.1)
ALBEDOS = (1.0, 0.9)
G = 0.75
OPTICAL_DEPTHS = np.linspace(0.0, 10.0, 201)
OUTPUT_DIR = Path(__file__).resolve().parent / "results"


def reference_coefficients(omega0, mu0, g, method):
    """Evaluate Table 1 independently of the production coefficient code."""
    if method == "eddington":
        gamma1 = (7.0 - omega0 * (4.0 + 3.0 * g)) / 4.0
        gamma2 = (omega0 * (4.0 - 3.0 * g) - 1.0) / 4.0
        gamma3 = (2.0 - 3.0 * g * mu0) / 4.0
    elif method == "quadrature":
        root_three = np.sqrt(3.0)
        gamma1 = root_three * (2.0 - omega0 * (1.0 + g)) / 2.0
        gamma2 = root_three * omega0 * (1.0 - g) / 2.0
        gamma3 = (1.0 - root_three * g * mu0) / 2.0
    else:
        raise ValueError(f"method must be one of {METHODS}")
    return gamma1, gamma2, gamma3, 1.0 - gamma3


def _derivative(depth, state, omega0, mu0, gamma):
    """Return derivatives for homogeneous and source-driven flux solutions."""
    gamma1, gamma2, gamma3, gamma4 = gamma
    homogeneous_plus, homogeneous_minus, particular_plus, particular_minus = state
    beam_source = omega0 * np.exp(-depth / mu0) / mu0
    return np.array(
        [
            gamma1 * homogeneous_plus - gamma2 * homogeneous_minus,
            gamma2 * homogeneous_plus - gamma1 * homogeneous_minus,
            gamma1 * particular_plus
            - gamma2 * particular_minus
            - gamma3 * beam_source,
            gamma2 * particular_plus
            - gamma1 * particular_minus
            + gamma4 * beam_source,
        ]
    )


def reference_reflect_transmit(
    tau,
    omega0=1.0,
    mu0=0.5,
    g=0.75,
    method="eddington",
    max_step=0.001,
):
    """Solve Eqs. (12)-(13) numerically with the two boundary conditions.

    The incident direct flux is normalized so that ``pi F mu0 = 1``. A
    homogeneous solution and a beam-forced particular solution are integrated
    together. Their linear combination is selected so that the downward
    diffuse flux vanishes at the top and the upward diffuse flux vanishes at
    the black lower boundary.
    """
    requested = np.asarray(tau, dtype=float)
    if np.any(~np.isfinite(requested)) or np.any(requested < 0.0):
        raise ValueError("tau must be finite and nonnegative")
    if max_step <= 0.0:
        raise ValueError("max_step must be positive")

    gamma = reference_coefficients(omega0, mu0, g, method)
    flat = requested.ravel()
    unique_depths, inverse = np.unique(flat, return_inverse=True)
    reflected_unique = np.empty_like(unique_depths)
    transmitted_unique = np.empty_like(unique_depths)

    # Columns 1-2 are the homogeneous solution; columns 3-4 are particular.
    state = np.array([1.0, 0.0, 0.0, 0.0])
    current_depth = 0.0
    for index, target_depth in enumerate(unique_depths):
        while current_depth < target_depth:
            step = min(max_step, target_depth - current_depth)
            k1 = _derivative(current_depth, state, omega0, mu0, gamma)
            k2 = _derivative(
                current_depth + step / 2.0,
                state + step * k1 / 2.0,
                omega0,
                mu0,
                gamma,
            )
            k3 = _derivative(
                current_depth + step / 2.0,
                state + step * k2 / 2.0,
                omega0,
                mu0,
                gamma,
            )
            k4 = _derivative(
                current_depth + step,
                state + step * k3,
                omega0,
                mu0,
                gamma,
            )
            state += step * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
            current_depth += step

        (
            homogeneous_plus,
            homogeneous_minus,
            particular_plus,
            particular_minus,
        ) = state
        reflected = -particular_plus / homogeneous_plus
        transmitted = (
            np.exp(-target_depth / mu0)
            + particular_minus
            + reflected * homogeneous_minus
        )
        reflected_unique[index] = reflected
        transmitted_unique[index] = transmitted

    reflected = reflected_unique[inverse].reshape(requested.shape)
    transmitted = transmitted_unique[inverse].reshape(requested.shape)
    if requested.ndim == 0:
        return float(reflected), float(transmitted)
    return reflected, transmitted


def run_validation():
    """Compare the analytical and independent ODE solutions over the sweep."""
    rows = []
    largest = {"R": (0.0, None), "T": (0.0, None)}
    for omega0 in ALBEDOS:
        for mu0 in MU_VALUES:
            for method in METHODS:
                analytical = reflect_transmit(
                    OPTICAL_DEPTHS, omega0=omega0, mu0=mu0, g=G, method=method
                )
                numerical = reference_reflect_transmit(
                    OPTICAL_DEPTHS, omega0=omega0, mu0=mu0, g=G, method=method
                )
                errors = [np.abs(a - n) for a, n in zip(analytical, numerical)]
                maxima = [float(np.max(error)) for error in errors]
                locations = [
                    float(OPTICAL_DEPTHS[np.argmax(error)]) for error in errors
                ]
                rows.append(
                    (
                        omega0,
                        mu0,
                        G,
                        method,
                        maxima[0],
                        locations[0],
                        maxima[1],
                        locations[1],
                    )
                )
                for quantity, error, location in zip(("R", "T"), maxima, locations):
                    if error > largest[quantity][0]:
                        largest[quantity] = (error, (omega0, mu0, method, location))
    return rows, largest


def main():
    rows, largest = run_validation()
    OUTPUT_DIR.mkdir(exist_ok=True)
    output = OUTPUT_DIR / "hw3_q1_validation_errors.csv"
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "omega0",
                "mu0",
                "g",
                "method",
                "max_abs_error_R",
                "tau_at_max_R",
                "max_abs_error_T",
                "tau_at_max_T",
            ]
        )
        writer.writerows(rows)

    print("Independent RK4 shooting validation")
    print("-----------------------------------")
    for quantity in ("R", "T"):
        error, case = largest[quantity]
        omega0, mu0, method, depth = case
        print(
            f"Maximum |delta {quantity}|: {error:.3e} "
            f"(omega0={omega0}, mu0={mu0}, {method}, tau={depth})"
        )
    print(f"Saved: {output}")


if __name__ == "__main__":
    main()
