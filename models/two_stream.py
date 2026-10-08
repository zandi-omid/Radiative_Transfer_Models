"""Meador-Weaver two-stream reflectivity and transmissivity.

Implements the Eddington and quadrature closures for a homogeneous layer over
a black lower boundary. The solar incidence cosine ``mu0`` is positive.
"""

import numpy as np


METHODS = ("eddington", "quadrature")


def coefficients(omega0: float, mu0: float, g: float, method: str):
    """Return gamma1..gamma4, alpha1, alpha2, and k for a closure."""
    if not 0.0 <= omega0 <= 1.0:
        raise ValueError("omega0 must be within [0, 1]")
    if not 0.0 < mu0 <= 1.0:
        raise ValueError("mu0 must be within (0, 1] (positive downward cosine)")
    if not -1.0 <= g <= 1.0:
        raise ValueError("g must be within [-1, 1]")

    if method == "eddington":
        gamma1 = (7.0 - omega0 * (4.0 + 3.0 * g)) / 4.0
        gamma2 = -(1.0 - omega0 * (4.0 - 3.0 * g)) / 4.0
        gamma3 = (2.0 - 3.0 * g * mu0) / 4.0
    elif method == "quadrature":
        gamma1 = np.sqrt(3.0) / 2.0 * (2.0 - omega0 * (1.0 + g))
        gamma2 = np.sqrt(3.0) / 2.0 * omega0 * (1.0 - g)
        gamma3 = (1.0 - np.sqrt(3.0) * g * mu0) / 2.0
    else:
        raise ValueError(f"method must be one of {METHODS}")

    gamma4 = 1.0 - gamma3
    alpha1 = gamma1 * gamma4 + gamma2 * gamma3
    alpha2 = gamma1 * gamma3 + gamma2 * gamma4
    k_squared = gamma1**2 - gamma2**2
    if k_squared < -1.0e-14:
        raise ValueError(f"closure produced negative k^2: {k_squared}")
    k = np.sqrt(max(0.0, k_squared))
    return gamma1, gamma2, gamma3, gamma4, alpha1, alpha2, k


def reflect_transmit(
    tau,
    omega0=1.0,
    mu0=0.5,
    g=0.75,
    method="eddington",
):
    """Return reflectivity and transmissivity for optical depth ``tau``.

    The nonconservative branch uses exponent-scaled forms of Meador and Weaver
    (1980), Eqs. (14)-(15), with coefficients from Eqs. (16)-(18) and Table 1.
    The conservative branch uses the analytical ``k -> 0`` limit (Eq. 24).
    Scalar input produces scalar output; array-like input produces arrays.
    """
    optical_depth = np.asarray(tau, dtype=float)
    if np.any(~np.isfinite(optical_depth)) or np.any(optical_depth < 0.0):
        raise ValueError("tau must be finite and nonnegative")

    gamma1, _, gamma3, gamma4, alpha1, alpha2, k = coefficients(
        omega0,
        mu0,
        g,
        method,
    )
    direct_beam = np.exp(-optical_depth / mu0)

    if omega0 == 1.0:
        # Exact conservative-scattering limit; no artificial absorption.
        reflected = (
            gamma1 * optical_depth
            + (gamma3 - gamma1 * mu0) * (-np.expm1(-optical_depth / mu0))
        ) / (1.0 + gamma1 * optical_depth)
        transmitted = 1.0 - reflected
    else:
        attenuation = np.exp(-k * optical_depth)
        attenuation_squared = attenuation**2
        scaled_d = (k + gamma1) + (k - gamma1) * attenuation_squared
        resonance_factor = 1.0 - (k * mu0) ** 2
        if abs(resonance_factor) < 1.0e-14:
            raise ValueError("singular two-stream case: 1 - (k*mu0)^2 is zero")
        denominator = resonance_factor * scaled_d

        reflected_numerator = (
            (1.0 - k * mu0) * (alpha2 + k * gamma3)
            - (1.0 + k * mu0)
            * (alpha2 - k * gamma3)
            * attenuation_squared
            - 2.0
            * k
            * (gamma3 - alpha2 * mu0)
            * direct_beam
            * attenuation
        )
        transmission_correction = (
            direct_beam
            * (
                (1.0 + k * mu0) * (alpha1 + k * gamma4)
                - (1.0 - k * mu0)
                * (alpha1 - k * gamma4)
                * attenuation_squared
            )
            - 2.0 * k * (gamma4 + alpha1 * mu0) * attenuation
        )
        reflected = omega0 * reflected_numerator / denominator
        transmitted = direct_beam - omega0 * transmission_correction / denominator

    if optical_depth.ndim == 0:
        return float(reflected), float(transmitted)
    return reflected, transmitted
