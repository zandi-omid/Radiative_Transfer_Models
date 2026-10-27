"""Liou (1974) two-stream discrete-ordinate reflectivity and transmissivity.

The implementation follows Eqs. (5)-(13) for a homogeneous layer with no
incoming diffuse radiation and a black lower boundary. Optical depth increases
downward, and ``mu0`` is the positive cosine of the solar zenith angle.
"""

import numpy as np


DISCRETE_COSINE = 1.0 / np.sqrt(3.0)


def liou_coefficients(omega0: float, g: float):
    """Return the first two phase moments and the Liou eigenvalue.

    For the first-order Henyey-Greenstein expansion in Liou Eq. (2),
    ``omega_bar_0 = omega0`` and ``omega_bar_1 = 3 * omega0 * g``.
    """
    if not 0.0 <= omega0 <= 1.0:
        raise ValueError("omega0 must be within [0, 1]")
    if not -1.0 <= g <= 1.0:
        raise ValueError("g must be within [-1, 1]")

    omega_bar_0 = omega0
    omega_bar_1 = 3.0 * omega0 * g
    eigenvalue_squared = (
        (1.0 - omega_bar_0)
        * (1.0 - omega_bar_1 * DISCRETE_COSINE**2)
        / DISCRETE_COSINE**2
    )
    if eigenvalue_squared < -1.0e-14:
        raise ValueError(f"discrete ordinate produced negative k^2: {eigenvalue_squared}")
    eigenvalue = np.sqrt(max(0.0, eigenvalue_squared))
    return omega_bar_0, omega_bar_1, eigenvalue


def liou_reflect_transmit(tau, omega0=1.0, mu0=0.5, g=0.75):
    """Return Liou two-stream reflectivity and total transmissivity.

    The nonconservative calculation evaluates Liou (1974), Eqs. (5)-(13),
    with incident irradiance ``F0 = 1``. Total transmission includes both the
    diffuse downward flux and the attenuated direct solar beam in Eq. (6).
    The conservative branch uses the exact zero-eigenvalue limit.
    """
    optical_depth = np.asarray(tau, dtype=float)
    if np.any(~np.isfinite(optical_depth)) or np.any(optical_depth < 0.0):
        raise ValueError("tau must be finite and nonnegative")
    if not 0.0 < mu0 <= 1.0:
        raise ValueError("mu0 must be within (0, 1] (positive downward cosine)")

    omega_bar_0, omega_bar_1, eigenvalue = liou_coefficients(omega0, g)
    mu1 = DISCRETE_COSINE
    direct_beam = np.exp(-optical_depth / mu0)

    if omega0 == 1.0:
        # Exact k -> 0 limit of the Liou discrete-ordinate equations. These
        # coefficients are the Gaussian-quadrature closure written in flux form.
        gamma1 = np.sqrt(3.0) * (1.0 - g) / 2.0
        gamma3 = (1.0 - np.sqrt(3.0) * g * mu0) / 2.0
        reflected = (
            gamma1 * optical_depth
            + (gamma3 - gamma1 * mu0) * (-np.expm1(-optical_depth / mu0))
        ) / (1.0 + gamma1 * optical_depth)
        transmitted = 1.0 - reflected
    else:
        resonance = 1.0 - (eigenvalue * mu0) ** 2
        if abs(resonance) < 1.0e-14:
            raise ValueError("singular Liou case: 1 - (k*mu0)^2 is zero")

        def eigenfunction(direction):
            return (
                omega_bar_0
                - omega_bar_1
                * (1.0 - omega_bar_0)
                * direction
                / eigenvalue
            ) / (1.0 + direction * eigenvalue)

        def beam_particular(direction):
            return (
                mu0
                * (direction - mu0)
                * (
                    omega_bar_0
                    - omega_bar_1
                    * (1.0 - omega_bar_0)
                    * direction
                    * mu0
                )
                / (4.0 * mu1**2 * resonance)
            )

        w_up = eigenfunction(mu1)
        w_down = eigenfunction(-mu1)
        z_up = beam_particular(mu1)
        z_down = beam_particular(-mu1)
        layer_attenuation = np.exp(-eigenvalue * optical_depth)

        a_plus = w_down + w_up * layer_attenuation
        a_minus = w_down - w_up * layer_attenuation
        c_plus = -(z_down + z_up * direct_beam)
        c_minus = -(z_down - z_up * direct_beam)
        l_plus = 0.5 * (c_plus / a_plus + c_minus / a_minus)
        l_minus_scaled = 0.5 * (c_plus / a_plus - c_minus / a_minus)

        # Eq. (5) at tau=0. The common factor pi*mu0*F0 is normalized out.
        upward_top = 2.0 * mu1 * (
            l_plus * w_up + l_minus_scaled * layer_attenuation * w_down + z_up
        )

        # Eq. (6) at tau=tau_N. l_minus_scaled omits the exp(-k*tau_N)
        # factor from Eq. (8), preventing an unnecessary exp(+k*tau_N).
        diffuse_downward_bottom = 2.0 * mu1 * (
            l_plus * w_down * layer_attenuation + l_minus_scaled * w_up
            + z_down * direct_beam
        )
        reflected = upward_top / mu0
        transmitted = direct_beam + diffuse_downward_bottom / mu0

    if optical_depth.ndim == 0:
        return float(reflected), float(transmitted)
    return reflected, transmitted
