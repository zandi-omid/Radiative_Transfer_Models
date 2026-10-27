"""Liou two-stream adding method for a Lambertian reflecting surface."""

from dataclasses import dataclass

import numpy as np

from models.liou_two_stream import (
    DISCRETE_COSINE,
    liou_coefficients,
    liou_reflect_transmit,
)


@dataclass(frozen=True)
class SurfaceResult:
    """Reflectivity is a fraction; fluxes use the supplied incident flux units."""

    toa_reflectivity: object
    surface_downward_flux: object
    surface_upward_flux: object
    surface_net_flux: object
    atmospheric_absorption: object
    reflection_denominator: object


def diffuse_reflect_transmit(tau, omega0=1.0, g=0.75):
    """Return bottom-incident diffuse R and T from Liou's homogeneous system.

    With U'=a U-b D and D'=b U-a D, the transfer matrix is
    cosh(k*tau) I + sinh(k*tau) M/k. Imposing D(0)=0 and U(tau)=1
    gives R=b sinh(k*tau)/(k cosh(k*tau)+a sinh(k*tau)) and
    T=k/(k cosh(k*tau)+a sinh(k*tau)). The scaled evaluation below
    avoids growing exponentials and uses the exact k=0 limit.
    """
    depth = np.asarray(tau, dtype=float)
    if np.any(~np.isfinite(depth)) or np.any(depth < 0.0):
        raise ValueError("tau must be finite and nonnegative")
    w0, w1, k = liou_coefficients(omega0, g)
    a = (1.0 - (w0 + w1 * DISCRETE_COSINE**2) / 2.0) / DISCRETE_COSINE
    b = (w0 - w1 * DISCRETE_COSINE**2) / (2.0 * DISCRETE_COSINE)
    if k == 0.0:
        q = 2.0 * depth
    else:
        q = -np.expm1(-2.0 * k * depth) / k
    denominator = 2.0 + (a - k) * q
    reflected = b * q / denominator
    transmitted = 2.0 * np.exp(-k * depth) / denominator
    if depth.ndim == 0:
        return float(reflected), float(transmitted)
    return reflected, transmitted


def surface_transport(
    tau, surface_albedo, omega0=1.0, mu0=0.9, g=0.75, incident_flux=1.0
):
    """Sum all cloud-surface reflections with consistent direct/diffuse fluxes.

    ``incident_flux`` is the horizontal direct solar flux pi*F0*mu0.
    Its default of one reports flux fractions. Set it to mu0 for pi*F0=1.
    Optical depth and surface albedo may be broadcast-compatible arrays.
    """
    albedo = np.asarray(surface_albedo, dtype=float)
    if np.any(~np.isfinite(albedo)) or np.any((albedo < 0.0) | (albedo > 1.0)):
        raise ValueError("surface_albedo must be finite and within [0, 1]")
    if not np.isfinite(incident_flux) or incident_flux <= 0.0:
        raise ValueError("incident_flux must be finite and positive")
    if not 0.0 < mu0 <= 1.0:
        raise ValueError("mu0 must be within (0, 1]")
    diffuse_r, diffuse_t = diffuse_reflect_transmit(tau, omega0, g)
    # Liou's eigenfunction representation degenerates at zero scattering.
    if omega0 == 0.0:
        direct_r = np.zeros_like(np.asarray(tau, dtype=float))
        direct_t = np.exp(-np.asarray(tau, dtype=float) / mu0)
    else:
        direct_r, direct_t = liou_reflect_transmit(tau, omega0, mu0, g)

    denominator = 1.0 - albedo * diffuse_r
    if np.any(~np.isfinite(denominator)) or np.any(denominator <= 0.0):
        raise ValueError("cloud-surface reflection series does not converge")
    downward = incident_flux * direct_t / denominator
    upward = albedo * downward
    reflected = direct_r + diffuse_t * upward / incident_flux
    net = (1.0 - albedo) * downward
    # Absorb the initial beam and each subsequent upward diffuse pass.
    absorption = incident_flux * (1.0 - direct_r - direct_t)
    absorption += upward * (1.0 - diffuse_r - diffuse_t)
    values = (reflected, downward, upward, net, absorption, denominator)
    return SurfaceResult(
        *(float(value) if np.ndim(value) == 0 else value for value in values)
    )
