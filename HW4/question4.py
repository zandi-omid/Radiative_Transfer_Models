"""Replot and assess Question 4 from the saved HW4 CSV without rewriting it."""

import csv
from types import SimpleNamespace

import numpy as np

from HW4.hw4 import CASES, G, MU0, OUTPUT_DIR, SURFACE_ALBEDOS, save_energy_partition
from HW4.tests.test_hw4 import reflecting_boundary_reference
from models.liou_two_stream import liou_reflect_transmit


def load_results():
    """Load all 44 saved cases and normalize fluxes by pi*F0*mu0=mu0."""
    with (OUTPUT_DIR / "hw4_results.csv").open(newline="") as handle:
        rows = [{key: float(value) for key, value in row.items()}
                for row in csv.DictReader(handle)]
    if len(rows) != len(CASES) * len(SURFACE_ALBEDOS):
        raise ValueError("Expected 44 HW4 rows")
    if not all(row["g"] == G and row["mu0"] == MU0 for row in rows):
        raise ValueError("Saved HW4 parameters do not match the assignment")
    results = {}
    for case in CASES:
        group = sorted(
            (row for row in rows if (row["tau"], row["omega0"]) == case),
            key=lambda row: row["surface_albedo"],
        )
        np.testing.assert_allclose(
            [row["surface_albedo"] for row in group], SURFACE_ALBEDOS,
            rtol=0.0, atol=1e-14,
        )
        results[case] = SimpleNamespace(
            toa_reflectivity=np.array([row["toa_reflectivity"] for row in group]),
            surface_downward_flux=np.array([row["surface_downward_flux"] / MU0 for row in group]),
            surface_net_flux=np.array([row["surface_net_flux"] / MU0 for row in group]),
            atmospheric_absorption=np.array([row["atmospheric_absorption"] / MU0 for row in group]),
        )
    return results


def main():
    results = load_results()
    arithmetic_residual = 0.0
    stored_absorption_discrepancy = 0.0
    independent_errors = np.zeros(4)
    print("| tau | omega0 | R(As=0) | R(As=0.5) | R(As=1) | Delta R |")
    print("|---:|---:|---:|---:|---:|---:|")
    for (tau, omega0), result in results.items():
        inferred = 1.0 - result.toa_reflectivity - result.surface_net_flux
        fractions = np.array([result.toa_reflectivity, result.surface_net_flux, inferred])
        if not np.all(np.isfinite(fractions)):
            raise RuntimeError("Nonfinite saved energy fraction")
        arithmetic_residual = max(
            arithmetic_residual, float(np.max(np.abs(np.sum(fractions, axis=0) - 1.0)))
        )
        stored_absorption_discrepancy = max(
            stored_absorption_discrepancy,
            float(np.max(np.abs(inferred - result.atmospheric_absorption))),
        )
        r = result.toa_reflectivity
        print(f"| {tau:g} | {omega0:g} | {r[0]:.9f} | {r[5]:.9f} | {r[-1]:.9f} | {r[-1]-r[0]:.9f} |")
        black_r, black_t = liou_reflect_transmit(tau, omega0, MU0, G)
        np.testing.assert_allclose(
            [r[0], result.surface_net_flux[0]], [black_r, black_t], rtol=0, atol=1e-13
        )
        np.testing.assert_allclose(result.surface_net_flux[-1], 0.0, rtol=0, atol=1e-14)
        if omega0 == 0.99999:
            if not np.all((inferred >= -1e-12) & (inferred < 0.001)):
                raise RuntimeError("Nearly conservative absorption exceeds expected scale")
        elif not np.all(inferred > 0.2):
            raise RuntimeError("Absorbing cloud does not show appreciable absorption")
        print(f"  Absorption range: {np.min(inferred):.9g} to {np.max(inferred):.9g}")
        if np.any(fractions < -1e-12) or np.any(fractions > 1.0 + 1e-12):
            print(f"  WARNING: nonphysical fractions retained; min R={np.min(r):.9g}")
        # Independent numerical boundary solve and volume absorption integral.
        for index, albedo in enumerate(SURFACE_ALBEDOS):
            reference = reflecting_boundary_reference(tau, albedo, omega0, MU0, G)
            actual = [r[index], result.surface_downward_flux[index],
                      result.surface_net_flux[index], inferred[index]]
            independent_errors = np.maximum(
                independent_errors, np.abs(np.asarray(actual) - reference)
            )
    if arithmetic_residual > 1e-12 or stored_absorption_discrepancy > 1e-12:
        raise RuntimeError("Saved energy fractions fail closure checks")
    if np.max(independent_errors) > 2e-10:
        raise RuntimeError("Saved results disagree with the independent ODE solver")
    print(f"Arithmetic closure residual: {arithmetic_residual:.9e}")
    print(f"Stored absorption vs inferred absorption: {stored_absorption_discrepancy:.9e}")
    print("Independent ODE absolute errors (R, down, net, absorption):", independent_errors)
    save_energy_partition(results, SURFACE_ALBEDOS)
    print(f"Saved PNG and PDF: {OUTPUT_DIR / 'figures' / 'energy_partition_three_panels'}")


if __name__ == "__main__":
    main()
