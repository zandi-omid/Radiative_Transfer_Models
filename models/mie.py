"""Batch access to the course-supplied DMIESS homogeneous-sphere calculation."""

import io
from pathlib import Path
import shutil
import subprocess

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCES = (ROOT / "HW5/fortran/mie_batch.f90", ROOT / "HW5/fortran/dmiess.f")
EXECUTABLE = ROOT / "HW5/build/mie_batch"
COLUMNS = ("x", "radius_um", "qext", "qscat", "qabs", "g", "omega0")


def mass_specific_cross_section(efficiency, radius_um, density_g_cm3):
    """Convert dimensionless Q to cross-section per mass in m^2/g.

    C=pi*r^2*Q and mass=(4/3)*pi*r^3*rho give 3Q/(4*rho*r).
    With r in micrometers and rho in g/cm^3, the unit conversions cancel
    numerically, giving 0.75*Q/(rho*r). Preserve signed roundoff residuals.
    """
    efficiency = np.asarray(efficiency, dtype=float)
    radius = np.asarray(radius_um, dtype=float)
    if np.any(~np.isfinite(efficiency)):
        raise ValueError("efficiency must be finite")
    if np.any(~np.isfinite(radius)) or np.any(radius <= 0):
        raise ValueError("radius_um must be positive and finite")
    if not np.isfinite(density_g_cm3) or density_g_cm3 <= 0:
        raise ValueError("density_g_cm3 must be positive and finite")
    specific = 0.75 * efficiency / (density_g_cm3 * radius)
    return float(specific) if specific.ndim == 0 else specific


def build_mie():
    """Compile when needed, with bounds checks and floating-point traps."""
    if EXECUTABLE.exists() and all(
        source.stat().st_mtime_ns <= EXECUTABLE.stat().st_mtime_ns
        for source in SOURCES
    ):
        return EXECUTABLE
    compiler = shutil.which("gfortran")
    if compiler is None:
        raise RuntimeError("Install/load gfortran to compile the supplied Mie routine")
    EXECUTABLE.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [compiler, "-std=legacy", "-ffixed-line-length-none", "-O0", "-g",
         "-fcheck=all", "-ffpe-trap=invalid,zero,overflow",
         *(str(source) for source in SOURCES), "-o", str(EXECUTABLE)],
        check=True, capture_output=True, text=True,
    )
    return EXECUTABLE


def mie_sweep(size_parameters, real_index, absorption_index, wavelength_um=0.55):
    """Return arrays of efficiencies and g for homogeneous spheres in air.

    Positive ``absorption_index`` is passed to DMIESS, which constructs
    m=real_index-i*absorption_index. x must be strictly positive: at x=0
    there is no particle and the ratio defining g is undefined.
    """
    x = np.atleast_1d(np.asarray(size_parameters, dtype=float))
    if x.ndim != 1 or x.size == 0 or np.any(~np.isfinite(x)) or np.any(x <= 0):
        raise ValueError("size_parameters must be a nonempty positive finite 1D array")
    for name, value in (("real_index", real_index), ("wavelength_um", wavelength_um)):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be positive and finite")
    if not np.isfinite(absorption_index) or absorption_index < 0:
        raise ValueError("absorption_index must be nonnegative and finite")
    input_text = "".join(
        f"{value:.17g} {real_index:.17g} {absorption_index:.17g} {wavelength_um:.17g}\n"
        for value in x
    )
    process = subprocess.run(
        [str(build_mie())], input=input_text, text=True, capture_output=True,
        check=True, timeout=120,
    )
    output = np.loadtxt(io.StringIO(process.stdout), ndmin=2)
    if output.shape != (x.size, len(COLUMNS)) or not np.all(np.isfinite(output)):
        raise RuntimeError("DMIESS returned incomplete or nonfinite output")
    np.testing.assert_allclose(output[:, 0], x, rtol=1e-14, atol=0)
    if np.any(output[:, 2:4] < 0) or np.any(output[:, 4] < -1e-10):
        raise RuntimeError("DMIESS produced negative efficiencies; results were not clipped")
    if np.any(np.abs(output[:, 5]) > 1 + 1e-10):
        raise RuntimeError("DMIESS produced g outside [-1,1]; results were not clipped")
    return dict(zip(COLUMNS, output.T))


def mie_phase_function(radius_um, real_index, absorption_index, wavelength_um=0.55,
                       angles_degrees=None):
    """Return unpolarized Mie P(theta), using the supplied 4*pi normalization.

    Angles are relative to the incident ray: zero is forward, 180 is backward.
    DMIESS computes angles in [0,90] and their reflected counterparts. Up to
    183 distinct folded angles are supported per call. Default sampling is
    0 through 180 degrees inclusive in one-degree steps.
    """
    if not np.isfinite(radius_um) or radius_um <= 0:
        raise ValueError("radius_um must be positive and finite")
    for name, value in (("real_index", real_index), ("wavelength_um", wavelength_um)):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be positive and finite")
    if not np.isfinite(absorption_index) or absorption_index < 0:
        raise ValueError("absorption_index must be nonnegative and finite")
    angles = np.arange(181, dtype=float) if angles_degrees is None else np.atleast_1d(
        np.asarray(angles_degrees, dtype=float)
    )
    if (angles.ndim != 1 or angles.size == 0 or np.any(~np.isfinite(angles))
            or np.any((angles < 0) | (angles > 180))):
        raise ValueError("angles_degrees must be a nonempty finite 1D array in [0,180]")
    folded, inverse = np.unique(np.minimum(angles, 180.0 - angles), return_inverse=True)
    if folded.size > 183:
        raise ValueError("DMIESS supports at most 183 distinct folded angles per call")
    x = 2 * np.pi * radius_um / wavelength_um
    input_text = (
        f"{x:.17g} {real_index:.17g} {absorption_index:.17g} {wavelength_um:.17g}\n"
        f"{folded.size}\n" + " ".join(f"{angle:.17g}" for angle in folded) + "\n"
    )
    process = subprocess.run(
        [str(build_mie()), "--phase"], input=input_text, text=True,
        capture_output=True, check=True, timeout=120,
    )
    lines = process.stdout.splitlines()
    if len(lines) != folded.size + 1:
        raise RuntimeError("DMIESS returned an incomplete phase-function table")
    metadata = np.fromstring(lines[0], sep=" ")
    table = np.loadtxt(io.StringIO("\n".join(lines[1:])), ndmin=2)
    if (metadata.shape != (7,) or table.shape != (folded.size, 3)
            or not np.all(np.isfinite(metadata)) or not np.all(np.isfinite(table))):
        raise RuntimeError("DMIESS returned malformed or nonfinite phase-function output")
    np.testing.assert_allclose(table[:, 0], folded, rtol=0, atol=1e-12)
    phase = np.where(angles > 90, table[inverse, 2], table[inverse, 1])
    if np.any(phase < 0):
        raise RuntimeError("DMIESS returned a negative phase function; values were not clipped")
    return {
        **dict(zip(COLUMNS, metadata)),
        "angle_deg": angles,
        "phase_4pi": phase,
        "phase_per_sr": phase / (4 * np.pi),
    }
