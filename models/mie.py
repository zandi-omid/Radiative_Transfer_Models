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
