# Radiative Transfer Models

This repository contains Monte Carlo photon-transport models and experiments
for ATMO 656A: Atmospheric Radiation and Remote Sensing. The code follows
individual photons through a homogeneous, plane-parallel atmospheric layer and
estimates the fractions that are reflected, transmitted, or absorbed.

The same transport model is shared by the HW1 and HW2 experiments. HW1 uses
isotropic, nonabsorbing scattering. HW2 adds Henyey-Greenstein anisotropic
scattering, absorption through the single-scattering albedo, and a computational
comparison between the isotropic and anisotropic models.

## Physical Model

The bottom of the layer is at `z = 0`, and the top is at
`z = tau * l0`. The vertical coordinate is positive upward, so the incident
solar direction cosine is negative for a downward-moving photon.

The distance to the next interaction is exponentially distributed:

```text
l = -l0 ln(1 - xi),  xi ~ U(0, 1)
```

The photon then moves vertically by:

```text
delta_z = mu * l
```

A photon crossing the top boundary is reflected. A photon crossing the bottom
boundary is transmitted. At an interaction inside the layer, the photon
scatters with probability `omega0` and is absorbed with probability
`1 - omega0`. The model therefore checks:

```text
R + T + A = 1
```

For nonabsorbing simulations, `omega0 = 1`, `A = 0`, and `R + T = 1`.

### Scattering Direction

For isotropic scattering (`g = 0`), the new vertical direction cosine is
uniformly sampled from `[-1, 1]`.

For Henyey-Greenstein scattering, the cosine of the scattering angle relative
to the photon's current direction is sampled by inverse CDF:

```text
q = (1 - g^2) / (1 - g + 2 g xi)
cos(theta) = (1 + g^2 - q^2) / (2 g)
```

A uniform azimuth `phi` is then used to rotate that scattering angle into the
vertical coordinate system:

```text
mu_new = mu_old cos(theta)
         + sqrt(1 - mu_old^2) sin(theta) cos(phi)
```

Here, `cos(theta)` is relative to the incoming photon direction, while `mu_old`
and `mu_new` are relative to the vertical. The implementation explicitly uses
the isotropic limit when `g` is zero.

## Repository Structure

```text
Radiative_Transfer_Models/
|-- models/
|   `-- monte_carlo.py
|-- HW1/
|   |-- HW1.docx
|   |-- hw1_monte_carlo_photons.py
|   |-- hw1_problem2_independent_runs.py
|   |-- hw1_problem3_incident_direction.py
|   `-- hw1_problem4_optical_depth.py
|-- HW2/
|   |-- hw2_part1_hg_scattering.py
|   |-- hw2_full_experiments.py
|   `-- compare_isotropic_hg_performance.py
|-- tests/
|   `-- test_monte_carlo.py
`-- requirements.txt
```

`models/monte_carlo.py` defines the shared `MonteCarloTransport` class and the
HG sampling and direction-rotation functions. Each homework script configures
that model for a particular experiment.

## Installation

Python 3 with NumPy and Matplotlib is required.

```bash
git clone https://github.com/zandi-omid/Radiative_Transfer_Models.git
cd Radiative_Transfer_Models
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Run all commands below from the repository root so Python can import the
shared `models` package.

## HW1: Isotropic Scattering

HW1 uses `g = 0` and `omega0 = 1`.

```bash
MPLBACKEND=Agg python -m HW1.hw1_monte_carlo_photons
MPLBACKEND=Agg python -m HW1.hw1_problem2_independent_runs
MPLBACKEND=Agg python -m HW1.hw1_problem3_incident_direction
MPLBACKEND=Agg python -m HW1.hw1_problem4_optical_depth
```

These scripts cover:

1. Convergence with increasing photon number.
2. Ten independent simulations and sample standard deviations.
3. Reflectivity and transmissivity for `mu0 = -0.1, ..., -0.9`.
4. Reflectivity, transmissivity, scattering count, and computation time for
   optical depths `tau = 1, ..., 8`.

## HW2: HG Scattering and Absorption

The tested baseline uses `N = 10,000`, `tau = 4`, `mu0 = -0.7`, `g = 0.75`,
and `omega0 = 1`:

```bash
python -m HW2.hw2_part1_hg_scattering
```

With seed 42 and 24 available workers, the reference result is:

```text
R = 0.423300
T = 0.576700
A = 0.000000
R + T + A = 1.000000
```

Run the complete HW2 experiments with:

```bash
python -m HW2.hw2_full_experiments
```

This uses `N = 10,000` and `g = 0.75` throughout and evaluates
`omega0 = 1.00, 0.95, 0.85, 0.75`. It performs:

- Ten independent simulations for each albedo.
- An incident-angle sweep from `mu0 = -0.1` through `-0.9`.
- An optical-depth sweep from `tau = 1` through `8`.
- Conservation checks for `R + T + A`.
- High-resolution comparison figures for reflectivity, transmissivity, and
  absorptivity.

Numerical tables and figures are written to `HW2/results/`.

## Computational Comparison

The performance benchmark compares the original isotropic model (`g = 0`) and
the nonabsorbing HG model (`g = 0.75`) under the same conditions:

```bash
python -m HW2.compare_isotropic_hg_performance
```

The benchmark uses exactly `N = 10,000`, `mu0 = -0.7`, `omega0 = 1`, and
`tau = 1, ..., 8`. Each model is timed five times per optical depth with the
same seeds and worker count. Execution order alternates between models to
reduce ordering bias. The output includes raw timings, timing means and sample
standard deviations, and a wall-clock-time comparison figure with error bars.

## Multiprocessing and Reproducibility

Simulations use up to 24 worker processes:

```text
workers_used = min(requested_workers, available_CPU_count, photon_count)
```

`numpy.random.SeedSequence` creates an independent random-number stream for
each worker. A run is reproducible when its seed and worker count are unchanged.
The benchmark reports the actual worker count in its terminal output and figure
title.

## Tests

Run the automated checks with:

```bash
python -m unittest discover -s tests -v
```

The tests verify photon conservation, reproducibility, input validation, the
mean cosine of the HG samples, the direction rotation, and conservation when
absorption is enabled.

Generated PNG and numerical result files are excluded from Git and can be
recreated by running the corresponding experiment scripts.
