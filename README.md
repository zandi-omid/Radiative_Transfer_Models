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

#### Henyey-Greenstein Phase Function

The Henyey-Greenstein (HG) phase function describes the angular distribution
of scattered radiation:

$$
P(\Theta) = \frac{1}{4\pi}
\frac{1-g^2}{\left(1+g^2-2g\cos\Theta\right)^{3/2}}.
$$

The scattering angle $\Theta$ is measured relative to the photon's **current
incoming direction**, rather than relative to the atmospheric vertical. The
asymmetry factor $g$ controls the preferred scattering direction:

- $g=0$ gives isotropic scattering.
- $g>0$ favors forward scattering, with directions near the incoming photon
  direction receiving greater probability.
- HW2 uses $g=0.75$, representing strongly forward-peaked scattering.

#### Inverse-CDF Sampling

Define the scattering-angle cosine as

$$
\mu_s = \cos\Theta.
$$

The function `sample_hg_cos_theta()` samples $\mu_s$ by inverse-transform
sampling. It first draws $\xi\sim U(0,1)$ and computes

$$
q = \frac{1-g^2}{1-g+2g\xi},
$$

followed by

$$
\cos\Theta = \frac{1+g^2-q^2}{2g}.
$$

This is conceptually the same inverse-transform method used to sample the
exponential free path,

$$
l=-l_0\ln(1-\xi),
$$

but the HG inverse CDF samples a scattering angle instead of a path length.
For $g=0$, `sample_hg_cos_theta()` handles the isotropic limit separately by
sampling $\cos\Theta$ uniformly from $[-1,1]$. In the actual isotropic
transport path, `scatter_direction()` directly samples the new vertical
direction cosine from this uniform distribution, preserving the original HW1
behavior.

#### Random Azimuth

Sampling $\Theta$ fixes the opening angle of a cone around the incoming photon
direction, but the photon may scatter anywhere around that cone. The HG phase
function does not impose a preferred azimuthal direction, so the code
independently draws

$$
\phi=2\pi\xi_2, \qquad \xi_2\sim U(0,1).
$$

#### Rotation into the Vertical Coordinate System

Three direction cosines must be distinguished:

- $\cos\Theta$ is relative to the **old photon direction**.
- $\mu_{\mathrm{old}}$ is the old photon direction cosine relative to the
  atmospheric vertical.
- $\mu_{\mathrm{new}}$ is the new photon direction cosine relative to the
  atmospheric vertical.

The function `scatter_direction()` rotates the locally sampled scattering
direction into the global vertical coordinate system using

$$
\mu_{\mathrm{new}} =
\mu_{\mathrm{old}}\cos\Theta
+ \sqrt{1-\mu_{\mathrm{old}}^2}\,\sin\Theta\cos\phi,
$$

or equivalently,

$$
\mu_{\mathrm{new}} =
\mu_{\mathrm{old}}\cos\Theta
+ \sqrt{1-\mu_{\mathrm{old}}^2}
  \sqrt{1-\cos^2\Theta}\cos\phi.
$$

Only $\mu_{\mathrm{new}}$ must be retained because this plane-parallel model
tracks photon position only along the vertical $z$ coordinate. The horizontal
$x$ and $y$ positions and direction components do not affect whether a photon
leaves through the top or bottom boundary.

#### Scattering Algorithm

For each interaction inside the layer, the model performs these steps:

1. Sample the exponential free path.
2. Move the photon using $\Delta z=\mu l$.
3. Check whether the photon exits through the top or bottom boundary.
4. If absorption is enabled, use $\omega_0$ to decide whether the interaction
   causes scattering or absorption.
5. If scattering occurs, sample the HG scattering angle $\Theta$.
6. Sample an independent random azimuth $\phi$.
7. Rotate the local scattering direction to obtain $\mu_{\mathrm{new}}$.
8. Continue the photon trajectory using the new vertical direction cosine.

### Single-Scattering Albedo

The single-scattering albedo $\omega_0$ determines the outcome of an
interaction inside the layer:

$$
P(\text{scattering})=\omega_0,
\qquad
P(\text{absorption})=1-\omega_0.
$$

An absorbed photon is terminated and added to the absorption count. Therefore,

$$
R+T+A=1.
$$

When $\omega_0=1$, every interaction scatters, absorption is zero, and the
nonabsorbing conservation relation becomes

$$
R+T=1.
$$

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
