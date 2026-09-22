# HW1 Monte Carlo Photon Transport

Four homework analyses share the photon transport model in `monte_carlo.py`.
The layer is plane-parallel, nonabsorbing, and isotropically scattering.
Each photon exits through either the top (reflectivity) or bottom
(transmissivity). Multiprocessing uses up to 24 workers and independent,
reproducible random streams for each worker.

## Setup

```bash
cd HW1
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On a cluster or remote machine without a display, set `MPLBACKEND=Agg`.

## Run

```bash
MPLBACKEND=Agg python hw1_monte_carlo_photons.py
MPLBACKEND=Agg python hw1_problem2_independent_runs.py
MPLBACKEND=Agg python hw1_problem3_incident_direction.py
MPLBACKEND=Agg python hw1_problem4_optical_depth.py
```

Problem 1 studies convergence with photon count. Problem 2 compares ten
independent runs. Problem 3 varies the incident direction. Problem 4 varies
optical depth and reports mean scatterings and runtime. Tables and figures
are saved beside the scripts and ignored by Git.

## Test

```bash
python -m unittest discover -s tests -v
```
