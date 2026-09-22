# HW1 Monte Carlo Photon Transport

Four homework analyses in `HW1/` share the photon transport model in
`models/monte_carlo.py`. The assignment is in `HW1/HW1.docx`.
The layer is plane-parallel, nonabsorbing, and isotropically scattering.
Each photon exits through either the top (reflectivity) or bottom
(transmissivity). Multiprocessing uses up to 24 workers and independent,
reproducible random streams for each worker.

## Setup

```bash
cd Radiative_Transfer_Models
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On a cluster or remote machine without a display, set `MPLBACKEND=Agg`.

## Run

```bash
MPLBACKEND=Agg python -m HW1.hw1_monte_carlo_photons
MPLBACKEND=Agg python -m HW1.hw1_problem2_independent_runs
MPLBACKEND=Agg python -m HW1.hw1_problem3_incident_direction
MPLBACKEND=Agg python -m HW1.hw1_problem4_optical_depth
```

Problem 1 studies convergence with photon count. Problem 2 compares ten
independent runs. Problem 3 varies the incident direction. Problem 4 varies
optical depth and reports mean scatterings and runtime. Tables and figures
are saved beside the scripts in `HW1/` and ignored by Git.

## Test

```bash
python -m unittest discover -s tests -v
```
