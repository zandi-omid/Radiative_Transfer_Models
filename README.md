# Radiative Transfer Models

This repository contains models and experiments for studying radiative transfer.
The current Monte Carlo model follows individual photons through a
plane-parallel, nonabsorbing layer with isotropic scattering. It estimates
reflectivity and transmissivity from the fraction of photons leaving through
the top and bottom of the layer. Independent random streams allow the
simulations to run in parallel and remain reproducible.

The shared transport model is in `models/`. The `HW1/` folder contains the
first homework assignment and analyses of photon-count convergence,
run-to-run variability, incident direction, and optical depth.
