# HW5 Question 4: Presentation Notes

## Copy-And-Paste Slide Notes

- The figure shows how strongly light is scattered into each direction for
  particles of radius 4 micrometers at wavelength 0.55 micrometers.
- Blue represents nonabsorbing particles ($m=1.5$); red represents absorbing
  particles ($m=2-0.66i$).
- Both curves have a strong forward peak near 0 degrees, showing that
  scattering is concentrated near the incoming light direction.
- The blue curve has pronounced oscillations caused by interference between
  scattered-wave contributions and rises again near the backward direction.
- Absorption damps internally propagated contributions in the red case,
  producing a smoother pattern at larger angles. The different real
  refractive index also changes the angular structure.
- The red case is more forward-directed overall: $g\approx0.865$, compared
  with $g\approx0.781$ for blue, consistent with Q3 at the same radius.
- The logarithmic vertical axis reveals weak side- and backward scattering
  alongside the much stronger forward peak.
- The phase function describes the directions of scattered light, not the
  fraction of incident energy transmitted or absorbed by the particle.

## Key Values

| Case | $P(0\,\mathrm{deg})$ | $P(90\,\mathrm{deg})$ | $P(180\,\mathrm{deg})$ | $g$ |
|---|---:|---:|---:|---:|
| Blue: $m=1.5$ | 1066.028 | 0.245471 | 0.165521 | 0.780943 |
| Red: $m=2-0.66i$ | 1932.451 | 0.130073 | 0.120389 | 0.865109 |

The blue curve's backward-angle peaks occur near, but not necessarily
exactly at, 180 degrees; the value at the endpoint alone does not represent
the full backward angular structure.

## Figure Caption

Unpolarized Mie phase functions versus scattering angle for homogeneous
spheres of radius 4 micrometers at wavelength 0.55 micrometers
($x=45.6959$). Blue and red curves use $m=1.5-0i$ and $m=2-0.66i$,
respectively, in the supplied Fortran sign convention. The 181 samples span
0-180 degrees in 1-degree increments and are connected without smoothing.
Scattering angle is measured relative to the incoming ray: zero is forward
and 180 degrees is backward. The ordinate is logarithmic. Phase functions
use the supplied driver's normalization $(4\pi)^{-1}\int P\,d\Omega=1$;
they are angular distributions conditional on scattering, not incident-energy
fractions. No particle-density assumption is required.

## Normalization And Validation Notes

- $P$ is normalized to $4\pi$, so its average over solid angle is one.
  Values greater than one are allowed and indicate angular concentration.
- A separate CSV column gives the unit-integral distribution $P/(4\pi)$ in
  inverse steradians. Different normalization conventions should not be mixed
  when comparing phase-function plots.
- All 362 plotted values agree with an independent spherical-Bessel and
  Legendre-polynomial calculation within $2.1\times10^{-11}$ absolute error.
- Separate Gauss-Legendre integration verifies normalization and recovers
  the coefficient-derived $g$ within $9\times10^{-12}$.
- The requested 1-degree grid is too coarse to integrate the narrow forward
  peak accurately: trapezoidal integration gives 0.9862 for blue and 0.9750
  for red. The plotted values were not rescaled to force these sums to one.
- A single-precision pi literal in the legacy angular conversion was replaced
  with double precision; the original reference source remains preserved.

For background on phase-function conventions, see the author's
[Mie angular-scattering documentation](https://miepython.readthedocs.io/en/latest/03_angular_scattering.html).
Our plotted curves come from the supplied Fortran routine, not that package.

## Files

- [Figure PNG](figures/question4_phase_function.png) and
  [PDF](figures/question4_phase_function.pdf).
- [Angular data CSV](results/question4_phase_function.csv).
- [Normalization and moment audit](results/question4_moments.csv).
- [Independent comparison](results/question4_validation.csv).
