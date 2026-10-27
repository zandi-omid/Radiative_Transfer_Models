# HW5 Question 3: Presentation Notes

## Plot Setup

- Two homogeneous-sphere cases: **red, absorbing** ($m=2-0.66i$), and
  **blue, nonabsorbing** ($m=1.5-0i$).
- Radius ranges from $10^{-3}$ to $10$ micrometers at wavelength 0.55 micrometers.
- Specific extinction and absorption are cross-sections per particle mass,
  measured in m^2/g; they differ from the dimensionless efficiencies in Q2.
- Density is explicitly assumed to be **1 g/cm^3 for both cases** because
  the assignment does not specify density. A different density rescales the
  specific cross-sections inversely but does not change $g$.

## Slide 1: Specific Extinction And Absorption

- **The absorbing material has finite absorption per gram even at very small
  radii.** Its small-particle specific absorption approaches about
  7.156 m^2/g under the stated wavelength and density assumptions.
- **The small-particle absorption curve is nearly flat:** absorption
  cross-section and particle mass both scale approximately with $r^3$.
- **The nonabsorbing material's extinction rises rapidly at small radius:**
  Rayleigh scattering cross-section scales with $r^6$, so scattering per unit
  mass scales with $r^3$.
- **The strongest extinction per gram occurs at intermediate sizes:** about
  19.41 m^2/g at 0.102 micrometers for the absorbing case and 10.01 m^2/g at
  0.310 micrometers for the nonabsorbing case, on the sampled grid.
- **Specific extinction decreases at large radius:** cross-section grows
  roughly with area ($r^2$), but mass grows with volume ($r^3$), giving an
  approximate $1/r$ dependence.
- **The nonabsorbing case has zero physical absorption.** Zero cannot appear
  on a logarithmic axis; its absorption is identified in the legend instead.
- **Absorbing particles do not always have larger extinction per gram:**
  extinction combines scattering and absorption, and the real refractive
  indices also differ between these two materials.

## Slide 2: Asymmetry Factor

- **Small particles have $g$ close to zero:** forward and backward scattering
  balance in the Rayleigh limit; the angular pattern is not isotropic.
- **Larger particles generally scatter more strongly forward**, although
  interference produces dips and oscillations.
- **The nonabsorbing material has pronounced angular structure.** Internal
  wave contributions remain important and create the oscillatory blue curve.
- **The absorbing material has a smoother $g$ curve.** Internal absorption
  damps these contributions; its different real index also affects the pattern.
- **At radius 10 micrometers, $g\approx0.866$ for the absorbing material and
  $0.809$ for the nonabsorbing material.** These are directionality measures,
  not fractions of incident light that are transmitted.

## Main Conclusion

- **Q2 shows efficiency per geometric area; Q3 shows effectiveness per mass.**
  Large particles can have $Q_{\mathrm{ext}}$ near two while their extinction
  per gram decreases. Size, refractive index, and the chosen mass density
  determine the magnitude and shape of the mass-specific curves.

## Figure Caption

Specific extinction and absorption (a) and scattering asymmetry factor (b)
versus radius for homogeneous spheres with $m=2-0.66i$ (red) and
$m=1.5-0i$ (blue), using the supplied Fortran sign convention. Solid and
dashed curves in panel (a) denote extinction and absorption. Calculations use
wavelength 0.55 micrometers, equal inner and outer radii, and an assumed
particle density of 1 g/cm^3 for both materials. The radius axis is logarithmic
over 0.001-10 micrometers. Panel (a) is also logarithmic vertically over
0.1-100 m^2/g to match the assignment. The nonabsorbing absorption curve is
zero within roundoff and is explicitly labeled rather than placed at an
artificial positive value. Its tiny signed numerical residuals and values
outside the displayed vertical range remain in the CSV. Curves use 5,451
radii per material without smoothing.

## Checks And Limitations

- All 10,902 cases were compared with an independent spherical-Bessel Mie
  calculation, including an explicit conversion from cross-section to mass.
- Maximum efficiency disagreement is below $1.5\times10^{-11}$; maximum
  specific-extinction disagreement is below $1.5\times10^{-12}$ m^2/g.
- Nonabsorbing $Q_{\mathrm{abs}}$ is zero within $1.34\times10^{-15}$;
  small signed residuals are retained without clipping.
- An exact homogeneous-sphere branch and a series convergence guard were
  needed to avoid two limitations in the supplied coated-sphere routine.
  The original supplied file is preserved; the fixes and regression checks
  are documented in [the HW5 README](README.md).
- Density is an explicit comparison assumption, not a value provided in the
  assignment. Absolute specific cross-sections should be rescaled if a
  different density is prescribed in class.
- The small-particle scaling follows the Rayleigh formulas documented in
  [PyMieScatt's reference equations](https://pymiescatt.readthedocs.io/en/latest/forward.html).

## Files

- [Combined figure, PNG](figures/question3_mie_radius.png) and
  [PDF](figures/question3_mie_radius.pdf).
- [Specific-cross-section figure](figures/question3_specific_cross_sections.png).
- [Asymmetry figure](figures/question3_asymmetry.png).
- [Full numerical results](results/question3_mie_results.csv).
- [Validation summary](results/question3_validation.csv).
