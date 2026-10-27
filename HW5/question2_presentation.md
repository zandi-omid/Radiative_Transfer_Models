# HW5 Question 2: Presentation Notes

## Setup

- We calculate Mie scattering by homogeneous spheres with inner radius equal
  to outer radius, following the supplied program and assignment.
- Particle size is expressed as $x=2\pi r/\lambda$; each curve contains
  1,000 calculations from $x=0.05$ to 50 at increments of 0.05.
- Blue represents $m=1.3873-0.0066i$; orange represents $m=1.2993-0.037i$.
  Positive absorption magnitudes are passed to the Fortran routine, which
  constructs negative imaginary refractive indices internally. The assignment
  graph uses the alternative positive-imaginary convention.

## Slide 1: Extinction Efficiency

- **Extinction is scattering plus absorption:**
  $Q_{\mathrm{ext}}=Q_{\mathrm{scat}}+Q_{\mathrm{abs}}$.
- **Small particles interact weakly:** extinction tends toward zero as
  $x\to0$ for these fixed material indices.
- **The first broad peaks occur in the Mie transition regime:** the blue
  case peaks at $Q_{\mathrm{ext}}\approx4.017$
  near $x=5.45$; the orange case peaks at about 3.373 near $x=6.65$.
- **The oscillations are wave-interference effects:** the relative phases of
  the scattered contributions change as particle size changes.
- **The blue curve has more pronounced peaks and ripples:** its weaker
  absorption allows internal wave contributions to persist. Stronger
  absorption in the orange material damps that structure. The real indices
  also differ and change phase shifts and peak locations; this comparison
  does not isolate absorption alone.
- **Both curves tend toward $Q_{\mathrm{ext}}\approx2$ at large size:**
  geometrical interception and forward diffraction both contribute to removal
  from the incident beam. At $x=50$ they are still about 2.075 and 2.149,
  respectively, rather than exactly two.
- **An efficiency above one is allowed:** it is a cross-section divided by
  geometric area, not an energy fraction or a reflection probability.

## Slide 2: Asymmetry Factor

- **$g=\langle\cos\Theta\rangle$ describes scattering directionality:**
  $g=0$ means zero mean scattering-angle cosine, and positive $g$ indicates a
  forward bias. The angle is relative to the incident ray, not the ground.
- **At small size, $g$ approaches zero:** the Rayleigh-limit angular pattern
  has balanced forward and backward contributions; it is not isotropic.
- **Increasing size generally increases forward scattering:** the forward
  diffraction peak narrows, although interference causes dips and ripples
  rather than a strictly monotonic increase in $g$.
- **The orange case is more forward-biased at the larger sizes shown:** at
  $x=50$, $g\approx0.9721$, compared with about 0.9112 for the blue case.
- **Absorption can remove internally propagated, wider-angle contributions:**
  the remaining scattered radiation can then be more strongly weighted toward
  the forward diffraction peak. The different real indices also affect the
  angular pattern, so absorption is not the only changing property.
- **High $g$ does not mean high transmission:** $g$ averages only the angular
  distribution of scattered light. It neither measures absorbed energy nor
  gives the probability of passing through a cloud unscattered.

## Numerical Anchors

| Case | Peak $Q_{\mathrm{ext}}$ | $x$ at peak | $Q_{\mathrm{ext}}$ at $x=50$ | $g$ at $x=50$ |
|---|---:|---:|---:|---:|
| Blue: $1.3873-0.0066i$ | 4.017045 | 5.45 | 2.075354 | 0.911235 |
| Orange: $1.2993-0.037i$ | 3.372965 | 6.65 | 2.148664 | 0.972076 |

Peak locations are sampled-grid maxima, not continuously optimized values.

## Take-Home Message

- **Size controls the broad transition from weak interaction to forward-peaked
  scattering; material properties control the detailed extinction and angular
  structure.** Stronger absorption smooths interference features, but a
  forward-peaked phase function must not be confused with absence of absorption.

## Figure Caption

Mie extinction efficiency (a) and asymmetry factor (b) as functions of size
parameter $x=2\pi r/\lambda$ for homogeneous spherical particles in air.
Blue and orange curves use $m=1.3873-0.0066i$ and $m=1.2993-0.037i$,
respectively, in the supplied Fortran program's sign convention. A fixed
wavelength of 0.55 micrometers maps $x$ to physical radius; inner and outer
radii are equal. Calculations cover $0.05\leq x\leq50$ in increments of
0.05. The horizontal dotted line indicates $Q_{\mathrm{ext}}=2$.
The phase distribution is represented here by its mean cosine $g$, not a
forward-scattering probability. Curves are unsmoothed calculations from the
course-supplied `DMIESS` routine.

## Checks And Sources

- All 2,000 computed cases are finite and have nonnegative extinction,
  scattering, and absorption efficiencies and $g$ within $[-1,1]$.
- A separate direct spherical-Bessel series agrees within
  $2.841\times10^{-7}$ for extinction efficiency and
  $6.644\times10^{-8}$ for $g$. These checks concern single-sphere Mie
  calculations, not atmospheric-layer reflectivity or transmission.
- Course slides: 19-23 (angular scattering), 24-35 (extinction), and 36-38
  (single-scattering albedo), plus the supplied HW5 Question 2 graph and
  Fortran routine.
- [PyMieScatt's documented Mie coefficient sums](https://pymiescatt.readthedocs.io/en/latest/forward.html)
  and [SciPy spherical Bessel functions](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.spherical_jn.html)
  provide the equations and special-function API used for the independent
  reference calculation. PyMieScatt itself is not required or installed.
- [NASA: Interpretation of Extinction in Gaussian-Beam Scattering](https://ntrs.nasa.gov/citations/19970022873)
  discusses the large-particle extinction efficiency near two.
