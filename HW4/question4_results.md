# HW4 Question 4: Results And Physical Interpretation

Question 4 asks how a reflecting or absorbing surface changes the radiation
leaving an atmospheric layer. These results use the existing Liou two-stream
discrete-ordinate model with all cloud-surface reflections. No transport
equations or physical assumptions were changed for this analysis.

## Parameters And Normalization

The optical depths are $\tau=1,4$, single-scattering albedos are
$\omega_0=0.99999,0.8$, solar cosine is $\mu_0=0.9$, asymmetry factor is
$g=0.75$, and Lambertian surface albedo is $0\leq A_s\leq1$. The assignment
sets $\pi F_0=1$, giving incident horizontal flux
$F_{\mathrm{incident}}=\pi F_0\mu_0=0.9$.

The figure and tables use fractions of incident horizontal solar flux:

$$
S=\frac{F_{\mathrm{net},s}}{F_{\mathrm{incident}}},\qquad
A_{\mathrm{atm}}=1-R_{\mathrm{TOA}}-S.
$$

Here $S$ is energy absorbed by the surface and $A_{\mathrm{atm}}$ is the
fraction absorbed by the atmosphere. Raw flux columns in `hw4_results.csv`
have the original flux units and must be divided by 0.9 before comparing with
these normalized quantities. The figure connects the 11 saved albedo values
for each case; no new model results are substituted for the existing CSV.

## Surface Influence On Reflectivity

The following values are calculated from `hw4_results.csv`, not estimated
from the plot. $\Delta R=R(A_s=1)-R(A_s=0)$ is the change in reflected
fraction when the surface changes from perfectly absorbing to perfectly
reflecting.

| $\tau$ | $\omega_0$ | $R(A_s=0)$ | $R(A_s=0.5)$ | $R(A_s=1)$ | $\Delta R$ |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.99999 | 0.023892451 | 0.464259186 | 0.999965273 | 0.976072822 |
| 4 | 0.99999 | 0.316088734 | 0.554667491 | 0.999836868 | 0.683748134 |
| 1 | 0.8 | -0.000643009 | 0.241998548 | 0.514129277 | 0.514772285 |
| 4 | 0.8 | 0.051516798 | 0.067077151 | 0.085785439 | 0.034268641 |

For comparison, the atmospheric absorption fractions are:

| $\tau$ | $\omega_0$ | $A_{\mathrm{atm}}(A_s=0)$ | $A_{\mathrm{atm}}(A_s=1)$ |
|---:|---:|---:|---:|
| 1 | 0.99999 | 0.000014161 | 0.000034727 |
| 4 | 0.99999 | 0.000074733 | 0.000163132 |
| 1 | 0.8 | 0.237311340 | 0.485870723 |
| 4 | 0.8 | 0.740204628 | 0.914214561 |

## Physical Interpretation

### Increasing Surface Albedo

A brighter surface returns a greater fraction of the downward radiation to
the atmosphere. Some of this radiation passes through the cloud and escapes
at the top, increasing TOA reflectivity. The surface absorbs a smaller fraction
of the radiation incident on it, so its net downward flux decreases. At
$A_s=1$, downward and upward surface fluxes are equal and net surface
absorption is exactly zero, although radiation still reaches the surface.

### Nearly Nonabsorbing Clouds

For $\omega_0=0.99999$, almost every interaction scatters rather than absorbs
radiation. Multiple cloud-surface reflections therefore redistribute energy
with little atmospheric loss. With a perfectly reflecting surface, almost
all incoming energy eventually escapes upward: TOA reflectivity reaches
0.999965 for $\tau=1$ and 0.999837 for $\tau=4$. The remaining atmospheric
absorption is small but nonzero because 0.99999 is not exactly one. Across
these cases it never exceeds about 0.0164% of the incident energy.

### Absorbing Clouds

For $\omega_0=0.8$, each additional passage through the layer provides another
opportunity for absorption. A bright surface sends radiation back into the
cloud, increasing both upward escape and atmospheric absorption. For
$\tau=1$, atmospheric absorption rises from about 23.7% above a black surface
to 48.6% above a perfectly reflecting surface. For $\tau=4$, it rises from
74.0% to 91.4%. Thus a perfectly reflecting surface does not make the whole
atmosphere-surface system perfectly reflecting when the atmosphere absorbs.

### Optical Depth And Surface Visibility

The thicker nearly conservative cloud already reflects more radiation above a
black surface, leaving a smaller change for the surface to produce. Its
reflectivity still responds strongly to surface albedo. In the absorbing
case, increasing optical depth also increases absorption along the paths
to and from the surface. At $\tau=4$, little surface-reflected energy reaches
space: $\Delta R=0.0343$, compared with 0.5148 at $\tau=1$. This change is
about fifteen times smaller for the thicker absorbing cloud.

### Absorbing And Reflecting Surfaces

A black surface ($A_s=0$) absorbs all radiation that reaches it, with no
surface-reflected contribution to TOA reflectivity. A perfectly reflecting
surface ($A_s=1$) absorbs none; energy is ultimately either reflected to space
or absorbed in the atmosphere. A nearly conservative atmosphere sends almost
all that energy to space, whereas the thick absorbing cloud absorbs most of
it. Atmospheric absorption and optical thickness therefore determine how
strongly the surface influences reflectivity observed from above.

## Checks And Numerical Limitation

All 44 saved cases (four atmospheric cases, 11 surface albedos) were checked.
The largest absolute residual in $R_{\mathrm{TOA}}+S+A_{\mathrm{atm}}-1$ is
**$1.110223025\times10^{-16}$**. Since the plotted absorption is inferred
from the other two fractions, this is an arithmetic closure check and does
not independently validate the solver.

The residual-inferred absorption also agrees with the CSV's separately
stored sum of absorption over direct and repeated diffuse passes, to
$2.220446049\times10^{-16}$. Both calculations share layer coefficients;
their agreement alone is not an independent physical validation.

A separate RK4 integration of the discrete-ordinate differential equations
was performed for all 44 cases, using a reflecting lower boundary and a
volume absorption integral. It uses neither the layer adding formula nor
absorption inferred from energy balance. Maximum absolute differences from
the saved normalized results were:

| Quantity | Maximum absolute difference |
|---|---:|
| TOA reflectivity | $1.703\times10^{-14}$ |
| Downward surface flux fraction | $1.155\times10^{-14}$ |
| Net surface flux fraction | $8.660\times10^{-15}$ |
| Atmospheric absorption fraction | $1.224\times10^{-14}$ |

Net surface flux is zero at $A_s=1$ in every case. At $A_s=0$, reflectivity
and normalized net surface flux agree with the HW3 black-boundary Liou
solution within $10^{-13}$. Nearly conservative absorption remains below
0.001, while the absorbing cases all exceed 0.2. The existing seven HW4
tests, including a 200-term reflection sum and the exact conservative limit,
also passed as part of the 31-test repository suite.

For $\tau=1$, $\omega_0=0.8$, and $A_s=0$, TOA reflectivity is
**-0.000643009**. This is unphysical and is neither clipped nor interpreted
as real negative reflected energy. The independent ODE solution reproduces
it. At this large $\mu_0$ and $g$, the first-order angular approximation gives
a negative upward direct-beam scattering source. This is a limitation of the
two-stream approximation, rather than a failure to sum multiple reflections
or a roundoff artifact. The corresponding $\Delta R$ includes this small
unphysical endpoint and should be read with that qualification. Numerical
agreement validates the solution of the approximate equations, not their
accuracy relative to a full angular radiative transfer solution.

## Figure Captions

**Three-panel figure.** Energy partition versus Lambertian surface albedo
$A_s$ from 0 to 1, calculated with Liou's two-stream discrete-ordinate model
including all cloud-surface reflections. Panel (a) shows top-of-atmosphere
reflectivity $R_{\mathrm{TOA}}$; panel (b) shows net downward surface flux
normalized by incident horizontal solar flux,
$F_{\mathrm{net},s}/F_{\mathrm{incident}}$; panel (c) shows atmospheric
absorption inferred from the energy budget,
$A_{\mathrm{atm}}=1-R_{\mathrm{TOA}}-F_{\mathrm{net},s}/F_{\mathrm{incident}}$.
Blue and orange denote optical depths $\tau=1$ and 4. Solid and dashed lines
denote single-scattering albedos $\omega_0=0.99999$ and 0.8. All cases use
$\mu_0=0.9$, $g=0.75$, and $\pi F_0=1$, so
$F_{\mathrm{incident}}=0.9$. Curves connect the saved albedo samples at
increments of 0.1. The common vertical range includes a small margin below
zero to retain the thin absorbing layer's unphysical negative reflectivity.
Liberation Serif is used because Times New Roman is unavailable on Rain.

**Individual reflectivity figure.** TOA reflected fraction versus surface
albedo for Liou two-stream layers with $\tau=1$ (blue), 4 (orange),
$\omega_0=0.99999$ (solid), and 0.8 (dashed). All cases have $\mu_0=0.9$,
$g=0.75$, and $\pi F_0=1$; all cloud-surface reflections are included.

**Individual net-flux figure.** Net downward surface flux divided by the
incident horizontal solar flux versus surface albedo for the same four
cases. Blue/orange denote $\tau=1/4$ and solid/dashed denote
$\omega_0=0.99999/0.8$, with $\mu_0=0.9$, $g=0.75$, and $\pi F_0=1$.
The net flux is the energy absorbed by the surface and vanishes at $A_s=1$.

## Outputs And Reproduction

- [Three-panel PNG, 300 DPI](figures/energy_partition_three_panels.png)
- [Three-panel PDF figure](figures/energy_partition_three_panels.pdf)
- [Preserved individual reflectivity figure](figures/reflectivity_vs_albedo.png)
- [Preserved individual net-flux figure](figures/surface_net_flux_vs_albedo.png)
- [Existing numerical CSV](hw4_results.csv)

From the repository root, run `python -m HW4.question4` to check the saved
results, print the comparison table, perform the independent numerical
comparison, and regenerate only the three-panel figure. Run
`python -m unittest discover -s tests -v` for the complete test suite.
The PNG/PDF outputs are figures only; no combined Questions 1-4 homework PDF
has been produced.
