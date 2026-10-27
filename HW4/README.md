# HW4: Two-Stream Transport With Surface Albedo

HW4 extends the HW3 Liou (1974) discrete-ordinate model to a Lambertian
surface beneath a homogeneous atmosphere. All reflections between the layer
and surface are included. The assignment supplied at
`/scratch/omidzandi/RTM/HW4/ATMO656A_2026_HW4.doc` was checked; implementation
and results live here, in the existing Git repository.

The cases use optical depths 1 and 4, single-scattering albedos 0.99999 and
0.8, solar cosine 0.9, asymmetry factor 0.75, and surface albedos from 0 to 1
in steps of 0.1. Plot curves use finer surface-albedo sampling.

## Question 4 Report

[Question 4 results](question4_results.md) gives the computed surface-albedo
sensitivity table, scientific interpretation, numerical checks, and captions.
The [three-panel PNG](figures/energy_partition_three_panels.png) and
[PDF figure](figures/energy_partition_three_panels.pdf) show reflectivity,
normalized surface absorption, and atmospheric absorption, using a Times-style
font. Run `python -m HW4.question4` to regenerate this figure directly from
the 44 saved CSV rows, preserving the individual figures. This version
connects the original albedo samples; `python -m HW4.hw4` instead plots a finer
sampling of the same solver. Panel (c) is inferred from the energy budget:
its closure is an arithmetic check, reported separately from the independent
reflecting-boundary ODE comparison. Atmospheric absorption strongly reduces
surface influence for the thicker cloud: at $\omega_0=0.8$, the reflectivity
change between black and perfectly reflecting surfaces is 0.0343 for
$\tau=4$, compared with 0.5148 for $\tau=1$.

## Equations And Boundary Conditions

Optical depth increases downward. The two Gaussian ordinates are
$\pm\mu_1$, with $\mu_1=1/\sqrt{3}$. With positive flux magnitudes $U$ and
$D$ for the upward and downward diffuse streams, the homogeneous equations are

$$
\frac{d}{d\tau}\begin{pmatrix}U\\D\end{pmatrix}
=\begin{pmatrix}a&-b\\b&-a\end{pmatrix}
\begin{pmatrix}U\\D\end{pmatrix},
$$

$$
a=\frac{1-\frac12(\bar\omega_0+\bar\omega_1\mu_1^2)}{\mu_1},
\qquad b=\frac{\bar\omega_0-\bar\omega_1\mu_1^2}{2\mu_1},
\qquad k=\sqrt{a^2-b^2},
$$

where $\bar\omega_0=\omega_0$ and $\bar\omega_1=3\omega_0g$ are the
first two Legendre moments. These follow from the Gaussian discretization of
Liou's transport equation and give the eigenvalue in his Eq. (11).
The direct source and its particular solution remain those of HW3,
Liou Eqs. (5)-(13).

For a diffuse unit upward flux injected at the lower boundary, impose
$D(0)=0$ and $U(\tau_L)=1$. Since the matrix above squares to $k^2I$, its
transfer matrix is $\cosh(k\tau_L)I+\sinh(k\tau_L)M/k$. Applying the two
boundary conditions gives

$$
R_{\mathrm{diff,bottom}}=
\frac{b\sinh(k\tau_L)}{k\cosh(k\tau_L)+a\sinh(k\tau_L)},
\qquad
T_{\mathrm{diff,up}}=
\frac{k}{k\cosh(k\tau_L)+a\sinh(k\tau_L)}.
$$

For $k=0$ the exact limits are $R_{\mathrm{diff}}=b\tau_L/(1+a\tau_L)$
and $T_{\mathrm{diff}}=1/(1+a\tau_L)$. The homogeneous layer is symmetric
under reversal of illumination, so its two diffuse reflection coefficients
are equal, as are its two diffuse transmission coefficients. This does not
make the diffuse coefficients equal to the direct-beam-illuminated ones.

`liou_reflect_transmit()` supplies $R_{\mathrm{dir}}$ and
$T_{\mathrm{dir}}$ for sunlight entering at $\mu_0$. The latter includes
both unscattered and downward-scattered sunlight. A Lambertian surface is
represented by $U(\tau_L)=A_s[D(\tau_L)+F_{\mathrm{beam}}(\tau_L)]$;
there is no incoming diffuse radiation at the top.

Define $\mathcal D=1-A_sR_{\mathrm{diff,bottom}}$. Summing the geometric
series in $A_sR_{\mathrm{diff,bottom}}$ gives

$$
R_{\mathrm{total}}=R_{\mathrm{dir}}+
\frac{A_sT_{\mathrm{dir}}T_{\mathrm{diff,up}}}{\mathcal D},
$$

$$
\frac{F_{\mathrm{down,surface}}}{F_{\mathrm{inc}}}
=\frac{T_{\mathrm{dir}}}{\mathcal D},\qquad
F_{\mathrm{up,surface}}=A_sF_{\mathrm{down,surface}},\qquad
F_{\mathrm{net,surface}}=(1-A_s)F_{\mathrm{down,surface}}.
$$

Atmospheric absorption is calculated by summing absorption on the initial
direct-illuminated pass and all subsequent diffuse passes:

$$
A_{\mathrm{atm}}=F_{\mathrm{inc}}(1-R_{\mathrm{dir}}-T_{\mathrm{dir}})
+F_{\mathrm{up,surface}}(1-R_{\mathrm{diff}}-T_{\mathrm{diff}}).
$$

## Flux Normalization

Liou's incident horizontal solar flux is $F_{\mathrm{inc}}=\pi F_0\mu_0$.
The assignment specifies $\pi F_0=1$, so $F_{\mathrm{inc}}=0.9$.
The CSV flux columns and printed table use these actual flux units;
reflectivity is dimensionless. The net-flux figure is normalized by
$F_{\mathrm{inc}}$, as labeled. Consequently, the CSV energy equation is

$$
0.9R_{\mathrm{total}}+F_{\mathrm{net,surface}}+A_{\mathrm{atm}}=0.9.
$$

Dividing by 0.9 gives the usual unit-incident-flux balance. The reusable
`surface_transport()` function defaults to unit horizontal incident flux;
the HW4 driver explicitly passes `incident_flux=0.9`.

## Run And Validate

From the repository root:

```bash
python -m HW4.hw4
python -m unittest discover -s HW4/tests -v
python -m unittest discover -s tests -v
```

Outputs are `HW4/hw4_results.csv` and PNG/PDF versions of
`HW4/figures/reflectivity_vs_albedo` and
`HW4/figures/surface_net_flux_vs_albedo`. PNG resolution is 300 DPI.
The additional `HW4/figures/energy_partition_three_panels.png` and `.pdf`
combine reflectivity, normalized net surface flux, and normalized atmospheric
absorption in one row with a shared legend. The CSV also includes explicit
normalized net-flux and absorption columns. The script checks all CSV and
plotted energy fractions and prints their maximum energy-budget residual.
Existing HW1-HW3 models and results are preserved.

The diffuse implementation scales away growing exponentials and uses
`expm1` near $k=0$. Exact conservation is separately tested at $\omega_0=1$;
the required $\omega_0=0.99999$ experiment retains its small real absorption.
The minimum reflection denominator in the four cases is approximately 0.536,
so all multiple-reflection sums converge comfortably. Individual downward
and upward surface fluxes can exceed incident flux because radiation can
cross the surface boundary repeatedly.

Tests check the black-surface HW3 limit, perfect-reflector boundary, energy
balance, flux normalization, exact conservative limit, and a finite sum of
200 reflections. An independent RK4 calculation solves the discrete-ordinate
equations directly with the reflecting boundary and integrates volume
absorption. Across all 11 albedos for both optical depths and albedos
1.0, 0.99999, and 0.8, maximum differences are below $1.4\times10^{-13}$
for reflectivity and $1.7\times10^{-14}$ for surface downward flux.

## Results And Interpretation

Values below use $\pi F_0=1$; $F_{\mathrm{net}}$ has flux units.

| Optical depth | Albedo | Surface albedo | TOA reflectivity | Net surface flux |
|---:|---:|---:|---:|---:|
| 1 | 0.99999 | 0 | 0.023892 | 0.878484 |
| 1 | 0.99999 | 0.5 | 0.464259 | 0.482146 |
| 1 | 0.99999 | 1 | 0.999965 | 0 |
| 4 | 0.99999 | 0 | 0.316089 | 0.615453 |
| 4 | 0.99999 | 0.5 | 0.554667 | 0.400704 |
| 4 | 0.99999 | 1 | 0.999837 | 0 |
| 1 | 0.8 | 0 | -0.000643 | 0.686999 |
| 1 | 0.8 | 0.5 | 0.241999 | 0.363177 |
| 1 | 0.8 | 1 | 0.514129 | 0 |
| 4 | 0.8 | 0 | 0.051517 | 0.187451 |
| 4 | 0.8 | 0.5 | 0.067077 | 0.102335 |
| 4 | 0.8 | 1 | 0.085785 | 0 |

A brighter surface increases top reflectivity while decreasing net energy
absorbed at the surface. It also returns more photons to the atmosphere,
increasing atmospheric absorption. An absorbing atmosphere therefore keeps
reflectivity below one even above a perfectly reflecting surface. At depth 4
and atmospheric albedo 0.8, most radiation is absorbed before completing the
trip to the surface and back; TOA reflectivity is weakly sensitive to surface
albedo. Nearly conservative layers allow much more of the surface-reflected
radiation to escape at the top.

For depth 1, atmospheric albedo 0.8, and a black surface, reflectivity is
approximately -0.000643. It is retained and flagged by the driver. The
independent reflecting-boundary solver reproduces it: the first-order angular
expansion has a negative upward direct source for this nearly vertical beam.
This is the same discrete-ordinate closure limitation documented in HW3.
No negative result is clipped. The uploaded homework references a class
derivation but does not supply it; the derivation here is checked against
Liou's equations and the independently solved boundary-value problem.
