# HW3 Question 1: Meador-Weaver Two-Stream Model

This experiment calculates layer reflectivity and transmissivity using the
Eddington and quadrature two-stream approximations from Meador and Weaver
(1980). The lower boundary is black, and the layer is homogeneous.

Unlike the HW1/HW2 Monte Carlo convention, the analytical equations use a
**positive** downward solar incidence cosine:

$$
0 < \mu_0 \leq 1.
$$

## Parameters

- Optical depth: $0 \leq \tau \leq 10$ at 201 evenly spaced values.
- Solar cosine: $\mu_0=0.9, 0.5, 0.1$.
- Single-scattering albedo: $\omega_0=1.0, 0.9$.
- Asymmetry factor: $g=0.75$.
- Closures: Eddington and quadrature.

## Two-Stream Coefficients

For the Eddington closure, Table 1 gives

$$
\gamma_1=\frac{7-\omega_0(4+3g)}{4}, \qquad
\gamma_2=-\frac{1-\omega_0(4-3g)}{4},
$$

$$
\gamma_3=\frac{2-3g\mu_0}{4}.
$$

For the quadrature closure,

$$
\gamma_1=\frac{\sqrt{3}}{2}\left[2-\omega_0(1+g)\right], \qquad
\gamma_2=\frac{\sqrt{3}}{2}\omega_0(1-g),
$$

$$
\gamma_3=\frac{1-\sqrt{3}g\mu_0}{2}.
$$

Both methods then use Eqs. (16)-(18):

$$
\gamma_4=1-\gamma_3,
$$

$$
\alpha_1=\gamma_1\gamma_4+\gamma_2\gamma_3, \qquad
\alpha_2=\gamma_1\gamma_3+\gamma_2\gamma_4,
$$

$$
k=\sqrt{\gamma_1^2-\gamma_2^2}.
$$

## Stable Reflectivity and Transmissivity

For $\omega_0<1$, `models/two_stream.py` evaluates exponent-scaled forms that
are algebraically identical to Meador and Weaver Eqs. (14)-(15). Define

$$
e=\exp(-k\tau), \qquad b=\exp(-\tau/\mu_0),
$$

$$
D_s=(k+\gamma_1)+(k-\gamma_1)e^2.
$$

The common denominator is

$$
\left[1-(k\mu_0)^2\right]D_s.
$$

The scaled implementation is

$$
R=\frac{\omega_0}{[1-(k\mu_0)^2]D_s}\left[(1-k\mu_0)(\alpha_2+k\gamma_3)-(1+k\mu_0)(\alpha_2-k\gamma_3)e^2-2k(\gamma_3-\alpha_2\mu_0)be\right],
$$

$$
T=b-\frac{\omega_0}{[1-(k\mu_0)^2]D_s}\left\{b\left[(1+k\mu_0)(\alpha_1+k\gamma_4)-(1-k\mu_0)(\alpha_1-k\gamma_4)e^2\right]-2k(\gamma_4+\alpha_1\mu_0)e\right\}.
$$

This scaling avoids the large intermediate exponentials in the direct forms,
especially for small $\mu_0$ and large optical depth.

### Conservative Limit

When $\omega_0=1$, $k=0$, so direct substitution into Eqs. (14)-(15) is
singular. The code uses the analytical $k\rightarrow0$ limit (Eq. 24):

$$
R=\frac{\gamma_1\tau+(\gamma_3-\gamma_1\mu_0)
[1-\exp(-\tau/\mu_0)]}{1+\gamma_1\tau},
$$

$$
T=1-R.
$$

No arbitrary epsilon or artificial absorption is introduced. At $\tau=0$,
both branches give $R=0$ and $T=1$.

For $\omega_0=0.9$, the absorbed fraction is diagnosed as

$$
A=1-R-T.
$$

## Run

From the repository root:

```bash
python -m HW3.hw3_question1
```

The script uses a noninteractive Matplotlib backend and writes four 400-DPI
figures and one CSV file to `HW3/results/`.

Run all repository tests with:

```bash
python -m unittest discover -s tests -v
```

The two-stream tests verify the zero-depth boundary, conservative energy
balance, absorbing energy partition, finite output over all requested cases,
and agreement between the stable implementation and the direct analytical
expressions where they are nonsingular.

## Interpretation and Limitation

Increasing optical depth generally increases reflectivity and decreases
transmissivity. The absorbing case also transfers energy into $A$. Eddington
and quadrature agree closely for many conditions but are still low-order
angular approximations.

For $g=0.75$ and nearly vertical illumination ($\mu_0=0.9$), both closures can
produce slightly negative reflectivity in optically thin layers. The largest
case in this sweep is approximately $-0.0134$ for quadrature. These values are
reported and plotted without clipping because they are a limitation of the
closure, rather than a numerical overflow or nonfinite result.
