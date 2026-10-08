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

## Independent Validation

The equations and coefficients were checked against the original Meador and
Weaver (1980) paper. The signs and exponential factors in Eqs. (14)-(15), the
definitions of $\alpha_1$, $\alpha_2$, and $k$ in Eqs. (16)-(18), the Table 1
Eddington and quadrature coefficients, and $\gamma_4=1-\gamma_3$ all match the
implementation.

The conservative expression can also be obtained without assuming energy
conservation. Normalize the incident direct flux to one and define
$D=I^+-I^-$. When $\omega_0=1$, $\gamma_1=\gamma_2$ and
$\gamma_3+\gamma_4=1$, so subtracting the two coupled flux equations gives

$$
\frac{dD}{d\tau}=-\frac{1}{\mu_0}\exp(-\tau/\mu_0),
\qquad
D(\tau)=R-[1-\exp(-\tau/\mu_0)].
$$

Substitution into the equation for $I^+$, followed by integration from the
top to the black lower boundary, gives

$$
R(1+\gamma_1\tau)=\gamma_1\tau+
(\gamma_3-\gamma_1\mu_0)[1-\exp(-\tau/\mu_0)],
$$

which is Eq. (24) and the limit used by the code. Thus $T=1-R$ follows for the
conservative black-boundary problem.

`HW3/validate_question1.py` supplies a separate numerical check. It integrates
the coupled two-stream differential equations (Eqs. 12-13) with a fourth-order
Runge-Kutta shooting method, rather than evaluating Eqs. (14)-(15). Run it with

```bash
python -m HW3.validate_question1
```

The comparison covers all 12 requested combinations over $0\leq\tau\leq10$
and writes maximum absolute errors to
`HW3/results/hw3_q1_validation_errors.csv`. The reported errors use a maximum
RK4 step of 0.001. The largest absolute differences are
$3.238\times10^{-12}$ for $R$ and $1.580\times10^{-12}$ for $T$. Halving the
step from 0.002 changes the reference solution by less than
$4.4\times10^{-11}$ for $R$ and $2.4\times10^{-11}$ for $T$.

The negative thin-layer reflectivity is reproduced by this independent ODE
solution. It is also anticipated by the paper's thin-layer analysis: the
standard Eddington closure becomes negative when $g>2/(3\mu_0)$, and the paper
notes a corresponding thin-layer limitation for quadrature. For
$g=0.75$ and $\mu_0=0.9$, that condition is satisfied. In a conservative
solution $T=1-R$, so a negative approximate $R$ necessarily produces $T>1$.
These values are mathematical artifacts of the low-order angular closures,
not clipping, floating-point cancellation, or a programming error.
