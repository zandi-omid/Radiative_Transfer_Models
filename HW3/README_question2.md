# HW3 Question 2: Liou Two-Stream Discrete Ordinates

Question 2 implements the analytic two-stream solution in Liou (1974),
Eqs. (5)-(13), and compares it with the Meador-Weaver quadrature closure used
in Question 1. The atmosphere is a homogeneous layer over a black lower
boundary, with no incoming diffuse radiation at either boundary. Optical depth
increases downward, and the solar cosine $\mu_0$ is positive.

## Discrete-Ordinate Model

The two-stream Gaussian quadrature uses one direction in each hemisphere,

$$
\mu_1=\frac{1}{\sqrt{3}}.
$$

For a first-order Henyey-Greenstein phase-function expansion, the two Legendre
moments in Liou's notation are

$$
\bar{\omega}_0=\omega_0, \qquad
\bar{\omega}_1=3\omega_0g.
$$

Liou Eq. (11) then gives the eigenvalue

$$
k_1=\frac{1}{\mu_1}
\sqrt{(1-\bar{\omega}_0)(1-\bar{\omega}_1\mu_1^2)}.
$$

The eigenfunction and direct-beam particular solution from Eqs. (12)-(13),
written for either discrete direction $X=\pm\mu_1$, are

$$
W_1(X)=\frac{\bar{\omega}_0-
\bar{\omega}_1(1-\bar{\omega}_0)X/k_1}{1+Xk_1},
$$

$$
Z(X)=\frac{\mu_0F_0(X-\mu_0)
[\bar{\omega}_0-\bar{\omega}_1(1-\bar{\omega}_0)X\mu_0]}
{4\mu_1^2(1-k_1^2\mu_0^2)}.
$$

`models/liou_two_stream.py` evaluates the boundary coefficients $A^\pm$,
$C^\pm$, $L_1$, and $L_{-1}$ from Eqs. (7)-(10), then evaluates the upward and
downward fluxes from Eqs. (5)-(6). Reflectivity is the upward flux at the top
divided by the incident flux $\pi\mu_0F_0$.

Total transmissivity includes both terms at the lower boundary:

$$
T=T_{\mathrm{diffuse}}+\exp(-\tau/\mu_0).
$$

The second term is the direct solar beam. Omitting it would make the
zero-optical-depth result incorrect.

## Conservative Limit

For $\omega_0=1$, Eq. (11) gives $k_1=0$, making Eqs. (12)-(13) singular when
used term by term. Liou suggested using $\omega_0=0.99999$ on the computers
available in 1974. This implementation instead uses the exact $k_1\to0$ limit
of the coupled discrete-ordinate equations:

$$
R=\frac{\gamma_1\tau+(\gamma_3-\gamma_1\mu_0)
[1-\exp(-\tau/\mu_0)]}{1+\gamma_1\tau},
$$

where

$$
\gamma_1=\frac{\sqrt{3}}{2}(1-g), \qquad
\gamma_3=\frac{1-\sqrt{3}g\mu_0}{2},
$$

and $T=1-R$. No artificial absorption or arbitrary epsilon is introduced.

## Experiment And Figures

The comparison uses $0\leq\tau\leq10$, $\mu_0=0.9,0.5,0.1$,
$\omega_0=1.0,0.9$, and $g=0.75$. Run it from the repository root:

```bash
python -m HW3.hw3_question2
```

Four 400-DPI figures and `hw3_question2_results.csv` are written to
`HW3/results/`. Colors identify $\mu_0$; solid lines show Meador-Weaver
quadrature and dashed lines show Liou.

The curves overlap because Liou's two-point Gaussian discrete ordinates and
the Meador-Weaver quadrature closure are algebraically equivalent. Across the
requested sweep, the largest computed differences are below $7\times10^{-16}$.

For $g=0.75$ and $\mu_0=0.9$, both formulations produce slightly negative
reflectivity in optically thin layers, reaching about $-0.0134$. In the
conservative case this also makes $T$ slightly greater than one. The values are
retained rather than clipped. Their identical appearance in both analytic
forms and in an independent RK4 solution shows that they are a limitation of
the low-order angular approximation, not an implementation or cancellation
error.
