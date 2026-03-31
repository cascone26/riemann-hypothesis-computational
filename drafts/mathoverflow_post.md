# MathOverflow Draft — Heat Kernel Dimension of Riemann Zeta Zeros

## Title
Scale-dependent effective dimension from the heat trace of Riemann zeta zeros

## Tags
riemann-zeta-function, spectral-theory, random-matrices, nt.number-theory

## Body

### Background

The Hilbert-Polya approach to the Riemann Hypothesis seeks a self-adjoint operator $H$ whose eigenvalues are the imaginary parts $\gamma_n$ of the non-trivial zeros. The heat trace of such an operator,
$$\Theta(t) = \sum_n e^{-t\gamma_n},$$
encodes spectral information about $H$, and its short-time asymptotics reveal the "dimension" of the underlying space via the Weyl law.

### Computation

Using the first 2,000,000 zeros from Odlyzko's tables, I computed $\Theta(t)$ for $t \in [10^{-4}, 10]$ and fitted the short-time behavior $\Theta(t) \sim A \cdot t^{\alpha}$.

**Results:**
- At $t \in [10^{-4}, 10^{-3}]$: effective dimension $d_{\mathrm{eff}} = -2\alpha \approx 2.36$
- At $t \in [10^{-3}, 10^{-2}]$: $d_{\mathrm{eff}} \approx 2.60$
- At $t \in [10^{-2}, 10^{-1}]$: $d_{\mathrm{eff}} \approx 3.52$

The effective dimension is **scale-dependent**, increasing from $\approx 2.36$ at short time scales to $\approx 3.52$ at larger scales. This is consistent with the counting function $N(E) \sim \frac{E}{2\pi}\log\frac{E}{2\pi}$ growing faster than any polynomial $E^{d/2}$.

The result is stable: 100K and 2M zeros give identical values ($d = 2.4630$ vs $d = 2.4631$).

### Question

1. Is this scale-dependent effective dimension known in the literature? I'm aware of the connection to the Weyl law and the log correction, but I haven't found explicit numerical computation of the running dimension from the heat trace of zeta zeros.

2. The non-integer, scale-dependent dimension is consistent with operators on noncommutative spaces (Connes' program) or fractal geometries (Lapidus' fractal membranes). Does this computation contribute any constraint beyond what the Weyl law already tells us, or is it a straightforward repackaging?

3. The Wu-Sprung potential reconstruction gives a fractal potential with box-counting dimension $\approx 1.5$. Is there a known relationship between the Wu-Sprung fractal dimension and the heat kernel effective dimension computed here?

### Code

All computations use Python with mpmath, numpy, and scipy. The heat trace is computed by direct summation $\sum_n e^{-t\gamma_n}$ using Odlyzko's tables (9 decimal places). The effective dimension at scale $t$ is computed as the local log-log slope: $d_{\mathrm{eff}}(t) = -2 \cdot \frac{d \log \Theta}{d \log t}$.

---

*Note: I've also formalized Robin's inequality statement in Lean 4 using Mathlib's existing infrastructure (`ArithmeticFunction.sigma`, `Real.eulerMascheroniConstant`, `RiemannHypothesis`). The formalization compiles against current Mathlib (v4.29.0). Happy to share details if there's interest in a Mathlib contribution.*
