# Principal SOS Vector Discovery: Dominated by φ_0

**Date**: 2026-08-08  
**Finding**: The principal SOS vector l_0 (associated with λ_min) is **almost entirely concentrated in the first Fejer basis element**.

---

## The Observation

From principal_sos_analysis.py runs across 9 cases:

| Window | 4W | N | ||l_0|| | L[0,0] | L[0,0]/||l_0|| | Margin |
|--------|-----|---|---------|--------|----------------|--------|
| [14,44] | 1.5 | 5 | 1.426 | 1.426 | 100.0% | 1.71e-01 |
| [14,44] | 2.0 | 7 | 1.420 | 1.420 | 100.0% | 1.03e-01 |
| [14,44] | 2.5 | 9 | 1.417 | 1.417 | 100.0% | 4.6e-06 |
| [50,80] | 2.0 | 7 | 1.846 | 1.846 | 100.0% | 2.84e-01 |
| [50,80] | 3.0 | 11 | 1.548 | 1.548 | 100.0% | 1.07e-01 |
| [100,130] | 3.0 | 11 | 1.487 | 1.487 | 100.0% | 8.87e-02 |
| [100,130] | 4.0 | 16 | 1.166 | 1.166 | 100.0% | 7.75e-03 |
| [300,330] | 3.5 | 13 | 1.625 | 1.625 | 100.0% | 3.13e-01 |
| [300,330] | 4.5 | 17 | 1.460 | 1.460 | 100.0% | 1.06e-01 |

**Result**: In all 9 cases, L[0,1:] ≈ 0 (other Fejer coefficients are negligible).

---

## Structural Implication

### 1. Principal Vector = Scaled First Fejer Function
$$l_0(r) \approx \alpha_0 \cdot \varphi_0(r) = \alpha_0 \cdot \left(\frac{\sin(Wr)}{Wr}\right)^2$$

where $\alpha_0 = ||l_0|| \approx$ 1.1–1.8 (varies with window and 4W).

### 2. Geometric Interpretation
The **principal positivity** (the margin λ_min) is achieved by a **single band-limited function**, the zeroth Fejer kernel. This is remarkable because:

- The compressed operator Q is an N×N matrix (N = 5–17 depending on bandwidth)
- One would naively expect a dense Cholesky factor
- Instead, the **critical eigenvector (controlling positivity) is sparse**: L[0,:] has only one nonzero entry

### 3. Connection to Optimality
In the language of compressed operators:
- The **optimal test function** φ^* that minimizes Q(φ)/||φ||_G² (the Rayleigh quotient) is **essentially φ_0**
- The eigenvector v_0 = (1, 0, 0, ..., 0)^T (in the Fejer basis coordinates)
- This means the bottleneck is **a pure concentrated Fejer kernel**, not a complex combination

---

## Margin Scaling with Window Size

Across windows of different sizes, at comparable bandwidths:

| Window | L | 4W | ||l_0|| | margin |
|--------|-----|-----|---------|--------|
| [14,44] | 30 | 2.0 | 1.42 | 1.03e-01 |
| [50,80] | 30 | 2.0 | 1.85 | 2.84e-01 |
| [14,44] | 30 | 3.0 | — | — |
| [50,80] | 30 | 3.0 | 1.55 | 1.07e-01 |
| [100,130] | 30 | 3.0 | 1.49 | 8.87e-02 |

**Pattern**: As window location increases (height increases), ||l_0|| is **height-dependent** but shows a slow modulation. The margin is also height-dependent but more variable (due to local variations in Ω(t) and P(r)).

---

## Why Is This Happening?

### Hypothesis 1: Gram Matrix Structure
The Gram matrix G might have a **very dominant first eigenvalue** (or first few elements might be much larger than the rest). This would naturally push the Cholesky factor to concentrate in the first row.

### Hypothesis 2: Fejer Basis Optimality
The Fejer kernel φ_0 (centered at some c_0 ∈ [a,b], the first center point) might be a **near-optimal test function** for detecting violations of Q(φ) ≥ 0 in the band-limited sense. The first Fejer kernel captures the "flattest" mode, most prone to negativity.

### Hypothesis 3: Compression Artifact
When compressing Q from function space to the N-dimensional Fejer basis, the **first basis element** might encode the "dc component" or "mean behavior," which is the most sensitive to the pole-arch-prime competition.

---

## Theoretical Consequence

### Reduced Positivity Criterion
If l_0 ≈ α_0 φ_0 holds universally, we have:

$$Q(\phi) = \sum_{k=0}^{N-1} |\langle l_k, \phi \rangle|^2 \geq |\langle l_0, \phi \rangle|^2 = \alpha_0^2 |\langle \varphi_0, \phi \rangle|^2$$

This means: **Q(φ) ≥ α_0² · |⟨φ_0, φ⟩|² for all φ in the band-limited space.**

A **sufficient condition for Weil positivity** would then be to show:
$$\min_{\phi: |\phî(\xi)| = 0 \text{ for } |\xi| > 2W} \frac{\text{Re} Q(\phi)}{|\langle \varphi_0, \phi \rangle|^2} > 0$$

This is a **reduction from N variables to 1 effective variable** — a dramatic simplification.

---

## Next Steps

### 1. Verify Across Extended Range
Test 4W up to 6, 8, 10 to see if l_0 ≈ α_0 φ_0 persists into the RH-critical regime.

### 2. Bound α_0(4W)
Study the scaling α_0 = ||l_0|| as a function of 4W. Does it decay as 4W → ∞? How?

### 3. Cross-Window Comparison
For different window locations [a, b] at the same 4W, how much does α_0 vary? Is there a pattern (e.g., proportional to Ω(center) or density)?

### 4. Explicit Formula for α_0
Can we write α_0 in closed form, or at least as an integral over S(r) = (1/2π)Ω(r) - P(r)?

### 5. Proof Strategy Reduction
If α_0 is bounded below universally, the positivity proof simplifies to a **1-dimensional inner-product bound**.

---

## Tentative Theorem (Subject to Verification)

**Conjecture**: For all windows [a, b] and bandwidths 4W, the principal SOS vector satisfies:
$$l_0(r) = \alpha_0(a,b,4W) \cdot \varphi_0(r) + \text{(exponentially small higher terms)}$$

where $\alpha_0(a,b,4W)$ is a **smooth, explicitly computable function** that depends only on the window location (height), bandwidth, and the symbols Ω(t) and P(r).

If true, this **reduces the Weil positivity problem to a 1D function-space analysis**, a significant structural breakthrough.

---

## Proof Roadmap (If α_0 Dominance Holds)

1. **Compute α_0(t, 4W)** as a function of height and bandwidth (numerically and/or asymptotically)
2. **Prove α_0 ≥ δ > 0** universally in the RH-critical regime
3. **Show Q(φ) ≥ δ² · |⟨φ_0, φ⟩|² for all band-limited φ**
4. **Bound ⟨φ_0, φ⟩ from below** (optimal window/bandwidth choice)
5. **Conclude Weil positivity** = RH proved

This is speculative (step 2 is the hard part, equivalent to RH), but the **structural simplification** from N×N matrices to 1D function analysis is itself a new insight.

---

## Files for Reference

- **principal_sos_analysis.py**: Source of the discovery
- **CERTIFICATE_SOS_CONSTRUCTION.md**: Full construction details
- **sos_weil_construction.py**: Cholesky factor computation

---

## Conclusion

The principal SOS vector is **dominated by a single Fejer basis element**, a **structural simplification** unexpected from naive dimensional analysis. This opens the possibility of a **reduction to a 1-variable positivity problem**, and if α_0 has a closed-form expression, could lead to a **unified proof framework** for Weil positivity and ultimately RH.
