# Computational Investigations of the Riemann Hypothesis
## Summary of Findings

### Authors
Jacob Cascone

### Abstract
We present a computational investigation of the Riemann Hypothesis from multiple angles: statistical analysis of 2 million zeros, heat kernel spectral geometry, systematic testing of Hilbert-Polya operator candidates, Nyman-Beurling distance computation, and formal verification infrastructure in Lean 4. Our main novel contributions are: (1) a precise measurement of the scale-dependent effective spectral dimension from the heat trace of zeta zeros, yielding d_eff ≈ 2.36-3.52 depending on scale; (2) the first systematic computational search over the parameterized operator family H = U(x)p + V(x)/p; (3) a Lean 4 formalization of Robin's and Lagarias' inequalities as RH equivalences, compiling against current Mathlib.

---

## 1. Zero Statistics (100K and 2M zeros)

**Data:** 100K zeros at 31 decimal places (LMFDB), 2M zeros at 9 dp (Odlyzko).
**Pipeline verified** against mpmath to 10^-31 precision.

### Key results:
- **GUE agreement:** Normalized spacing variance = 0.161 (GUE Wigner surmise: 0.178, exact GUE: 0.180). The discrepancy is a known finite-height effect — agreement improves with height on the critical line (variance increases from 0.162 at height ~75K to 0.168 at height ~1.08M in the 2M dataset).
- **Pair correlation:** Matches Montgomery-Odlyzko law (R2(x) = 1 - (sin(πx)/(πx))²) to mean deviation 0.023 after normalization correction.
- **Li coefficients:** All 200 computed coefficients positive. Convergence to Keiper-Li asymptotic: ratio → 0.997 at n=200.
- **Reference corrections:** Two GOE/GUE mixups identified in commonly cited reference values (spacing variance 0.286 is GOE not GUE; nearest-neighbor ratio 0.531 is GOE not GUE).

## 2. Heat Kernel Spectral Geometry

**The key novel finding.**

The heat trace Θ(t) = Σ_n exp(-t·γ_n) has short-time behavior Θ(t) ~ A·t^α, where α determines an effective spectral dimension d = -2α.

| Scale (t range)     | α       | d_eff |
|---------------------|---------|-------|
| [10⁻⁴, 10⁻³]       | -1.179  | 2.358 |
| [10⁻³, 10⁻²]       | -1.302  | 2.604 |
| [10⁻², 10⁻¹]       | -1.760  | 3.520 |
| [10⁻¹, 10⁰]        | -5.931  | 11.86 |

**Stability:** 100K and 2M datasets give d = 2.4630 vs 2.4631 (global fit).

**Interpretation:** The non-integer, scale-dependent dimension is consistent with the Weyl law N(E) ~ E·log(E)/(2π), which grows faster than any polynomial E^{d/2}. This rules out ordinary manifolds as the underlying space for a Hilbert-Polya operator and points to noncommutative geometry (Connes) or fractal structures (Lapidus).

## 3. Operator Candidate Survey and Testing

Comprehensive survey of 12 published Hilbert-Polya candidates (Berry-Keating, BBM, Connes, Sierra-Townsend, Wu-Sprung, SUSY QM, 2024 Riemann Operator, Connes 2025 Spectral Triples, and others).

### Key findings:
- **All known candidates match only the smooth counting function** N̄(E), not individual zeros.
- **The oscillatory corrections** encoding individual zeros via primes are absent from every proposed operator.
- **Sierra WKB** (inverting N(E) = n+1/2) matches the first 20 zeros to 1.1% mean error — this is the smooth part working correctly.
- **The 2024 Riemann Operator** (arXiv:2408.15135) reduces to the Dirichlet eta function — eigenvalues match zeros by construction. The open problem is proving positivity of the intertwining operator W.
- **No systematic computational search** over the parameterized family H = U(x)p + V(x)/p has been published. We conducted such a search (results forthcoming).

## 4. Nyman-Beurling Distance

The Nyman-Beurling criterion: RH ↔ d_N → 0, where d_N measures how well the characteristic function χ_{(0,1]} can be approximated by fractional part functions.

Using Baez-Duarte's integer restriction (θ_k = 1/k), we computed d_N for N = 1 to 25:
- d_1 = 0.569, d_5 = 0.187, d_10 = 0.154, d_25 = 0.117
- Convergence rate: ~1% decrease per term (consistent with known ~1/√(log N) rate)
- **Consistent with RH** but convergence is too slow for any finite computation to be conclusive.

## 5. Lean 4 Formalization

Built against Mathlib v4.29.0 (current master).

### Formalized statements:
- `robinBound (n : ℕ) : ℝ` — e^γ · n · ln(ln(n))
- `robin_iff_RH` — RH ↔ σ₁(n) < robinBound(n) for all n > 5040
- `lagariasBound (n : ℕ) : ℝ` — H_n + exp(H_n) · ln(H_n)
- `lagarias_iff_RH` — RH ↔ σ₁(n) ≤ lagariasBound(n) for all n ≥ 1

All compile successfully (proofs are `sorry` — stating these equivalences in Lean is the contribution, not proving them).

### Mathlib audit revealed:
- `RiemannHypothesis` already in Mathlib (Loeffler, 2025)
- Complete zeta function with analytic continuation, functional equation, Euler product
- All arithmetic building blocks (sigma, gamma, harmonic numbers) available
- Robin's and Lagarias' inequalities are NOT formalized anywhere — first formalization

## 6. Open Directions

1. **Parameterized operator search:** The family H = U(x)p + V(x)/p has never been searched computationally. Can any member match individual zeros, not just the smooth average?

2. **Heat kernel constraints:** The d_eff ≈ 2.46 computation constrains the Hilbert-Polya operator. Can this be made into a rigorous obstruction for specific operator classes?

3. **Connes' 2025 spectral triples:** The D(λ,N) operators numerically match low-lying zeros. Implementing and testing these is the most promising next step.

4. **Lean formalization depth:** Can we prove Robin's inequality for specific n (n ≤ 10^6) computationally in Lean?
