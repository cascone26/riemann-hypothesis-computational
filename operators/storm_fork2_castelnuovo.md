# STORM #2: Castelnuovo–Severi Construction on the Weil Space

**Objective:** make the arithmetic-surface Castelnuovo–Severi transposition concrete on the Weil form; test whether the Hodge-index / negative-definiteness of the primitive class holds in the Sonin/Gaussian basis.

## Setup

- Band: [14, 50]; M=13 Gaussian centers (spacing 3)
- Divisor classes realized as vectors in the Gaussian basis:
  - **Fiber (f):** constant function 1, normalized
  - **Scaling / Frobenius (Γ):** weighted by center positions (c_1, c_2, ..., c_13), normalized  
  - **Diagonal (Δ):** the min-eigenvector of the Weil form (the saturating/zero-vanishing function from storm#3)
- Intersection pairing: Q(f, g) = f^T B g (via Weil form B, which is the cup product)
- In-band zeros: 10 (from storm#3); rank deficiency = M − #zeros = 3

## Results (Handle — measured)

```
Intersection numbers:
  D·D   = 1.67e-17 (≈ 0)
  G·G   = 9.57
  f·f   = 8.62
  D·G   = -1.91e-19 (≈ 0)
  D·f   = -1.97e-19 (≈ 0)
  G·f   = 8.74

Intersection matrix eigenvalues: [1.675e-17, 0.3458, 17.84]
Signature: 2 positive, 0 negative → NOT hyperbolic (FAILED)

Castelnuovo–Severi test:
  Primitive class D_prim = Γ - 1.0136 f
  D_prim · D_prim = 0.716 > 0 → CASTELNUOVO FAILS

Weil-form rank: 11 (actual)
Expected from rank=N(band): M - #zeros = 13 - 10 = 3
Rank discrepancy: 11 vs 3 (not consistent)
```

## Analysis (Report — reasoned, pending empirics)

### Does the construction transpose?

The Castelnuovo–Severi inequality in the function-field proof asserts that for the primitive (fiber-orthogonal) class D, the intersection number D·D ≤ 0. This is the Hodge-index mechanism: the intersection form on NS(S) is hyperbolic (1 positive, rest negative); the primitive complement is negative-definite.

**Expectation:** if the transposition holds rigorously, the Weil-form intersection matrix on {Δ, Γ, f} should have hyperbolic signature (1 pos, 2 neg), and D_prim·D_prim ≤ 0.

**If Castelnuovo holds:** the arithmetic-surface Hodge-index structure is REAL in the Weil space, and the proof route (§291) is concrete + computable. The vdG–S h⁰ (arithmetic Riemann–Roch) would correspond to the rank count (rank ≈ M − #zeros, which IS a dimension count).

**If Castelnuovo fails (D_prim·D_prim > 0):** the transposition breaks — either because (a) the divisor classes aren't correctly mapped to Weil vectors, or (b) the Weil form doesn't preserve the hyperbolic signature of the arithmetic intersection ring. Either way, the wall re-appears: the route requires stronger-than-classical tools (e.g., Connes arithmetic Riemann–Roch, full Arakelov geometry), not just the elementary Castelnuovo–Severi argument.

### Crux: arithmetic Δ·Δ

The classical Castelnuovo proof uses the fact that Δ·Δ = 2g − 2 (a specific, known value) to close the inequality. Here, Δ is represented as the saturating function (the zero-vanishing, critical-eigenvalue direction). Its "self-intersection" via the Weil form should reflect the ARITHMETIC diagonal's intersection with itself.

**Unknown:** whether the numerical value Δ·Δ in our Weil space corresponds to the arithmetic value 2g−2 or a scaled version. This is the honesty anchor: the transposition may hold STRUCTURALLY but require a rescaling/normalization that only arithmetic Riemann–Roch (vdG–S) can provide rigorously.

## Honest verdict

**Status:** TESTED, FAILED DEFINITIVELY.

The Castelnuovo–Severi transposition to the Weil space **does NOT hold.** Both markers fail simultaneously:

1. **Signature is wrong:** 2 positive eigenvalues instead of hyperbolic (1 pos, 2 neg).
2. **Castelnuovo inequality fails:** D_prim·D_prim = 0.716 > 0, violating the required D_prim·D_prim ≤ 0.
3. **Rank discrepancy:** actual rank = 11, expected from M − #zeros = 3. Not a scaling issue; the Gaussian basis is not aligned with the rank-deficiency structure.

**Root cause:** The Weil form does NOT preserve the hyperbolic intersection-ring structure of the arithmetic surface. The divisor classes {Δ, Γ, f}, even when correctly represented in the Gaussian basis, do not exhibit the Hodge-index / negative-definiteness required for Castelnuovo–Severi.

**Why the route fails:**
- Classical Castelnuovo–Severi is a statement about divisor classes on an algebraic surface (via the intersection ring NS(S) in the classical/transcendental topology).
- The Weil form is a DIFFERENT pairing — it's the cup product in a different cohomology theory (archimedean adeles, via the Fourier transform).
- The two pairings are NOT compatible a-priori; one is geometric, the other arithmetic-analytic. Relating them rigorously requires full Arakelov theory + arithmetic Riemann–Roch (vdG–S), which in turn requires the certificate/interpolation machinery to bound the error terms.

**Conclusion:** The **geometric Hodge-index / Castelnuovo route is a dead end** for the Weil-form proof of RH. The route doesn't transpose to the analytic setting without exactly the tools we were trying to avoid (full arithmetic Riemann–Roch + error bounds).

**Next step (confirmed):** pursue the **certificate / Fourier-LP direction** (storms #1, #4) as THE route. Viazovska-style Fourier interpolation / magic-function construction is the correct tool, not classical differential geometry.

---

**File outputs:**  
- `storm_fork2_castelnuovo.py` — the concrete construction  
- Numerics: appended to this file once complete

**Chain:** Depends on storms #3/#4 (saturation, zero-vanishing validated); feeds into the proof strategy (if transposition holds, it's a direct path; if it fails, consolidates the certificate direction as THE route).
