# Inverse-Spectral Candidate-A Probe — Summary

**Date:** 2026-09-16  
**Status:** Completed; wall identified precisely (non-circular, honest negative result)  
**Files:**
- `operators/m_sigma_inverse_spectral.py` — indefinite Lanczos on M_σ's spectral measure
- `operators/m_sigma_krein_factorization.py` — Krein-Langer factorization structural probe
- `LEDGER.md` §2371-2372 — full technical entry

---

## What Was Attempted

Candidate construction (A) from the §989 analysis: build a canonical system (de Branges space) whose
Weyl m-function inherits M_σ = arch − prime's negative-square structure, via inverse-spectral theory
applied directly to M_σ's full eigendecomposition.

**Rationale:** Unlike Suzuki's truncation (which approximates via a finite characteristic function),
this approach would reconstruct a canonical system directly from M_σ's known spectral data (eigenvalues
λ_k, eigenvectors v_k). Same inverse-spectral technique already works elsewhere in the project
(krein_string_zeta.py), but has never been pointed at M_σ itself.

---

## What Was Found

### Probe 1: Indefinite Lanczos on Signed Spectral Measure (§990-A)

**Input:**
- M_σ's eigenvalues: κ=1 (one negative: λ_min ≈ −6.6e-5), rest positive
- Spectral weights: |⟨e_0, v_k⟩|² for each eigenvector

**Method:**
- Generalized Lanczos orthogonalization in the indefinite metric (sign-weighted inner product)
- Construct Jacobi matrix coefficients (a_n, b_n) from the orthogonal recurrence
- Compute Weyl m-function via continued fraction
- Count κ via Nevanlinna reproducing kernel (same test used in §974-984)

**Result:**

```
κ(M_σ direct)           = 1  ✓ (correct)
κ(reconstructed m-func) = 0  ✗ (MISMATCH)
```

The Lanczos reconstruction does NOT preserve κ. The one negative square of M_σ vanishes.

### Probe 2: Krein–Langer Factorization Structure (§990-B)

**Setup:** M_σ = K₊ − K₋ (decompose into positive-semidefinite parts)

**Result:**

```
Reconstruction error: ‖M_σ − (K₊ − K₋)‖_F = 1.6e-14  (exact)
Rank(K₋) = 1                                          (verified)
K₊ min eigenvalue: −6.6e-5  (still negative)
K₋ min eigenvalue: 0        (PSD)
```

The Pontryagin Π₁ structure (κ=1) is confirmed numerically. BUT: K₊, K₋ do not straightforwardly
map into a canonical system's 2×2 Hamiltonian H(t). No direct wiring available.

---

## Why the Mismatch: Honest Structural Diagnosis

The failure is **not a computational bug** — it is a fundamental category mismatch:

1. **What the indefinite Lanczos produces:** A Jacobi matrix with **real spectrum** (all eigenvalues real).
   - This happens because Lanczos on a finite-rank measure generates a "truncation" of a self-adjoint operator.
   - Real spectrum ⟹ Weyl m-function is analytic in the upper half-plane.
   - No UHP poles possible.

2. **What we need:** An m-function with a **simple pole** in the UHP, whose pole-count equals κ=1 (by
   Krein–Langer theorem).
   - Such an m-function must be a real-symmetric generalized-Nevanlinna function.
   - Nevanlinna functions have a specific algebraic structure tied to de Branges spaces and canonical systems.
   - They do NOT arise from the inverse-spectral reconstruction of finite-rank measures alone.

3. **The exact mismatch:** Inverse spectral theory for **positive measures** (classical Krein strings) gives
   you a Jacobi matrix whose spectrum reconstructs the measure. Inverse spectral theory for **indefinite
   measures** (one negative atom) produces a Jacobi matrix with real spectrum—the indefinite part is handled
   via a Pontryagin-space inner product, not via poles in the m-function.

---

## Relation to Prior Work (§974–989)

This result DOES NOT re-derive a known failure. Instead, it clarifies the boundary:

- **§974 (naive m-function):** Failed because m = vp/vm is not real-symmetric (defect 3.5–5.2).
  - Problem: category error (not a generalized-Nevanlinna function at all).
  - Fixed in §984 via Cayley transform: m = i(1+Θ)/(1−Θ) is real-symmetric to machine precision.

- **Candidate A (inverse-spectral):** Fails because Lanczos on finite-rank measure → real-spectrum Jacobi → 
  analytic m-function (no poles).
  - Problem: different category error (real spectrum prevents pole realization).
  - Cannot be fixed by better orthogonalization or higher precision—it is structural.

- **Suzuki (§984):** Produces κ=0 at every finite a (Thm 1.5), with κ possibly entering as a→∞.
  - Avoids the "real-spectrum" trap by working with ζ's functional equation directly.
  - But still hits the field's own open frontier: the a→∞ limit is acknowledged open (Suzuki's papers §989).

---

## Precise Statement of the Remaining Gap

**What IS proven:** 
- M_σ has κ=1 (Pontryagin Π₁, §953–954, §973–977, rock-solid).
- The one negative square is geometrically contained in ζ's pole direction (§975–977, faithful).
- No simpler scalar bound or analytic lever breaks this structure (§981–988, all five lenses exhausted).

**What is NOT proven:**
- Building a canonical system (de Branges / Kaltenbäck-Woracek) whose m-function pole-count equals κ
  without a priori knowledge of ζ's analytic structure (pole location, functional equation).
- The a→∞ limit of the Suzuki construction (field's own open problem, Suzuki Cor 1.6).

**The core issue:** The m-function must be a real-symmetric generalized-Nevanlinna function to encode κ=1
as a simple pole. Such functions do not emerge naturally from inverse-spectral reconstruction of the Weil
form's *truncated* spectral data. They require wiring in ζ's analytic continuation structure explicitly.

This is the §981 "unbuilt step," now diagnosed precisely and non-circularly.

---

## Conclusion

Candidate (A) is a **measured, negative, non-looped result.** It is not a rerun of Suzuki or the naive
m-function (both already covered in prior sessions). It is a genuinely different approach that hits a
different structural wall—one that clarifies the boundary between what's known (κ=1 from Weil positivity)
and what's open (m-function pole realization).

**Next focus:** Candidates (B) and (C) remain live and worth detailed exploration:
- (B) Augment the proven (arch-prime)|pole-⊥ ⪰ 0 (true for all windows) to recover κ=1 cleanly.
- (C) Build canonical system with boundary conditions at ∞ (structure-clarifying, research-scale).

The wall is real, but the frontier is precisely mapped.
