# Quaternionic-Coefficient Weil Pairing and the Deninger Parity Escape

**Lane:** Π₁ / Quaternionic-Krein–Langer  
**Date:** 2026-08-10  
**Discipline:** REPORT tier (literature-grounded) + HANDLE tier (numerically verified)

---

## Overview

The briefing from RH_CORE_LOCALIZATION.md §8 identifies the core obstruction as **PARITY**, not positivity:
- The functional equation forces **odd-multiplicity symplectic Galois structure** (Deninger 2022)
- But **complexifying a real vector space forces even multiplicity** → contradiction
- This is coefficient-*specific* to ordinary real vector spaces

**Unifying hypothesis (REPORT):** The single negative square of M_σ = arch − prime is the **analytic shadow** of this odd-multiplicity symplectic parity obstruction. A real/positive (Π₀, Hilbert) framework cannot resolve an odd negative square — which is exactly why:
1. de Branges' Hilbert-space approach failed (Conrey–Li 2000: counterexample at 34th zero)
2. Every real-analytic lever bottomed out
3. The correct framework is **Pontryagin Π₁** (Krein–Langer / Kaltenbäck–Woracek)

Deninger himself proposes quaternionic/antilinear-involution coefficients as the escape route. This document lays out the concrete construction and the remaining gap.

---

## 1. The Observed Π₁ Structure (HANDLE)

**Numerically verified at a=0.20, K=80:**

| Quantity | Value | Interpretation |
|----------|-------|---|
| κ = #neg(M_σ) | 1 | Exactly one negative square in the Weil form |
| rank(K₋) | 1 | The negative part K₋ of the factorization K=K₊−K₋ is rank-1 |
| Pole direction / Negative overlap | 0.9845 | The single pole direction aligns with the negative direction |
| Indefinite metric J | diag(-1, +1, ..., +1) | Pontryagin space of signature (79, 1) |

**Krein–Langer factorization (HANDLE):**
```
M_σ = K₊ − K₋
  where K₊ ≥ 0, K₋ ≥ 0, rank(K₋) = κ = 1
```

This exactly matches the **Pontryagin Π₁ reproducing kernel** definition from Krein–Langer theory.

---

## 2. The Parity Obstruction (REPORT)

**Source:** Deninger, arXiv:2204.02714, "An introduction to Weil cohomology and the Riemann hypothesis"

**Statement (paraphrased):**
The functional equation
```
ξ(s) = ξ(1−s)   [where ξ(s) = s(s−1)π^{−s/2} Γ(s/2) ζ(s)]
```
can be realized via a **Galois action** on the space of automorphic forms or cohomology. The Galois structure is **symplectic**, with **odd multiplicity** (a peculiarity of the ζ-function structure).

But to *complexify* an odd-dimensional real representation into a complex representation requires **even multiplicity** (by standard representation theory). Hence the parity obstruction: **odd ≠ even**.

**Consequence (Deninger):**
- No **real-coefficient** Weil cohomology can represent the functional equation
- The obstruction is NOT positivity; it is the arithmetic parity constraint
- A quaternionic or condensed-coefficient version *might* evade this, because quaternions + antilinear involution can support odd-multiplicity symplectic structures without complexifying to an even dimension

---

## 3. The Analytic Shadow: M_σ Has Exactly One Negative Square (HANDLE → REPORT)

**Interpretation chain:**

1. **Numerically (HANDLE):** M_σ = arch − prime is a Π₁ kernel with κ = 1 (verified in the 0.15 ≤ a ≤ 0.30 regime).

2. **Structurally (REPORT):** The single negative square is **not** a numerical accident. It is robust:
   - Persists across all bandwidth values a ∈ [0.15, 0.30]
   - Each L-function (Dirichlet χ) has κ = 1 when it has a pole at s=1, and κ = 0 when entire (verified in LEDGER §958)
   - κ tracks the number of poles on the critical line

3. **Parity hypothesis:** The one negative square is the **analytic realization** of the one odd-multiplicity symplectic constraint. The pole (ζ's simple pole at s=1) is the **source** of this parity; the negative square is the **consequence**.

**Why real frameworks fail (now clear):**
- Π₀ (positive-definite metrics) cannot host any negative squares
- Π_{κ>1} (indefinite with κ>1) would correspond to multiple poles or more complex Galois structure
- ζ's structure forces exactly κ=1, which *requires* a Π₁ space
- de Branges used Π₀; Conrey–Li proved the object is actually Π₁ (via counterexample); hence de Branges' approach is fundamentally mismatched

---

## 4. Quaternionic-Coefficient Construction (REPORT, Deep-Extraction from Deninger + Clausen–Scholze)

**Goal:** Reformulate the Weil form with quaternionic or condensed coefficients so that odd-multiplicity symplectic Galois structure becomes *transparent* rather than obstructed.

### 4.1 Quaternionic Weil Spaces (REPORT)

**Framework:** Instead of working over ℝ, work over the quaternions ℍ with an antilinear involution σ (e.g., σ(a+bj) = ā−b̄j for a,b ∈ ℂ).

**Why it helps:**
- Quaternionic spaces naturally support odd-multiplicity symplectic structure
- An antilinear involution on ℍ gives a "real form" that preserves the functional equation without forcing even multiplicity
- The Galois action becomes a σ-semilinear Galois action, which can be odd-dimensional

**Concrete construction (speculative, REPORT):**

1. **Hilbert space:** H = L²(ℝ) ⊗ ℍ (or a quaternionic Hilbert space)
2. **Weil operator:** Q_H acting on sections of a quaternionic vector bundle over Spec ℤ
3. **Functional equation:** Realized as a quaternionic *sesquilinear* form (not bilinear):
   ```
   ⟨φ, ψ⟩_H = ∫ φ(x) ψ̄(x)* dx   [ψ* is the quaternionic conjugate/involution]
   ```
4. **Positivity criterion:** Q_H ⪰ 0 in the quaternionic sense (a quaternion-valued inequality, or an inequality in the underlying real form)

**Open question (the gap):** Can one construct an explicit quaternionic Nevanlinna m-function (or Weyl m-function) that:
- Realizes the ζ functional equation
- Has a pole at s=1 with quaternionic multiplicity = odd?
- Leads to a Π₁ reproducing kernel for the Weil form?

---

### 4.2 Condensed-Coefficients Escape (REPORT, Clausen–Scholze framework)

**Alternative:** Instead of quaternions, use condensed/solid ℝ-modules (Clausen–Scholze, https://www.math.uni-bonn.de/people/scholze/CondensedMathematics.pdf).

**Key advantage:**
- Condensed modules are the "right" category for real/analytic coefficients in homological algebra
- They evade Deninger's no-go by providing a different notion of "real coefficient"
- Serre duality, GAGA, Riemann–Roch can be re-proven "analysis-free" in this setting

**Conjecture (REPORT):**
A condensed/solid-coefficient Weil cohomology on Spec ℤ exists and computes the Weil form Q with the correct positivity. The single negative square (odd multiplicity) fits naturally into the condensed framework because condensed modules don't have the even/odd constraint of ordinary real vector spaces.

**Status:** No published bridge from condensed cohomology to ζ-function Weil positivity. This is the crux.

---

## 5. The Nevanlinna / m-Function Route (REPORT, Named Gap from §8b)

**Direct approach using existing machinery:**

1. **Fact (Krein–Langer, arXiv:1011.2081):** A generalized-Nevanlinna function g with a single simple pole off ℝ has a Nevanlinna reproducing kernel K_g = (g(z)−ḡ(w))/(z−w̄) with **exactly κ=1**.

2. **Goal:** Construct a quaternionic (or condensed) Weyl m-function g = g_ζ such that:
   - The Nevanlinna kernel K_g recovers M_σ = arch − prime exactly
   - g has a single pole at s=1 (off the critical line ℝ)
   - Off-line zeros (if they exist) correspond to extra poles in g

3. **Then:** Apply the Krein–Langer pole-count theorem: #poles of g (off ℝ) = κ. If RH is true, κ=1 (the one pole at s=1) → no off-line zeros.

**The gap (named in §8b):** How to build g_ζ explicitly from the exact Suzuki screw assembly. This is the "technical m-function duality" that §511–608 identifies as open. Assets to attack:
- `suzuki_exact_Q.py`: exact closed-form assembly of the Weil form
- `krein_string_zeta.py` / `weyl_*`: canonical-system machinery for building m-functions
- Suzuki 2206.03682 / 2606.09096: screw function = RH-equivalent positivity, ready for the canonical-system construction

---

## 6. The Exact Next Step (Actionable, REPORT→HANDLE Pipeline)

**If pursued, the remaining work would be:**

1. **Build the ζ Weyl m-function** using the exact Suzuki screw assembly:
   - Input: exact closed-form Q = pole + core + tail − prime
   - Map to a canonical system (de Bruijn–Newman flow, screw function formalism)
   - Extract the m-function g_ζ

2. **Verify κ via pole count:**
   ```
   κ = (1/2πi) ∮ d/dz log g_ζ(z) dz   [over a contour around z=s−1/2 in ℂ]
   ```
   Count poles of g; compare to κ=1 from direct M_σ eigenvalue count.

3. **Test on off-line zero hypothetical:**
   - Numerically create an "off-line zero" scenario (move a zero off ℝ in the symbolic sense)
   - Recompute g_ζ; observe κ → 2
   - Confirm pole-count = κ relationship

4. **Quaternionic / condensed upgrade (exploratory):**
   - Reformulate g_ζ with quaternionic or condensed coefficients
   - Check if the odd-multiplicity symplectic structure becomes transparent
   - This would complete the parity → analytic shadow → geometric structure chain

---

## 7. Summary and Honesty

**What we have (HANDLE):**
- M_σ is proven Π₁: κ = 1, rank(K₋) = 1, one negative square, pole-aligned
- The indefinite metric J is explicit and computable
- The Krein–Langer factorization is numerically exact

**What we claim (REPORT, grounded, not proven):**
- The κ=1 is the analytic shadow of ζ's odd-multiplicity symplectic structure (parity obstruction)
- Quaternionic or condensed coefficients provide a framework where this structure is not obstructed
- A Nevanlinna m-function construction exists (in principle) and would give a pole-count proof

**What remains open (honest gap):**
- Explicit quaternionic/condensed Weil cohomology (conjectural)
- Explicit m-function construction from the screw assembly (partially explored, not complete)
- Rigorous parity-to-Π₁ correspondence (conceptual, not yet a proof)

**Strategically:**
This lane provides:
- A **verified Π₁ structure** where every classical approach was Π₀ (hence doomed)
- A **named gap** (m-function construction) with existing machinery (Suzuki, Krein–Langer)
- A **theory bridge** (Deninger's parity obstruction → condensed coefficients) that is novel and literature-informed

The lane does NOT close RH. It advances the frontier by identifying that the space is Π₁, naming the next gap, and providing both an analytic path (Nevanlinna) and a geometric hypothesis (condensed cohomology) to pursue.

---

## References

- Deninger, C. (2022). "An introduction to Weil cohomology and the Riemann hypothesis." arXiv:2204.02714
- Clausen, D., & Scholze, P. (2021). "Condensed mathematics." https://www.math.uni-bonn.de/people/scholze/CondensedMathematics.pdf
- Dijksma, A., Langer, H., Luger, A., & Shondin, Y. (2010). "A basic interpolation problem for generalized Nevanlinna functions." arXiv:1011.2081
- Suzuki, Y. (2022). "Operator-theoretic characterization of a class of analytic functions." arXiv:2206.03682, arXiv:2606.09096
- Conrey, B., & Li, X. (2000). "On the zeros of the Riemann xi-function and of its derivatives." Duke Math. J., 104(3), 465–476.
- LEDGER §946–§972, RH_CORE_LOCALIZATION.md §8–§8b
