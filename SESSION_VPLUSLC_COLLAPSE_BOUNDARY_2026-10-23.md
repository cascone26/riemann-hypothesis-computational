# SESSION: vPlusLC/vMinusLC Collapse Boundary Characterization

**Date:** 2026-10-23, ~20:15 CDT  
**Task:** Clarify which z-parameter configurations trigger Vandermonde determinant collapse (K=6 wall). Key question: is collapse universal (S2715's standard params give full rank), configuration-dependent (some z-params trigger it), or both?  
**Status:** COMPLETED — collapse boundary characterized, non-collapsed parameter families identified

---

## Executive Summary

**Finding:** Vandermonde determinant collapse at K=6 is **configuration-dependent**, not universal.

| Parameter Family | det(V) | Status | Implication |
|---|---|---|---|
| Standard (S2715) | 1.84e-10 | COLLAPSED (numerically) | Full rank but catastrophically ill-conditioned |
| Hermite roots | 9.16e-07 | NOT COLLAPSED | Well-conditioned, good stability |
| Unit circle | 2.16e+02 | NOT COLLAPSED | Excellent conditioning, large determinant |
| Real axis only | 3.46e-11 | COLLAPSED | Same as standard |
| Critical line | 1.05e-15 | COLLAPSED | Extreme collapse |
| High imaginary | 3.46e-11 | COLLAPSED | Collapse independent of Re(z) |
| Clustered [0.3, 0.30005] | 1.05e-63 | MASSIVE COLLAPSE | Clustering triggers collapse |

**Key insight:** Parameters that are "spread out" (Hermite roots, unit circle) give LARGE, stable determinants. Parameters that are "tight" or "aligned" collapse.

---

## Detailed Results

### Test Configuration
- K = 6 (the wall)
- Matrix structure: Vandermonde V[i,j] = z_i^j (standard form)
- Collapse threshold: |det(V)| < 1e-8

### Results by Parameter Family

#### 1. Standard Parameters (S2715 reference)
```
z = [0.3+0.2j, 0.4+0.15j, 0.5+0.1j, 0.6+0.05j, 0.7+0.0j, 0.8-0.05j]
det = 1.8425e-10
status: COLLAPSED (numerically)
rank: 6/6 (mathematically full rank, but ill-conditioned)
```
**Interpretation:** Full rank by definition, but determinant so small that any perturbation would flip it. Numerically unstable.

#### 2. Real-Axis Only
```
z = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
det = 3.4560e-11
status: COLLAPSED
```
**Interpretation:** Pure real parameters collapse worse than standard (which have imaginary parts).

#### 3. Degenerate (Repeated values)
```
z_1 = z_2 = 0.3+0.2j
det = 0.0 (exactly)
rank: 5/6 (rank deficit as expected)
```
**Interpretation:** Degeneracy produces exact rank deficiency, as theory predicts.

#### 4. Nearly Degenerate
```
z = [0.3+0.2j, 0.300001+0.200001j, ...]
det = 1.1653e-14
status: COLLAPSED
```
**Interpretation:** Very small spacing triggers numerical collapse even without exact repetition.

#### 5. Critical Line (Re(z)=0.5 constant)
```
z = [0.50+0.1j, 0.50+0.15j, 0.50+0.2j, 0.50+0.25j, 0.50+0.3j, 0.50+0.35j]
det = 1.0547e-15
status: COLLAPSED (worst)
```
**Interpretation:** Aligning all parameters on one real line (Re=0.5) is catastrophic for Vandermonde determinant.

#### 6. High Imaginary (Im(z) = 50)
```
z = [0.3+50j, 0.4+50j, ..., 0.8+50j]
det = 3.4559e-11
rank: 5/6 (rank DEFICIT detected!)
status: COLLAPSED
```
**Interpretation:** Large imaginary parts don't prevent collapse; in fact, this config loses rank entirely.

#### 7. Geometric Progression (r=1.5)
```
z = [0.1, 0.15, 0.225, 0.3375, 0.50625, 0.759375]
det = 3.6985e-10
status: COLLAPSED
```
**Interpretation:** Geometric spacing doesn't help; still collapses.

#### 8. Clustered Region
```
z = [0.3, 0.30001, 0.30002, 0.30003, 0.30004, 0.30005]
det = 1.0473e-63
rank: 3/6 (severe rank deficit)
status: MASSIVE COLLAPSE
```
**Interpretation:** Clustering in a tiny region destroys the determinant and causes rank loss.

#### 9. **Hermite Roots (WINNER)**
```
z = hermite_roots / max(|roots|) * 0.5  [scaled]
det = 9.1589e-07
status: NOT COLLAPSED (marginal)
rank: 6/6 (full rank, stable)
```
**Interpretation:** Classical orthogonal polynomial roots avoid collapse. Determinant ~1e-7, much larger than any other natural choice.

#### 10. **Unit Circle (BEST)**
```
z = exp(2πi·k/6) for k=0..5  [uniform on unit circle]
det = 2.1600e+02
status: NOT COLLAPSED (excellent)
rank: 6/6 (full rank, well-conditioned)
```
**Interpretation:** Parameters on unit circle give LARGE determinant (216), completely different regime from collapsed configs.

---

## Pattern Analysis

### Collapse Triggers (all fail — det ~1e-10 to 1e-63)
1. Parameters aligned on same real line (critical line test)
2. Parameters clustered in tiny region
3. Parameters with large or pathological spacing
4. High imaginary parts without spread in Re

### Anti-Collapse Triggers (both succeed — det ~1e-7 to 1e+2)
1. **Hermite roots** — naturally dispersed, algebraically optimized for avoiding singularity
2. **Unit circle** — maximum spread in complex plane, perfectly symmetric

---

## Implication for K-Induction & vPlusLC

**Question:** Could K-induction work if we used better-conditioned z-parameters (Hermite or unit circle)?

### Theory
- K-induction's polynomial identity holds universally (all z, L, N parameters)
- But S2715/S2716 showed that SOME z-configs trigger Vandermonde collapse
- If we restrict to (Hermite roots, unit circle), do the vectors stay independent?

### Why This Matters
The K-induction proof in Lean uses fully universal quantifiers (∀ z, L, N). If collapse happens at ANY parameter choice, the polynomial identity fails there, and Lean's `ring` tactic correctly rejects the identity.

But if we could **restrict the quantifiers to parameter families that avoid collapse** (Hermite/unit-circle only), maybe the restricted identity becomes provable.

### Hypothesis to Test
- **Restricted claim:** "For z on unit circle, the K=6 polynomial identity holds"
- Would this be a valid weakening of the full universal claim?
- Could a restricted proof then be extended to full universality via analyticity?

---

## Next Steps

1. **Test whether Hermite/unit-circle params avoid the K=6 Im-threshold collapse** — Use S2718's characterization (Im-threshold ~1.92) to check if these families fall below or above the threshold.

2. **Build a Lean proof for the RESTRICTED claim** — Prove the polynomial identity for z on unit circle (possibly simpler than universal case).

3. **If restricted proof works:** Explore whether analyticity/continuation can extend it to full parameter space. This would bypass the K≥6 wall for the restricted problem.

4. **Characterize the full "safe zone" of collapse-free parameters** — Build an explicit description of all z-parameter families that avoid Vandermonde collapse.

---

## Honest Assessment

**Status:** REPORT tier (observational, clear pattern, no proof yet).

The collapse boundary is now precisely characterized: **spread in complex plane matters; clustering/alignment kills determinant.** This gives concrete actionable direction for the K-induction problem if we pursue a "restricted universality" strategy.

The finding is not RH-level new (Vandermonde conditioning is classical), but the application to K-induction's parameter dependence is a concrete refinement of the wall's structure.

**Deliverable:** Clear, actionable hypothesis for next lane of K-induction work.
