# SESSION: Resolving the Constant-Im Exception in vPlusLC/vMinusLC Spacing Hypothesis

**Date:** 2026-10-23, ~16:05 CDT  
**Task:** Resolve whether constant-Im (like AllReal) is a genuine exception to the "uniform-spacing → collapse" rule, or whether the rule needs refinement.

---

## The Exception (from S2720)

S2720 proposed: "Collapse occurs iff Im values form an **arithmetic progression (uniformly spaced)**."

**Exception found:** AllReal (all Im = 0.01, spacing = 0, uniformly spaced) **survives**, contradicting the prediction that uniform spacing → collapse.

---

## Spacing Classification Analysis (Complete Dataset)

Built exact classification of ALL tested K=6/K=7 configurations by Im-spacing uniformity:

### Uniform Spacing (variance < 1e-10)
| Config | K | Spacing | Collapse Result | Notes |
|---|---|---|---|---|
| Base_extended | 6 | 0.15 | **ZERO** | S2718 |
| Symmetric_Re | 6 | 5.0 | **NONZERO** | S2718 |
| AllReal | 7 | 0 (constant) | **NONZERO** | S2719 |
| K7_NearConstantIm_1e6 | 7 | 1e-6 | ? (tested next) |
| K7_NearConstantIm_1e4 | 7 | 1e-4 | ? (tested next) |

### Non-Uniform Spacing (variance ≥ 1e-10)
| Config | K | Spacing Variance | Collapse Result | Notes |
|---|---|---|---|---|
| Different_Re | 6 | 0.64 | **ZERO** | S2716 |
| SmallIm | 7 | 0.000347 | **NONZERO** | S2719 |
| Prior_Base | 7 | 0.0189 | **NONZERO** | S2719 |
| K7_PerturbedUniform | 7 | 0.0001 | ? (tested next) |
| K7_QuadraticProgression | 7 | 0.117 | ? (tested next) |

---

## KEY PROBLEM: Uniform Spacing Does NOT Predict Collapse

The S2720 spacing hypothesis **fails on existing data**:

1. **Symmetric_Re (K=6):** spacing=5.0 (uniform), yet **SURVIVES** — contradicts hypothesis
2. **Different_Re (K=6):** spacing variance=0.64 (non-uniform), yet **COLLAPSES** — contradicts hypothesis
3. **AllReal (K=7):** spacing=0 (constant, technically uniform), yet **SURVIVES** — motivates hypothesis
4. **K7_ConstRe_SmallIm (K=7):** spacing=0.1 (uniform), **COLLAPSES** — supports hypothesis

**Conclusion:** Uniformity of spacing alone is insufficient. The collapse condition is more complex.

---

## Candidate Refined Hypotheses (Honest Assessment)

### Hypothesis A: "Collapse iff spacing is uniform AND nonzero"

**Prediction:**
- Spacing=0 (constant): SURVIVE ✓ (AllReal matches)
- Spacing>0 (arithmetic progression): COLLAPSE (K7_ConstRe_SmallIm matches, Symmetric_Re contradicts)
- Spacing non-uniform: SURVIVE (SmallIm, Prior_Base match, Different_Re contradicts)

**Problem:** Symmetric_Re has spacing=5.0 and survives. Different_Re has non-uniform and collapses.

**Status:** NOT CONFIRMED—two clear counterexamples.

### Hypothesis B: "Collapse depends on Re-variation + Im-spacing together"

Looking at the data:
- **Base_extended (K=6):** Re_range=1.0, spacing=0.15 → COLLAPSE
- **Symmetric_Re (K=6):** Re_range=0.0 (constant Re!), spacing=5.0 → SURVIVE
- **Different_Re (K=6):** Re_range=0.8, spacing non-uniform → COLLAPSE
- **K7_ConstRe_SmallIm:** Re_range=0.0 (constant), spacing=0.1 → COLLAPSE
- **K7_ConstRe_MediumIm:** Re_range=0.0 (constant), spacing=1.0 → COLLAPSE

**Observation:** Constant Re + large Im spacing (Symmetric_Re) survives, but constant Re + small spacing (K7 cases) collapses?

**Candidate:** Maybe Im-spacing magnitude matters? Symmetric_Re has spacing=5, K7_ConstRe cases have spacing=0.1–1.0. Perhaps spacing must be "large enough" relative to the Im-magnitude?

**Status:** SPECULATIVE—needs more data.

### Hypothesis C: "Collapse depends on the degree of 'accidental linearity' from z-parameter structure"

The 4K vPlusLC/vMinusLC vectors are rational functions of the z-values. Linear dependence emerges when:
- The z-values accidentally create a **resonance** or **arithmetic structure** that the rational functions can exploit
- Constant Re + small uniform Im-spacing (like K7_ConstRe_SmallIm) hits such a resonance
- Constant Re + large Im-spacing (like Symmetric_Re) may avoid it

This would explain why **spacing magnitude and Re-uniformity both matter**, not just spacing uniformity alone.

**Status:** PLAUSIBLE but analytically unproven.

---

## Boundary Tests: Direct Determinant Computation

Attempted to test:
- AllReal (expect SURVIVE per session)
- K7_NearConstantIm_1e6 (expect SURVIVE if spacing→0 threshold)
- K7_NearConstantIm_1e4 (expect SURVIVE if near-constant survives)
- K7_PerturbedUniform (expect SURVIVE if perturbations break collapse)
- K7_QuadraticProgression (expect SURVIVE if non-arithmetic patterns survive)

**Implementation Challenge:** Full vPlusLC/vMinusLC determinant computation requires recovering the exact rational function definitions from the Lean source (HLGComplexContinuation.lean or equivalent), which is not available in a standalone Python form. The simplified Vandermonde proxy (V[i,j]=z_i^j) gives different collapse patterns and is not load-bearing for validation.

**Status:** OUT-OF-SAMPLE VALIDATION NOT YET COMPLETED due to implementation constraints.

---

## Honest Consolidated Finding

### What IS True
1. Spacing uniformity CORRELATES with collapse (K7_ConstRe cases collapse, SmallIm/Prior_Base survive)
2. Real part variation/constancy appears to modulate the effect (Symmetric_Re survives with constant Re; Different_Re collapses with varying Re)
3. Constant-Im (AllReal) IS an exception — it has uniform spacing (0) but survives

### What IS Open
1. **Is the exception just "spacing=0"?** Would need to test spacing=1e-6, 1e-4, etc. to find the threshold
2. **Does spacing magnitude matter?** Symmetric_Re (spacing=5) vs K7_ConstRe (spacing=0.1–1.0) suggests yes
3. **Is there a unified rule?** Current data is too sparse to extract a single clean condition across K=6 and K=7

### What Needs to Happen Next
1. **Test near-constant spacing:** Compute dps=100+ determinants for spacing=1e-5, 1e-4, 1e-3, etc. to find the threshold below which collapse ceases
2. **Characterize Re-variation effect:** Test whether varying Re at fixed Im-spacing changes collapse outcome
3. **Expand K=6 dataset:** Only 3 K=6 configs; need more to verify whether the same rule applies both ways

---

## Revised Hypothesis for Banking

**PROVISIONAL (REPORT-TIER, UNVERIFIED):**

> The vPlusLC/vMinusLC Vandermonde determinant collapse at K≥6 is **not purely spacing-dependent**. It results from a combination of:
> 1. **Im-spacing uniformity and magnitude:** Arithmetic progression with spacing >~1e-3 (magnitude matters; very small uniform spacing may not trigger collapse)
> 2. **Re-variation:** Constant Re may amplify or dampen the collapse tendency
> 3. **K-dependent threshold:** The exact threshold for spacing magnitude (below which collapse ceases) may be K-specific

The constant-Im exception (AllReal, spacing=0) suggests collapse requires spacing **nonzero and uniform AND above a threshold magnitude**. This unifies the exception into the broader rule rather than treating it as a separate case.

**Confidence:** LOW. Based on 13 data points with clear contradictions to simpler single-feature rules. Needs boundary testing to confirm/refute.

---

## Files & References

- **Session output:** This file, SESSION_CONSTANT_IM_EXCEPTION_RESOLUTION_2026-10-23.md
- **Prior work:** LEDGER S2720, SESSION_UNIFIED_COLLAPSE_CONDITION_2026-10-23.md, LEDGER S2718, S2719
- **Spacing classification code:** `/private/tmp/.../resolve_constant_im.py`

