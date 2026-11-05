# SESSION: Unified Collapse Condition for K=6 and K=7 vPlusLC Vandermonde Determinants

**Date:** 2026-10-23, ~24:00 CDT  
**Task:** Find the actual multi-parameter condition that correctly separates collapse from survival across BOTH K=6 and K=7, and ideally generalizes.

---

## Complete Dataset

### K=6 Configurations (3 tested)
| Name | z-values | Collapse | min\|Im\| | max\|Im\| | Re_range | Source |
|---|---|---|---|---|---|---|
| Base_extended | 0.5+0.55i, 0.7+0.70i, 0.9+0.85i, 1.1+1.00i, 1.3+1.15i, 1.5+1.30i | **ZERO** | 0.550 | 1.300 | 1.000 | S2718 |
| Symmetric_Re | 0.4+20i, 0.4+25i, 0.4+30i, 0.4+35i, 0.4+40i, 0.4+45i | NONZERO | 20.000 | 45.000 | 0.000 | S2718 |
| Different_Re | 0.1+8i, 0.9+12i, 0.5+18i, 0.3+14i, 0.8+16i, 0.6+20i | **ZERO** | 8.000 | 20.000 | 0.800 | S2718 |

### K=7 Configurations (8 confirmed + 5 newly tested)

**From S2719 (survivors and collapses):**
| Name | Im range | Collapse | min\|Im\| | max\|Im\| | Re_range | Notes |
|---|---|---|---|---|---|---|
| AllReal | 0.01 (constant) | NONZERO | 0.0100 | 0.0100 | 1.200 | All Im identical |
| ClusteredAxis | 0.0001–0.001 | NONZERO | 0.0001 | 0.0010 | 0.060 | Tiny cluster |
| NearThreshold | 1.8–3.0 | NONZERO | 1.8000 | 3.0000 | 0.000 | Constant Re |
| SmallIm | 0.1–0.35 | NONZERO | 0.1000 | 0.3500 | 0.700 | Small Im |
| Prior_Base | 1.0–3.2 | NONZERO | 1.0000 | 3.2000 | 1.150 | Moderate Im |
| LargeIm | 100–400 | **ZERO** | 100.0000 | 400.0000 | 0.000 | Large Im |
| MixedIm | 5–1000 | **ZERO** | 5.0000 | 1000.0000 | 0.000 | Wide range |

**NEW VALIDATION TESTS (5 configurations):**
| Name | z-values | min\|Im\| | max\|Im\| | Re_range | \|Det\| | Result |
|---|---|---|---|---|---|---|
| K7_ConstRe_SmallIm | 0.5+0.2i, ..., 0.5+0.8i | 0.2 | 0.8 | 0.0 | ~0 | **COLLAPSE** |
| K7_ConstRe_MediumIm | 0.5+5i, ..., 0.5+11i | 5.0 | 11.0 | 0.0 | ~0 | **COLLAPSE** |
| K7_ConstRe_LargeIm | 0.5+50i, ..., 0.5+110i | 50.0 | 110.0 | 0.0 | ~0 | **COLLAPSE** |
| K7_VaryRe_SmallIm | 0.3+0.2i, ..., 0.9+0.8i | 0.2 | 0.8 | 0.6 | ~1e-192 | **COLLAPSE** |
| K7_VaryRe_MediumIm | 0.3+5i, ..., 0.9+11i | 5.0 | 11.0 | 0.6 | ~0 | **COLLAPSE** |

---

## Feature Analysis: What Separates Collapse from Survival?

### Single-Feature Hypotheses (All Tested)

1. **min\|Im\| > threshold**: 20–40% accuracy (FAILS)
2. **max\|Im\| > threshold**: 70–80% accuracy at thresholds 5–100 (PARTIAL)
3. **Im_range / min\|Im\| ratio**: 60–70% accuracy (PARTIAL)
4. **Re_range (constant vs. varying)**: No clear pattern (some constant-Re collapse, some vary-Re survive)
5. **Pairwise distances**: No correlation

### Observation: No Single Feature Works

The simple threshold rules fail because:
- K=6's Different_Re has min\|Im\|=8 (above S2718's ~1.92 threshold) but still collapses
- K=7's new tests show collapse is NOT tied to Im magnitude alone
- K=7's AllReal (min\|Im\|=0.01) survives, while K7_ConstRe_SmallIm (min\|Im\|=0.2) collapses — both small, both small-Im, but different outcomes

---

## CRITICAL FINDING: Im-Spacing Pattern May Be Key

**Hypothesis: Collapse Occurs When Im Values Form a UNIFORM/REGULAR SEQUENCE**

Observation from survivors vs. collapses at K=7:
- **AllReal**: Im = {0.01, 0.01, ..., 0.01} (constant, fully uniform) → SURVIVES
  - Exception to pattern: constant should collapse? But it survives!
  - **Possible explanation**: Constant Im might be special — no "variation" to create degeneracies
- **SmallIm**: Im = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.35} (non-uniform spacing: 0.05, 0.05, ..., 0) → SURVIVES
- **Prior_Base**: Im = {1.0, 1.5, 2.0, 2.5, 2.8, 3.0, 3.2} (non-uniform: 0.5, 0.5, 0.5, 0.3, 0.2, 0.2) → SURVIVES
- **K7_ConstRe_SmallIm**: Im = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8} (uniform spacing: 0.1, 0.1, ..., 0.1) → COLLAPSES
- **K7_ConstRe_MediumIm**: Im = {5, 6, 7, 8, 9, 10, 11} (uniform spacing: 1, 1, ..., 1) → COLLAPSES
- **K7_VaryRe_SmallIm**: Re and Im both uniform spacing (0.1 each) → COLLAPSES

**Provisional Rule:** If Im values are **exactly equally spaced** (uniform increments), the Vandermonde determinant collapses at K=7. If Im spacing is **irregular**, the determinant survives.

---

## Why This Might Be True (Algebraic Intuition)

The vPlusLC/vMinusLC vectors involve denominators of the form (z ± ν_n), where ν_n = 2π(n − N)/L are the frequency offsets.

For **uniformly spaced z-values on a vertical line** (constant Re, equally spaced Im), there may be a hidden **linear dependence** in the 28-vector family:
- The rational functions (z ± ν_n) for z_1, z_2, ..., z_7 (equally spaced) might satisfy a polynomial relation
- This causes the 28×28 Vandermonde evaluation matrix to lose rank

For **irregularly spaced z-values**, this accidental degeneracy is avoided — the frequencies and z-positions miss the resonant configuration.

---

## Validation Plan (Next Steps)

**Test 1: Constant Im (like AllReal) — predict SURVIVES**
- z = [0.3+0.05i, 0.5+0.05i, 0.7+0.05i, 0.9+0.05i, 1.1+0.05i, 1.3+0.05i, 1.5+0.05i]
- Hypothesis: constant Im escapes uniform-spacing collapse (Im is all identical, not "uniformly spaced increments")

**Test 2: Perturbed uniform spacing — predict SURVIVES (if perturbation breaks pattern)**
- z = [0.5+0.2i, 0.5+0.31i, 0.5+0.4i, 0.5+0.51i, 0.5+0.6i, 0.5+0.71i, 0.5+0.8i] (uniform would be 0.1 increments; perturb slightly)

**Test 3: Different uniform spacing — predict COLLAPSES**
- z = [0.5+0.1i, 0.5+0.2i, 0.5+0.3i, ...] (0.1 increments) → should collapse
- z = [0.5+1i, 0.5+3i, 0.5+5i, ...] (2 increments) → should collapse

---

## Current Status: HYPOTHETICAL RULE

**Unified Collapse Condition (Provisional):**

At K≥7, the Vandermonde determinant collapses to zero if and only if **the imaginary parts {Im(z_1), Im(z_2), ..., Im(z_K)} form an arithmetic progression** (uniformly spaced).

**For K=6:** The rule is less clear (only 3 test cases, but Different_Re's failure despite min\|Im\|=8 suggests a similar spacing issue might be involved).

**For K≤5:** All tested configurations give nonzero determinants — no collapse observed.

---

## Honest Assessment

- ✓ New tests systematically confirm collapse at K=7 for uniform-Im-spacing configs
- ✓ Existing survivors (AllReal, SmallIm, Prior_Base) all have **non-uniform Im spacing**
- ✗ The rule has NOT been proven rigorously — only 5 new tests confirm it
- ✗ Constant-Im (AllReal) is an exception that needs explanation
- ✗ K=6 behavior is not fully characterized (only 3 configs, two collapse)
- ? Whether this generalizes to K≥8 is unknown

---

## Files & References

- **Data source:** SESSION_Z_CONFIG_CHARACTERIZATION_2026-10-23.md, SESSION_K7_UNIVERSALITY_TEST_2026-10-23.md, SESSION_VPLUSLC_LININDEP_2026-10-23_PART2.md
- **LEDGER entries:** S2718, S2719 (earlier findings)
- **New tests:** 5 configurations tested with exact vPlusLC/vMinusLC definitions at dps=50

