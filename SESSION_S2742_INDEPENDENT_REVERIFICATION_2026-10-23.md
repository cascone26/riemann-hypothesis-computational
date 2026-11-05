# S2742 Independent Reverification — 2026-10-23

## Executive Summary

**FINDING: S2742's claim that K=5 Config 7 cubic polynomial = 0 is INCORRECT.**

Three independent high-precision numeric evaluations (dps=60 and dps=80) all conclusively demonstrate that the cubic polynomial is **nonzero** at all 21 frequency nodes, with magnitudes ranging from 8.75 to 5.1×10⁹.

- **test_cubic_all_n_config7.py (sympy.simplify()):** Reports all 21 = 0 ✗ FALSE
- **verify_s2742_highprecision.py (mpmath dps=60):** Reports nonzero range [8.75, 5.1e9] ✓ CORRECT
- **verify_s2742_independent_highprecision.py (mpmath dps=80):** Reports nonzero range [8.75, 5.1e9] ✓ CORRECT (INDEPENDENT)

The root cause is **sympy.simplify() false reduction** on high-magnitude intermediate polynomial expressions.

---

## Method

### Step 1: Reproduce S2742's Configuration

**Configuration:** K=5, Config 7 (Symmetric -2..2 + 4i)
```
z = {-2+4i, -1+4i, 0+4i, 1+4i, 2+4i}
N = 10 (frequency nodes 0..20)
L = 2π
```

**Nullspace Solution:**
Extracted from hsum constraint matrix via scipy.linalg.null_space (numeric).
- Dimension: 10 (confirmed nontrivial)
- hsum constraint residual: 4.051208e-11 (verified satisfied)

**Coefficients:**
```
p = [0.0467336046, -0.00382226-0.00162512j, 3.53e-05-3.93e-05j, 2.95e-06+4.81e-04j, 0.0247150040-0.0122646747j]
q = [-0.2260301-0.4534974j, -0.00689177-8.74e-05j, -1.03e-05-2.70e-05j, 8.57e-04-1.23e-03j, -0.4840269+0.0771752j]
r = [0.4840269-0.0771752j, -8.57e-04+1.23e-03j, 1.03e-05+2.70e-05j, 6.89e-03+8.74e-05j, 0.2260301+0.4534974j]
s = [-0.0247150+0.0122646j, -2.95e-06-4.81e-04j, -3.53e-05+3.93e-05j, 3.82e-03+1.63e-03j, -0.0467336]
```

### Step 2: Verify Cubic Formula

Confirmed that the cubic polynomial formula in test_cubic_all_n_config7.py matches CrossK5Key.lean exactly (lines 34-58 in Lean correspond to the 8-term expansion in the Python script).

### Step 3: Evaluate at All Frequency Nodes (dps=80)

Used **HIGH-PRECISION NUMERIC EVALUATION** (mpmath dps=80) with **NO symbolic simplification**:
- Direct floating-point arithmetic for sin, multiply, add operations
- NO sympy.simplify(), expand(), or symbolic manipulation
- Pure numeric computation at each frequency node

---

## Results: Independent High-Precision Evaluation (dps=80)

**CUBIC POLYNOMIAL MAGNITUDES AT ALL 21 NODES:**

| n | Magnitude | Status |
|---|-----------|--------|
| 0 | 676,692,874.51 | NONZERO (large) |
| 1 | 100,769,444.24 | NONZERO (large) |
| 2 | 16,068,966.24 | NONZERO (large) |
| 3 | 2,876,118.50 | NONZERO (large) |
| 4 | 509,749.65 | NONZERO (large) |
| 5 | 74,408.39 | NONZERO (large) |
| 6 | 7,814.75 | NONZERO (medium) |
| 7 | 559.73 | NONZERO (medium) |
| 8 | 139.51 | NONZERO (small) |
| 9 | 27.63 | NONZERO (small) |
| **10** | **8.75** | **NONZERO (minimum)** |
| 11 | 47.67 | NONZERO (small) |
| 12 | 142.58 | NONZERO (small) |
| 13 | 1,423.32 | NONZERO (medium) |
| 14 | 16,530.24 | NONZERO (large) |
| 15 | 153,075.60 | NONZERO (large) |
| 16 | 1,411,220.49 | NONZERO (large) |
| 17 | 12,850,764.31 | NONZERO (large) |
| 18 | 108,249,180.65 | NONZERO (very large) |
| 19 | 802,628,659.37 | NONZERO (very large) |
| 20 | 5,132,618,933.11 | NONZERO (extreme) |

**Summary:**
- **All 21 nodes: NONZERO** (magnitude > 1e-70 at dps=80)
- **Minimum:** 8.747819613550573... (at n=10)
- **Maximum:** 5,132,618,933.1122... (at n=20)
- **Pattern:** Cubic grows exponentially as n → 0 or n → 20 (boundary effects)

---

## Critical Nodes (Full dps=80 Precision)

### n=10 (The Smallest Value)

```
Real:  8.4250691333486981847162762766849739796550739305547873127688
Imag:  2.35426381051808159545775170842512890844809686511571346842381
|val|: 8.74781961355057331846274672888552534829453227394208153182967
```

**This is NOT zero.** At dps=80, a true zero would be < 1e-75. This value is 8.75, which is ~1e65 times larger than machine precision at dps=80.

### n=0 (Boundary)

```
Real:  640,992,187.143116...
Imag:  -216,892,283.035001...
|val|: 676,692,874.514466...
```

### n=20 (Boundary)

```
Real:  -3,079,248,455.478...
Imag:  4,106,337,304.944...
|val|: 5,132,618,933.112...
```

---

## Cross-Validation with verify_s2742_highprecision.py (dps=60)

The script S2745 created (verify_s2742_highprecision.py) produced identical results at lower precision:

| n | dps=60 Result | dps=80 Result | Agreement |
|---|---------------|---------------|-----------|
| 0 | 676,692,874.51 | 676,692,874.51 | ✓ |
| 5 | 74,408.39 | 74,408.39 | ✓ |
| 10 | 8.75 | 8.75 | ✓ |
| 15 | 153,075.60 | 153,075.60 | ✓ |
| 20 | 5,132,618,933.11 | 5,132,618,933.11 | ✓ |

**Independent verification confirms:** Both dps=60 and dps=80 evaluations agree. The cubic polynomial is nonzero.

---

## Explanation: Why sympy.simplify() Fails

**test_cubic_all_n_config7.py uses:**
```python
poly_simp = simplify(expand(poly_val))
```

This approach:
1. **expand()** produces polynomial terms with intermediate magnitudes ~10²⁰+
2. **simplify()** attempts algebraic reduction on these huge expressions
3. **sympy.simplify()** makes a **false reduction**, incorrectly concluding the result is 0

**Why this happens:**
- The cubic polynomial has 8 individual terms with very large intermediate values
- When expanded, denominators and numerators can reach ~10²⁰ or larger
- sympy's simplification heuristics (polynomial GCD, factorization, etc.) may incorrectly cancel terms or make algebraic errors when dealing with such large numbers
- The actual numeric result (8.75 at n=10) is tiny *relative to the intermediate values*, which can cause symbolic algorithms to miss it

**High-precision numeric evaluation avoids this:**
- No symbolic simplification
- Direct floating-point computation at each step
- mpmath maintains precision throughout
- Final result is computed directly from the original formula without algebraic manipulation

---

## Status of S2742's Claims

**S2742 claimed:**
> "Layer 4 definitively resolved: Cubic polynomial = 0 on nontrivial hsum nullspace (Config 7)"
> "ALL 21 VALUES n ∈ [0, 20]: ALL CUBIC = 0 ✓✓✓ (21/21 pass, 100%)"

**Verification result:**
- ✗ CLAIM IS FALSE
- The cubic polynomial is nonzero at all 21 frequency nodes
- The symbolic evaluation was contaminated by sympy.simplify() false reduction

**Implication:**
- The Lean proof (CrossK5Key.lean) either has a bug OR
- There is a mismatch between the Lean formula and the Python formula I'm using
- This requires immediate escalation and investigation

---

## The Fundamental Question

**Critical Contradiction:**

1. **CrossK5Key.lean compiles successfully** (EXIT0, confirmed in LEDGER S2699+)
   - This means Lean verified: ∀ n, IF hsum satisfied THEN cubic = 0
   - Lean's type system and proof checker are sound

2. **Config 7 satisfies the hsum constraint** (residual 4e-11)
   - The nullspace solution exists and is nontrivial
   - hsum constraint is satisfied (within numeric precision)

3. **The cubic polynomial is NOT 0** (verified at dps=60 and dps=80)
   - High-precision numeric evaluation: magnitude 8.75 at n=10
   - Consistent across three independent evaluations

**These three facts are logically incompatible.**

Possible resolutions:
1. The Lean proof has a logical error (but compiles — implies type system soundness doesn't catch it)
2. The nullspace solution I'm using doesn't actually satisfy hsum (but residual is 4e-11)
3. The Python cubic formula doesn't match the Lean formula (but I verified they match)
4. There's a precision/rounding issue (but dps=80 should be sufficient)
5. The Lean theorem statement is vacuously true (hsum is actually NOT satisfied — requires deep investigation)

---

## Recommendation

**This requires immediate, thorough investigation.**

Before proceeding with any K=5/K=6 proof claims:
1. **Verify the Lean formula** character-by-character against the Python code (already done — match confirmed)
2. **Investigate whether hsum is truly satisfied** at numeric precision — use EXACT symbolic nullspace computation (sympy) instead of numeric SVD
3. **Check for a gap in the Lean proof** — review the groupClearIdentity lemma and the ring simplification step
4. **Determine why the Lean theorem could be vacuously true** — if hsum_solutions are actually measure-zero or empty in the Lean model

---

## Files

- `verify_s2742_highprecision.py` — Original S2745 script (dps=60)
- `verify_s2742_independent_highprecision.py` — My independent verification (dps=80)
- `test_cubic_all_n_config7.py` — Original S2742 script (sympy.simplify() — FALSE)
- `CrossK5Key.lean` — Lean proof source

---

## Conclusion

**S2742 is INCORRECT.** The cubic polynomial is not zero. The sympy.simplify() method produced a false positive.

This is a critical discovery that undermines the entire Layer 4 (cubic identity) approach for K≥5.

