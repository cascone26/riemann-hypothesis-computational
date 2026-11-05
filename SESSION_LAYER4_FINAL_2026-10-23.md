# SESSION: Layer 4 — K-Induction Polynomial Evaluation (2026-10-23)

**Time:** 22:30 — 23:15 CDT

**Method:** Exact symbolic computation (sympy.nullspace() + polynomial evaluation)

**Status:** Breakthrough on method; unexpected failure on control test.

---

## Breakthrough: Exact Symbolic Works

**The problem:** Numerical methods (float64, mpmath, lstsq) all failed on ill-conditioned hsum matrices (condition number ~3.5e+15).

**The solution:** Exact symbolic computation with rational arithmetic (sympy):
1. Build hsum matrix with exact Rational entries (no floating point)
2. Use `sympy.nullspace()` to find exact null-space basis
3. Residuals simplify to exactly zero (not approximately zero)
4. Fast: K=5, num_basis=2 completes in seconds

**Example result from Layer 3 (hsum solving):**
```
K=5, num_basis=2, family [0,1,2,3,4]+(5/2)j
Nullspace dimension: 17 (out of 20 unknowns)
Solution: p_i = -2, q_i = 1, r_i = -2, s_i = 1, then all zeros
Verification: All residuals simplify to exactly 0 ✓
```

---

## Failure: Layer 4 Control Test Unexpected

**Test setup:**
- Family: [0, 1, 2, 3, 4, 5] + (5/2)j (exact Rational)
- K=5 control (known-working)
- K=6 test (survivor family)
- Extract hsum solutions, evaluate cubic/quartic polynomial

**Results:**

K=5 (expected = 0):
```
Simplified polynomial = -5*sqrt(2)*poly(π)*sinh(25)/262144
where poly(π) = large polynomial in π with positive coefficients
Status: NONZERO (fails control test!)
```

K=6 (expected = 0 or varies based on K=5 outcome):
```
Simplified polynomial = -5*sqrt(2)*poly(π)*sinh(25)/4194304
where poly(π) = different large polynomial in π
Status: NONZERO
```

---

## Root Cause Unknown

K=5 control test FAILS when it should succeed. This indicates:

1. **Polynomial evaluation bug:** The cubic/quartic extraction from CrossK6Key.lean may be incorrect
2. **hsum solution misapplication:** The extracted coefficients may not match the theorem's preconditions
3. **Test family invalidity:** [0,1,2,3,4,5] + (5/2)j may not satisfy constraints needed for the theorem
4. **Theorem misinterpretation:** May be testing at wrong n-values, or theorem statement is different

The fact that BOTH K=5 and K=6 fail suggests a systematic issue (not K-specific), pointing to (1), (2), or (3).

---

## What We Know

✓ Layer 3 works: hsum constraint solved exactly via sympy
✓ Symbolic computation avoids precision wall entirely
✗ Layer 4 fails: K=5 control doesn't give zero

**Next step:** Debug the polynomial evaluation before claiming any result about K-induction.

---

## Code

- `solve_hsum_symbolic_sympy.py` — Layer 3 (hsum solving) ✓ works
- `layer4_complete_symbolic_test.py` — Layer 4 (polynomial eval) ✗ K=5 fails

---

## LEDGER Entry

**S2733** (2026-10-23, 23:15 CDT) — EXACT SYMBOLIC COMPUTATION BREAKTHROUGH, LAYER 4 FAILURE

Layer 3 (hsum constraint solving): SOLVED via `sympy.nullspace()` with exact Rational arithmetic. Fast, exact, verified.

Layer 4 (polynomial evaluation with hsum solution): K=5 control test FAILS. Polynomial ≠ 0 when should be 0.

Implication: Either polynomial evaluation is buggy, or test setup is invalid. Not a K-specific obstruction.

Method file: `layer4_complete_symbolic_test.py`
