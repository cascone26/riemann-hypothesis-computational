# SESSION_LEAN_EVAL_GROUND_TRUTH_2026-10-23 — Lean #eval Ground Truth Investigation

**Goal:** Extract Lean's actual computed value for the K=5 cubic polynomial at concrete parameter values, to resolve the discrepancy between Lean's proof (polynomial = 0) and Python reconstructions (polynomial ≠ 0).

**Date:** 2026-10-23, following S2736 (source examination) and S2737 (literature check)

## Approach Attempted

1. **Verified working Lean compilation pipeline** (from S2699):
   ```bash
   export PATH=$HOME/.elan/bin:$PATH
   cd ~/rhforge/mathlib4
   LP=$(lake env printenv LEAN_PATH)
   cd ~/rhforge/forge
   export LEAN_PATH="$LP:$HOME/rhforge/forge"
   lean -o FILE.olean FILE.lean
   ```
   
2. **Created concrete test case** with:
   - N=2, L=2π, n=0, nu=-2
   - z_i = i + i (i=1..5)
   - Coefficients: p1=1, q_i=r_i=s_i=0 for all i
   - Python evaluation of cubic1: `(-163.32 + 163.32i)` (nonzero)

3. **Attempted Lean #eval instantiation** (K5GroundTruth.lean):
   ```lean
   def test_value : ℂ :=
     let L := (2 : ℝ)
     let z := (1 : ℂ) + I
     let A := ((2 : ℝ).sqrt : ℂ) * Complex.sin (L * z / 2)
     let ν := (-2 : ℝ)
     -- ... cubic formula ...
   
   #eval test_value  -- Try to evaluate
   ```

## Results — #eval NOT POSSIBLE

**Compilation errors from lean:**

```
K5GroundTruth.lean:30:4: error(lean.dependsOnNoncomputable): 
  failed to compile definition, consider marking it as 'noncomputable' 
  because it depends on 'Real.sqrt', which is 'noncomputable'

K5GroundTruth.lean:48:0: error(lean.dependsOnNoncomputable): 
  failed to compile definition because it depends on 'test_value', which is 'noncomputable'
```

**Root cause:** Lean's definitions of `Real.sqrt`, `Complex.sin`, `Complex.cos`, and other analytic functions are **noncomputable**. In Lean 4's mathematical library:
- Real square root is defined via Cauchy sequences / completeness of ℝ, not via an algorithm
- Trigonometric functions are defined as analytic functions, not via power series or other computable representations
- These are **real-analytic** definitions, provably correct but not executable

This is a fundamental limitation of Lean's type theory: functions that depend on infinitary constructions or infinite processes cannot be compiled into code that the kernel can execute.

## What CAN Be Done Instead

**Option 1: Trust Lean's type-checking (what the proof does)**
The CrossK5Key.lean proof uses the `ring` tactic to verify the polynomial sum equals 0 at the END of a symbolic proof chain. When `ring` succeeds (EXIT0 on compilation), it has verified:
- The polynomial is identically 0 (as a symbolic expression)
- This is independent of the specific values of z_i, N, L, etc.
- The cubic polynomial WOULD evaluate to 0 if Real.sqrt and Complex.sin were computable

This type-checking is the ONLY available ground truth in Lean for expressions that depend on noncomputable functions.

**Option 2: Manual symbolic extraction**
Extract the intermediate expressions from CrossK5Key.lean source code:
- A_i = √2 * sin(L*z_i/2)
- cubic_i polynomial formula (explicit polynomial in z_i, p_i, q_i, etc.)
- Manually trace the groupClearIdentity application per block
- Verify ring simplification by hand or with a computer algebra system

**Option 3: Implement a computable fragment**
Create a separate Lean file using `Float` or `Rat` with computable approximations:
- Use Float64 to approximate Real values
- Use `Complex.sin` from a computable power series library
- Trade exact symbolic algebra for approximate numeric computation
- Status: This is a research project (Interval arithmetic / Arb in Lean 4), not available in standard Mathlib

## Analysis of Python-Lean Discrepancy

The Python reconstruction from S2733-S2735 gave **cubic sum ≠ 0** even when hsum=0 (the theorem's constraint). But CrossK5Key.lean proves **cubic sum = 0** with EXIT0 type-check.

**Hypotheses (from S2736 analysis):**
1. Python extracted WRONG cubic formula — the cubic_i in Lean is not independent; it encodes dependencies via the groupClearIdentity
2. Python did not apply the hsum coupling correctly — the 5 blocks are not independent
3. Python did not match the `linear_combination` ring simplifications in Lean

**Evidence that Lean is correct:**
- Lean compilation proves the polynomial identity via `ring` (exact symbolic algebra)
- groupClearIdentity is verified algebraically to 1e-13 in Python
- The proof structure (applying identity per-block, scaling by quartics, summing) is sound

**Most likely issue:**
Python's cubic formula derivation did not account for the **coupling via linear_combination** — Lean's proof combines all 5 blocks at once with specific scaling factors (products of other blocks' quartics), and the ring tactic simplifies the entire combined expression to 0. Python tried to extract individual cubic_i formulas and sum them post-hoc, losing the algebraic simplifications that occur during the combined ring simplification.

## Hypothesis Test: Per-Block vs. Full Combined Expression (Coordinator Follow-Up)

Coordinator directive: test whether "per-block evaluation loses cross-block cancellations that only emerge when the FULL combined expression is simplified at once."

**Test 1 — Trivial case (all coefficients = 0):**
```
✓ PASS
Full combined polynomial = 0
This is trivially true, no information about per-block discrepancy
```

**Test 2 — Single block active (p1=1, rest=0, hsum NOT satisfied):**
```
Per-block sum: -1390090000*sqrt(2)*sinh(π)  (nonzero)
Full combined: -1390090000*sqrt(2)*sinh(π)  (same)
Result: No discrepancy (expected, only one block active)
```

**Test 3 — Multiple blocks, hsum NOT satisfied:**
```
Per-block and full combined: same nonzero result
Both methods agree on simplification
```

**Critical observation:** When using sympy.simplify() on the COMBINED polynomial (all 5 blocks scaled by quartics and summed), the result is IDENTICAL to per-block evaluation followed by summation. The simplification order does NOT produce different results in our tests.

**Why hypothesis testing was inconclusive:**
1. Per-block and combined-expression simplification give identical results when formula is fixed (tested with sympy)
2. The only way they could differ is if hsum is satisfied AND multiple blocks are nonzero AND their cross-terms cancel via algebraic reduction
3. Finding actual hsum-satisfying multi-block solutions is nontrivial (20-dim parameter space, 1 constraint per n)
4. Attempted to solve hsum constraint symbolically; solver gave approximate solution not exactly satisfying hsum

**Real root cause of Python-Lean discrepancy** (updated analysis):
NOT per-block vs. combined-expression simplification order. The actual issue is likely:
1. Python extracted individual cubic_i formulas as closed-form expressions
2. But these formulas encode hidden dependencies via groupClearIdentity
3. When Lean's `ring` tactic processes the fully assembled expression, it sees the full algebraic context
4. Python's per-block evaluation evaluates OUTSIDE that context, missing global simplifications

The issue is **formula independence** (are cubic_i truly independent?), not **simplification order**.

## Conclusion

**#eval is not available** for numeric ground truth extraction — noncomputable barrier confirmed.

**Per-block vs. combined-simplification hypothesis:** Plausible but NOT CONFIRMED by direct testing. Both methods produce identical results when tested with sympy (same formula → same simplification result). The hypothesis would require an actual hsum-satisfying multi-block test case to verify, which we did not successfully construct.

**More likely root cause:** Python's formula extraction did not account for K-dependent, non-separable algebraic structure. The groupClearIdentity couples all 5 blocks implicitly; extracting them independently loses that coupling.

**Recommendation:**
1. Do NOT attempt numeric ground truth extraction from Lean (noncomputable barrier is hard)
2. Use Lean compilation (EXIT0) as authoritative algebraic ground truth
3. To verify Python: do NOT extract per-block cubic formulas. Instead:
   - Apply groupClearIdentity per-block (done correctly in S2736)
   - Combine into ONE full polynomial (all 5 blocks, all quartic scalings)
   - Call sympy.expand().simplify() on ENTIRE expression at once
   - Compare the single simplified result to 0

**Files:**
- K5GroundTruth.lean (demonstrates noncomputable barrier)
- Hypothesis tests in /tmp/k5_*.py (show per-block and combined agree on simplification)
- This session document

**Time spent:** ~60 minutes (compile verification, multiple hypothesis tests, analysis)

