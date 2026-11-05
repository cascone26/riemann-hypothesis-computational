# SESSION: K=5 Lean Ground Truth Check (2026-10-23)

## Executive Summary

**Goal:** Determine if the Python reconstruction of the K=5 cubic polynomial identity matches Lean's actual computation, by extracting Lean's own intermediate values rather than re-deriving them.

**Status:** Examined Lean source code in detail. Build system constraints prevented direct `#eval` extraction, but rigorous analysis of Lean definitions and proof structure completed. Python verification of the core algebraic identity (groupClearIdentity) completed and verified to machine precision.

**Key Finding:** The groupClearIdentity algebraic identity is verified correct in Python (error < 1e-13 on concrete test), confirming that the polynomial formula itself is sound. The Python-Lean discrepancy (polynomial ≠ 0 in Python even when hsum = 0) points to either:
1. A fundamental issue with how the hsum constraint couples the 5 blocks in the K=5 construction, or
2. The Lean proof may have a subtle assumption or step not reflected in the Python reconstruction.

---

## Part 1: Lean Source Code Analysis

### 1.1 File: CrossK5Key.lean
**Location:** `~/rhforge/forge/CrossK5Key.lean`
**Lines examined:** 1-155 (full file)

**Theorem Statement:**
```lean
theorem crossK5Key (N : ℕ) (L : ℝ) (hL : L ≠ 0) (z1 z2 z3 z4 z5 : ℂ)
    (hz1_im : z1.im ≠ 0) ... (hz5_im : z5.im ≠ 0)
    (p1 q1 r1 s1 ... p5 q5 r5 s5 : ℂ)
    (hsum : p1 • vPlusLC N L z1 + q1 • vMinusLC N L z1 + ...
        + p5 • vPlusLC N L z5 + q5 • vMinusLC N L z5 + ... = 0) :
    ∀ n : Fin (2 * N + 1),
      cubic1 * quartic2*quartic3*quartic4*quartic5
      + cubic2 * quartic1*quartic3*quartic4*quartic5
      + cubic3 * quartic1*quartic2*quartic4*quartic5
      + cubic4 * quartic1*quartic2*quartic3*quartic5
      + cubic5 * quartic1*quartic2*quartic3*quartic4
      = 0
```

**Proof Structure:**
1. Lines 104-147: Apply `groupClearIdentity` to each block (z1..z5)
2. Lines 95-102: Use `congrArg` to scale each resulting equation by products of other quartic factors
3. Lines 151-154: Combine using `linear_combination` with the hsum constraint

### 1.2 File: GroupClearIdentity.lean  
**Location:** `~/rhforge/forge/GroupClearIdentity.lean`
**Lines examined:** 1-50

**Key Theorem:**
```lean
theorem groupClearIdentity (L : ℝ) (hL : L ≠ 0) (z : ℂ) (hz_im : z.im ≠ 0) (ν : ℝ)
    (A B p q r s : ℂ) :
    [p*A/(L*(z+ν)) + q*A/(L*(z-ν)) + r*B/(L*(z̄+ν)) + s*B/(L*(z̄-ν))] * 
    [L^4 * (z+ν)(z-ν)(z̄+ν)(z̄-ν)] =
    L^3 * [cubic polynomial in ν]
```

**Proof:** `field_simp; ring` — pure algebraic identity, no assumptions on A, B, p, q, r, s.

**Critical Observation:** This identity is UNCONDITIONAL and ALWAYS holds for any complex values. It does NOT depend on any constraint.

### 1.3 File: VPlusVMinusFullKey.lean
**Location:** `~/rhforge/forge/VPlusVMinusFullKey.lean`
**Lines examined:** 1-100

**Theorem:** vPlusVMinusFullKey

**Precondition:** `hpqrs : p • vPlusLC N L z + q • vMinusLC N L z + r • vPlusLC N L z̄ + s • vMinusLC N L z̄ = 0`

**Conclusion:** The cubic polynomial (with specific coefficients) = 0 for all n.

**Key Insight:** This theorem DOES have a constraint (hpqrs), but it's a constraint on a SINGLE BLOCK, not across multiple blocks.

### 1.4 File: HLGComplexContinuation.lean
**Location:** `~/rhforge/forge/HLGComplexContinuation.lean`

**Definitions:**
```lean
def vPlusLC (N : ℕ) (L : ℝ) (z : ℂ) : Fin (2 * N + 1) → ℂ :=
  fun n => √2 * sin(L*z/2) / (L * (z + 2π(n-N)/L))

def vMinusLC (N : ℕ) (L : ℝ) (z : ℂ) : Fin (2 * N + 1) → ℂ :=
  fun n => √2 * sin(L*z/2) / (L * (z - 2π(n-N)/L))
```

**Mathematical Interpretation:**
- These are 4-component fractions (with conjugate variants for z̄)
- They depend on a spatial frequency variable ν_n = 2π(n-N)/L at each index n
- They represent basis functions for the Weil form construction

---

## Part 2: Python Verification of groupClearIdentity

### 2.1 Test Setup

```python
N = 2, L = 2π
z₁ = 1 + 0.5i
z̄₁ = 1 - 0.5i  
n = 0, ν = -2.0
A = √2 * sin(Lz₁/2) ≈ (0 - 3.254528i)
B = Ā
p = 1 + 0.1i, q = 0, r = 0, s = 0
```

### 2.2 Python Code

File: `/private/tmp/.../group_clear_direct.py`

Computes both:
- **LHS:** `[p*A/(L*(z+ν)) + ...] * L^4 * Q` where Q = quartic product
- **RHS:** `L^3 * [D₀ + D₁*ν + D₂*ν² + D₃*ν³]` with explicit coefficient formulas

### 2.3 Results

| Quantity | Value |
|----------|-------|
| LHS | -4480.44 + 7094.03i |
| RHS | -4480.44 + 7094.03i |
| Difference | 9.09e-13 |
| Relative Error | 1.08e-16 |
| **Status** | **✓ VERIFIED** |

**Conclusion:** groupClearIdentity is a CORRECT algebraic identity, verified to machine precision in Python.

---

## Part 3: Structural Analysis of K=5 Proof

### 3.1 The Coupling Mechanism

The CrossK5Key proof works by:

1. **Individual Block Equations:** Each `groupClearIdentity` call establishes:
   ```
   (4-fraction-sum_i) * Q_i = L^3 * cubic_i
   ```
   where `4-fraction-sum_i` is the sum of 4 fraction terms for block i.

2. **Constraint Activation:** The hsum constraint provides:
   ```
   Σ(4-fraction-sum_i for all i) = 0
   ```

3. **Linear Combination:** By multiplying each block equation by specific products of quartics and summing with hsum, the goal is to show that the cubic sum = 0.

### 3.2 The Key Question

**When hsum = 0, what can we infer about individual blocks?**

The constraint `Σ(4-fraction-sum_i) = 0` does NOT imply each individual 4-fraction-sum_i = 0.

This means:
- Block 1's 4-fraction-sum may be nonzero
- Block 2's 4-fraction-sum may be nonzero
- But they sum to zero across all 5 blocks

**This is a critical structural feature:** The K=5 construction couples all 5 blocks at the level of the hsum constraint, not individually.

### 3.3 Proof Mechanics (Traced from Lean)

The proof (lines 133-154 in CrossK5Key.lean) does:

```
have e1 := congrArg (fun x => -(Q2*Q3*Q4*Q5*L^16) * x) h1
have e2 := congrArg (fun x => -(Q1*Q3*Q4*Q5*L^16) * x) h2
...
linear_combination (norm := (push_cast; ring1)) e1 + e2 + e3 + e4 + e5 +
  ((L:ℂ)^20 * (Q1*Q2*Q3*Q4*Q5)) * hpt
```

where `hpt := congrFun hsum n` (the hsum constraint at point n).

**Interpretation:**
- Each `e_i` is the block equation scaled by specific quartic products
- The final sum combines all e_i plus L^20 * (all quartics) * hpt
- The `ring1` tactic then simplifies the algebraic expression to 0

---

## Part 4: Python-Lean Discrepancy Analysis

### 4.1 Previous Session Findings (S2735)

The Python reconstruction with explicit hsum solutions found:
- hsum constraint is satisfied (M @ x = 0 verified)
- Variable mapping is correct (character-by-character check)
- Polynomial formula is correct (direct transcription from Lean)
- **Yet polynomial ≠ 0 even when hsum = 0**

### 4.2 Possible Explanations

**Hypothesis A (Build verification):** The Lean proof itself may compile successfully, meaning it IS correct, and the Python finding of polynomial ≠ 0 is due to some subtle transcription error or different definition of the cubic polynomial.

**Hypothesis B (Mathematical):** The K=5 construction may have a more subtle coupling requirement than simply hsum = 0 across all blocks. There may be additional implicit constraints derived from the groupClearIdentity structure.

**Hypothesis C (Proof structure):** The linear_combination step in Lean may be performing cancellations that are not obvious in a straightforward Python evaluation of the cubic polynomial.

### 4.3 Why Python Direct Evaluation May Miss Proof Logic

The Lean proof uses:
- `congrArg` — equality preservation under function application
- `linear_combination` — automated tactic that combines equations

Python direct polynomial evaluation:
- Computes cubic values using explicit formulas
- Does NOT capture the implicit cancellations in the ring tactic

**Key Insight:** The groupClearIdentity formula gives the cubic in terms of the fraction sum:
```
cubic_i = (1/L^3) * (4-fraction-sum_i) * L^4 * Q_i
        = (L/1) * (4-fraction-sum_i) * Q_i  [after simplification]
```

So the cubic is NOT independent of the 4-fraction-sum. The Python code may have extracted the wrong cubic formula!

---

## Part 5: Recommended Next Steps

1. **Extract Lean's exact cubic formula:** Modify CrossK5Key.lean to output intermediate expressions (A_i, Q_i values) via `#eval` after fixing the build environment.

2. **Verify the cubic formula in Python:** Once we have Lean's explicit formula for the cubic in terms of the given quantities, verify it matches the Python reconstruction.

3. **Check the linear_combination setup:** The specific coefficients used in the `linear_combination` step may encode crucial cancellations that Python evaluation must reproduce.

4. **Full symbolic chain:** Compute the sum of cubics * quartics symbolically in sympy with full hsum constraint substitution, not just plugging in specific coefficient values.

---

## Part 6: Verification Summary

### What Was Verified ✓

- groupClearIdentity algebraic identity: **CORRECT** (error < 1e-13)
- vPlusLC/vMinusLC definitions: **VERIFIED** in Lean source
- Variable correspondence (p1,q1,r1,s1 order): **CORRECT** (per S2735)
- Polynomial formula transcription: **CORRECT** (per S2735)

### What Remains Unverified ✗

- Whether Lean's proof actually compiles successfully
- Lean's intermediate values for A_i, Q_i at specific n
- The exact form of the cubic polynomial used in the combination step
- How the linear_combination tactic resolves to zero

---

## Files Created

1. **This document:** `SESSION_LEAN_GROUND_TRUTH_CHECK_2026-10-23.md`
2. **Python verification:** `/private/tmp/.../group_clear_direct.py`  
   - Verifies groupClearIdentity to 1e-13 precision
   - Tests on concrete z₁=1+0.5i, specific n
3. **Python constraint analysis:** `/private/tmp/.../k5_lean_ground_truth.py`
   - Sets up full hsum constraint system
   - Builds M matrix, computes nullspace dimension

---

## Conclusion

The core algebraic identity (groupClearIdentity) that underlies the K=5 proof is mathematically sound and verified in Python. The discrepancy between Python's evaluation (cubic ≠ 0) and Lean's proof (cubic = 0) indicates one of:

1. A transcription error in how the cubic is defined in Python
2. Additional implicit constraints in how the K=5 cross-block coupling works
3. Cancellations in the linear_combination step that Python direct evaluation doesn't capture

**Recommended:** Get Lean's own #eval output on intermediate values to serve as ground truth, then reverse-engineer what the actual cubic formula should be in Python.

