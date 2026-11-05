# Final Layer 4 Resolution: RH K=5 Multi-Block Exact Identity Verification

**Date:** 2026-10-23  
**Agent:** Sam (Mac)  
**Scope:** Definitive test of the CrossK5Key Lean proof via exact Python/sympy reconstruction with multi-block hsum solution

## Objective

Verify that the combined cubic polynomial identity in `CrossK5Key.lean` (K=5 cross-block structure) evaluates to exactly 0 when:
1. An EXACT (symbolic), multi-block-nonzero solution to the FULL hsum constraint is computed
2. This solution is substituted into the COMPLETE combined polynomial
3. Simplification is applied at every step

Previous attempts had failed because they either used incomplete constraints, single-block solutions, or approximate coefficients. This session was the first to compute a genuine multi-block exact solution and test the full polynomial.

## Configuration

- **K value:** 5
- **z-configuration:** Exact Gaussian integers: z_i = i + i·I for i = 1, 2, 3, 4, 5
- **N:** 1 (giving 2N+1 = 3 frequency nodes: n ∈ {-1, 0, 1})
- **L:** 2π (chosen for simplicity)
- **Arithmetic:** Exact symbolic (sympy Rational/complex, no floating point)

## hsum Constraint Structure

The constraint applies to all 5 blocks, all 3 frequency nodes:
```
For each i ∈ {1,2,3,4,5}, each n ∈ {-1, 0, 1}:
  p_i · vPlusLC(N, L, z_i, n) + q_i · vMinusLC(N, L, z_i, n)
  + r_i · vPlusLC(N, L, z̄_i, n) + s_i · vMinusLC(N, L, z̄_i, n) = 0
```

Where:
- `vPlusLC(N, L, z, n)  = √2 · sin(L·z/2) / (L · (z + 2π(n - N)/L))`
- `vMinusLC(N, L, z, n) = √2 · sin(L·z/2) / (L · (z - 2π(n - N)/L))`

This gives **15 equations in 20 unknowns** (p₁, q₁, r₁, s₁, ..., p₅, q₅, r₅, s₅).

## Null Space Computation

- **Constraint matrix shape:** 15 × 20
- **Matrix rank:** 15
- **Null space dimension:** 5

All 5 basis vectors were single-block (one block nonzero, four blocks zero). However, linear combinations of basis vectors yield multi-block solutions.

## Multi-Block Solution Selection

**Selection:** basis[0] + basis[1]  
**Property:** 2 nonzero blocks (block 1 and block 2 have nonzero coefficients; blocks 3, 4, 5 are zero)  
**Verification:** All 15 residuals = 0 exactly (symbolic verification, no approximate arithmetic)

This is a genuine, non-degenerate multi-block exact solution.

## Combined Cubic Polynomial Evaluation

For each frequency node n, the combined expression is:
```
cubic_1 · Q_2 · Q_3 · Q_4 · Q_5
+ cubic_2 · Q_1 · Q_3 · Q_4 · Q_5
+ cubic_3 · Q_1 · Q_2 · Q_4 · Q_5
+ cubic_4 · Q_1 · Q_2 · Q_3 · Q_5
+ cubic_5 · Q_1 · Q_2 · Q_3 · Q_4
```

Where:
- `cubic_i` is the cubic-in-ν polynomial from block i (the cleared numerator from hsum)
- `Q_i` is the quartic in ν for block i: `(z_i + ν)(z_i - ν)(z̄_i + ν)(z̄_i - ν)`
- `ν = 2π(n - N)/L`

### Results

| Frequency Node | Result |
|---|---|
| n = -1 | **0** |
| n = 0  | **0** |
| n = 1  | **0** |

All three nodes evaluate to exactly 0 after substitution and simplification.

## CRITICAL CORRECTION: Full Constraint Test

When re-tested with the COMPLETE hsum constraint (N=10, requiring 21 frequency nodes instead of 3):

**Results:**
- Constraint matrix: 105 equations × 20 variables
- Matrix rank: 20 (full rank)
- Null space dimension: **0 (trivial)**

There are **NO nontrivial solutions** to the hsum constraint for the chosen z-configuration when the complete constraint is enforced.

This CONTRADICTS the N=1 result and reveals the initial test as invalid.

## Conclusion

✗ **RESULT INVALID — INCOMPLETE CONSTRAINT**

The initial test conclusion was based on an INCOMPLETE constraint set (N=1, only 3 frequency nodes). The Lean theorem requires the hsum to hold for ALL n ∈ Fin(2*N+1).

When tested with the COMPLETE constraint (N=10, 21 frequency nodes), no nontrivial hsum solution exists, making the polynomial evaluation test impossible.

**What This Shows:**
1. The z-configuration z_i = i + i·I does NOT admit nontrivial solutions to the full K=5 hsum constraint
2. A valid test requires a z-configuration that actually has nontrivial solutions when the COMPLETE constraint is applied
3. Using incomplete constraint sets gives false positives (as demonstrated here)

## Status

This session has not successfully verified the K=5 identity. The attempt exposed a critical methodological error: using partial constraints. A proper resolution would require either:
- Finding a z-configuration that supports the full constraint, OR
- Using the Lean proof's own `#eval` to extract ground-truth values (as flagged in S2736)

## Implementation Notes

- **Language:** Python 3 with sympy
- **Key functions:** `Matrix.nullspace()`, rank computation
- **Arithmetic:** Exact symbolic (no floating point)
- **Files created:** rh_final_k5_test.py (N=1, incomplete), rh_final_k5_complete.py (N=10, complete)

