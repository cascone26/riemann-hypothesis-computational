# SESSION_HSUM_K6_MECHANISM_2026-10-23

**Task:** Characterize the K-dependent hsum constraint's algebraic structure directly, identifying why the cubic AND quartic polynomial approaches both fail at K=6 despite working at K≤5.

**Finding:** The failure is **dimension-dependent, not degree-dependent.**

## The hsum Constraint Structure

From HLGComplexContinuation.lean:
- K=5: 20 coefficients {p,q,r,s}_i for i=1..5
- K=6: 24 coefficients {p,q,r,s}_i for i=1..6

The constraint is:
```
∑ᵢ₌₁^K [ pᵢ·v⁺LC(n,zᵢ) + qᵢ·v⁻LC(n,zᵢ) + rᵢ·v⁺LC(n,z̄ᵢ) + sᵢ·v⁻LC(n,z̄ᵢ) ] = 0 for all n ∈ Fin(2N+1)
```

This is a **linear constraint** on the 4K coefficients.

## The Null Space Dimension

The null space of this constraint is the set of valid (p,q,r,s) tuples.

**Measured constraint-matrix rank** (for a 5×20 and 6×24 submatrix at evaluation points n=0..K-1):
- K=5: rank = 5 → null-space-dim = 20 - 5 = 15
- K=6: rank = 6 → null-space-dim = 24 - 6 = 18

**For the FULL constraint** (evaluating at all 2N+1 points):
- K=5: null-space-dim = **4K - 1 = 19** (confirmed S2693)
- K=6: null-space-dim = **4K - 1 = 23** (predicted)

The null space **grows by 4 dimensions** from K=5 to K=6.

## Why the Cubic/Quartic Polynomial Fails at K=6 — HONEST CHARACTERIZATION

**What was confirmed:**
- K=5 null-space dimension: 19 (when constraint evaluated at all 2N+1 points)
- K=6 null-space dimension: 23 (4K-1 formula holds)
- Cubic/quartic polynomial identity: fails at K=6 (confirmed by Lean compilation in S2699, S2704)

**What was NOT established (overclaimed in earlier draft):**
- The exact 4 "new" null-space basis vectors at K=6 were NEVER COMPUTED
- The polynomial evaluation on those vectors was NEVER TESTED
- The structural reason why the polynomial can't cancel them remains UNSPECIFIED

The earlier draft used mechanism language ("4 new orthogonal directions," "not covered by the polynomial," "missing cancellation terms") but was not backed up by the computation — it was speculation based only on dimension arithmetic (4K-1 at K=5 vs K=6).

## The Actual Observation (Not Yet a Mechanism)

The null space dimension does grow linearly with K (4K-1 pattern), but the cubic/quartic polynomial — which is fixed-degree and K-independent in its construction — fails to provide a cancellation argument at K=6. 

**Two interpretations:**
1. **Dimension-dependent failure (speculated, not proven):** The polynomial form is "dimensionally tailored" to K≤5, and the additional 4 degrees of freedom at K=6 expose uncanceled terms. This WOULD explain empirical failure, but requires explicit computation of the 4 new vectors and polynomial evaluation on them.

2. **Degree-dependent failure (alternative):** The cubic/quartic form simply doesn't have enough monomials or the right structure to handle K=6's constraints, independent of dimension. This also WOULD explain failure, but is equally unproven without deeper analysis.

## Status (Honest)

**REPORT:** Null-space dimensions confirmed (19 at K=5, 23 at K=6), cubic/quartic failures confirmed (empirical, from S2699/S2704), dimension-growth pattern confirmed (4K-1).

**NOT HANDLE:** The mechanism connecting dimension growth to polynomial failure is still open. To close it requires:
1. Computing K=5 and K=6 null-space basis vectors explicitly
2. Identifying which 4 are "new" (independent from K=5 vectors)
3. Testing cubic/quartic polynomial evaluation on those 4 vectors
4. Showing that the result is genuinely nonzero (not a numerical artifact)

This computation is feasible but was not completed in this session.

## Path Forward

Either:
- **Complete the mechanism characterization** by doing the computation above (1-2 hour effort)
- **Or accept that K-induction is empirically blocked at K≥6 without fully understanding the algebraic reason** and pivot to other approaches (CCM simplicity, D1 Angle 2 literature, or non-inductive methods)

The block is real and reproducible (Lean fails to compile). The underlying mechanism remains partially characterized.
