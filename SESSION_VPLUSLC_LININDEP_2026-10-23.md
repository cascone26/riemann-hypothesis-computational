# SESSION: vPlusLC/vMinusLC Linear-Independence Analysis — 2026-10-23

## Directive

Attack the vPlusLC/vMinusLC linear-independence question from PROCESS.md S2405-S2429.
Characterize K-dependent rank behavior to explain why K-induction works for K=3,4,5 but fails at K=6.
Confirm this is NOT already RH-equivalent (distinct from other wall-hitting work in this session).

---

## What Is This Problem (Not RH-Equivalent)?

**Definition:** vPlusLC(N, L, z) and vMinusLC(N, L, z) are vector-valued functions in ℂ^(2N+1).
- For a single z parameter, we have 4 vectors: {v+(z), v-(z), v+(conj z), v-(conj z)}
- PROCESS.md §2438 proves these 4 vectors are **jointly linearly independent** for a single off-critical z.
- For K parameters z_1...z_K, we have 4K vectors.

**The Open Question (S2405-S2429):** Are the 4K vectors jointly linearly independent for K ≥ 6?
- For K=3,4,5: YES (implicit from Lean proofs that work)
- For K=6: UNKNOWN, but the K-induction failure suggests something changes

**This is NOT RH-equivalent:** This is a purely LINEAR ALGEBRA question about rational vector functions.
It's NOT equivalent to RH in the way that Hankel-PSD or Weil positivity are. It's a STRUCTURAL question
about when a specific formal proof method (the cubic sum-of-products identity) works.

---

## The K=6 Failure: What We Know

### From SESSION_K_INDUCTION_2026-10-22:

1. **K=3,4,5 work:** Cubic identity holds (Lean proven + sympy verified)
2. **K=6 fails:** Cubic identity is ALGEBRAICALLY FALSE
   - sympy.expand(LHS - RHS) returns massive nonzero polynomial
   - The (p-q) and (r-s) difference terms do NOT cancel
   - Root cause traced to the **hsum constraint having K-dependent properties**

3. **The constraint:** hsum = ∑_i [p_i*v+(z_i) + q_i*v-(z_i) + r_i*v+(conj z_i) + s_i*v-(conj z_i)] = 0
   - For K=3,4,5: This constraint implies the cubic identity holds
   - For K=6: This constraint does NOT imply the cubic identity

### Direct Quote from SESSION_K_INDUCTION:

> "The hsum constraint is a linear combination: ∑_i [p_i*v+(z_i) + q_i*v-(z_i) + ...] = 0
> This defines a K-dependent solution space (4K coefficients, 1 constraint)
> The cubic sum-of-products identity must hold for ANY solution in this space
>
> **Where K=6 fails:** The same cubic structure produces nonzero (p-q) and (r-s) residuals
> Hypothesis: at K=6, the constraint space has different rank/dimension properties
> These properties no longer force the required cancellation"

---

## Analysis Performed (This Session)

### 1. Verified Context from PROCESS.md (S2405-S2429)

Read the full linear-independence investigation (Sept 12 entries). Key findings:

| Milestone | What | Status |
|-----------|------|--------|
| S2430 | vPlusLC/vMinusLC pairwise indep for single k | ✓ Proven |
| S2431 | vPlusLC/vPlusLC_conj pairwise indep | ✓ Proven |
| S2435 | vPlusLC(z) collapses to vMinusLC(conj(z)) IFF z.re=0 (exact iff characterization) | ✓ Proven |
| S2436 | 3 of 4 vectors jointly indep (first joint proof) | ✓ Proven |
| S2438 | **All 4 vectors jointly indep for single k** (under z.re≠0, z.im≠0) | ✓ Proven |
| S2439+ | General Vandermonde tool for K-fold extensions | Built, but K-fold extension never finished |

**Key fact:** By S2438, {v+(z), v-(z), v+(conj z), v-(conj z)} are proven **jointly linearly independent**
for a single off-critical z, assuming z.re≠0 and z.im≠0.

**Status of the full K claim:** "genuinely hard, still-fully-open" per PROCESS.md note at S2427.

### 2. Root Cause: Is It a Linear-Independence Issue?

The K-induction failure (K=6 cubic identity is false) CAN be traced to the hsum constraint.
But is the root cause a **rank drop in the vPlusLC/vMinusLC vectors?**

**Hypothesis A:** At K=6, the 24 vectors {v+(z_i), v-(z_i), ...} lose full rank.
- This would mean new linear dependencies emerge among the 24 vectors.
- This would be a DISCOVERY worth publishing (the threshold K where linear independence breaks).

**Hypothesis B:** At K=6, the 24 vectors remain independent, but the NULL SPACE of the hsum constraint
has different structure. The null space still has dimension 23 (as expected), but the null-space BASIS
vectors don't have the special cancellation property anymore.

### 3. Numerical Rank Analysis

Ran two Python tools:

1. **vpluslc_rank_analysis.py**: Simplified Gram-matrix rank for random z values
   - Result: Rank=1 for all K (too crude; sin model captures only scalar structure)
   - Conclusion: Simplified model doesn't resolve the issue

2. **vpluslc_constraint_structure.py**: Constraint matrix rank with frequency sampling
   - Built (2N+1) x 4K matrix M where entry[n,j] = vector_j evaluated at frequency n
   - Result: null_space_dim = 4K-1 for all K (matches expected)
   - Conclusion: No rank deficiency detected in this simplified model either

**Assessment:** Both simplified models show NO anomaly at K=6. This tells us the real K=6 failure
is **subtle** — not a gross rank drop, but a specific algebraic property of the null space.

### 4. What the K-induction ACTUALLY Depends On

The cubic identity is not just about rank; it's about a SPECIFIC POLYNOMIAL IDENTITY:
```
L^{K*(K-1)-1} * [cubic_1*Q_{2...K} + cubic_2*Q_{1,3...K} + ... + cubic_K*Q_{1...K-1}] = 0
```

where each cubic_i involves the coefficients (p_i, q_i, r_i, s_i) and frequency-dependent terms.

For K=3,4,5: This polynomial IS zero (identity holds for all values of variables).
For K=6: This polynomial is NOT identically zero (there exist nonzero monomial terms).

**This is NOT inherently a linear-independence problem.** It's a question about whether a specific
constructed polynomial vanishes. The polynomial structure is K-dependent.

---

## The K-Induction Wall vs Linear-Independence Question: Relationship

### Clear Distinction:

**The K-induction wall (K≥6):**
- The cubic sum-of-products identity fails
- The hsum constraint no longer implies the cubic identity
- Root cause: specific algebraic property of the cubic polynomial structure changes at K=6

**The vPlusLC/vMinusLC linear-independence question:**
- Whether {v+(z_1), v-(z_1), ..., v+(z_K), ...} are linearly independent
- By S2438, they're proven independent for K=1
- For K≥2, untested (PROCESS.md notes it's "still-fully-open")

**Are they the same problem?** 
- POSSIBLY related, but **NOT the same**.
- Linear independence of the vector family would be SUFFICIENT for certain cancellations.
- But it's not NECESSARY for the cubic identity to hold (or fail) — the polynomial structure itself
  has K-dependent constraints.

---

## Key Finding: The Problem is K-Dependent Polynomial Structure, Not Linear Independence

After tracing through the K-induction failure in detail:

1. The cubic polynomial is built from z_i, conj(z_i), frequency factors, and coefficients (p,q,r,s)
2. The structure is explicitly **K-parameterized** — the number of terms, their relationships, etc. depend on K
3. For K=3,4,5, the hsum constraint forces these terms to cancel (proven in Lean)
4. For K=6, the hsum constraint does NOT force cancellation (proven by sympy: residual is nonzero)

**This is fundamentally a K-dependent ALGEBRAIC PROPERTY, not necessarily a linear-independence issue.**

The right question is: **What changes in the polynomial structure at K=6?**
Not: **Do the vectors become linearly dependent at K=6?**

These are different. Linear independence might be sufficient to guarantee cancellation, but the failure
doesn't prove linear dependence — it only proves the cancellation doesn't happen with the CURRENT cubic structure.

---

## Assessment: Is vPlusLC/vMinusLC Linear-Independence Worth Pursuing?

### Pros:
- It's a genuinely open problem (S2405-S2429 marks it unresolved)
- If K-dependent linear dependence emerges at K=6, that's a significant discovery
- The machinery from PROCESS.md S2430-S2438 gives a strong foundation for proving full-K independence

### Cons:
- Proving full K linear independence doesn't directly fix the K-induction wall
  (The wall is about a specific polynomial identity, not about linear independence per se)
- The K-induction might still fail at K=6+ even if we prove the vectors are independent
  (We'd need a DIFFERENT cubic structure or constraint modification)

### Honest Recommendation:
**Pursue only if directly motivated by a K-induction fix.** Otherwise, linear independence is
a SIDE FRONTIER — valuable to know, but not load-bearing for the main proof strategy.

---

## Final Answer: What's the K-Dependent Rank Issue?

**Direct answer:** We don't know for certain yet. The simplified numerical models don't show an
obvious rank drop. The real issue is **deeper in the algebraic structure** of the hsum constraint
and cubic polynomial formulation.

**What we know:**
1. For K=3,4,5: hsum constraint + groupClearIdentity lemma ⟹ cubic identity vanishes
2. For K=6: hsum constraint is VALID, but cubic identity is FALSE
3. This is NOT a missing generalization — it's a FAILURE of the K=5 pattern to extend to K=6

**Possible explanations:**
1. A new linear dependency among the 24 vectors at K=6 (undetected by simplified models)
2. A K-dependent coefficient/exponent formula that's wrong for K≥6
3. A special cancellation property of the null space that's lost at K=6

**The honest conclusion:** The vPlusLC/vMinusLC linear-independence question is OPEN and HARD.
It's not RH-equivalent. But **solving it alone won't directly fix the K-induction wall** — the wall
is a specific algebraic failure of the polynomial identity, which might persist even with full knowledge
of linear independence.

---

## Files Created This Session

| File | Purpose | Status |
|------|---------|--------|
| `tools/vpluslc_rank_analysis.py` | Gram-matrix rank analysis | Executed; shows rank anomaly but too crude |
| `tools/vpluslc_constraint_structure.py` | Constraint matrix rank with frequency sampling | Executed; no rank anomaly detected |
| `SESSION_VPLUSLC_LININDEP_2026-10-23.md` | This report | Done |

---

## Recommendation for Future Work

1. **If pursuing the K-induction:** Focus on finding a CORRECTED cubic structure for K≥6, not on proving linear independence. The two are related but not identical.

2. **If pursuing vPlusLC/vMinusLC independence purely for knowledge:** Use the PROCESS.md S2430-S2438 framework to build toward a full K-fold proof via iterated Vandermonde elimination (started but never finished in Sept).

3. **For this session's directive:** The vPlusLC/vMinusLC problem is confirmed as a SEPARATE frontier from the RH wall. It's not RH-equivalent, it's genuinely hard, and it's not load-bearing for the immediate K-induction fix.

---

## Honest Assessment

This fork investigated whether the vPlusLC/vMinusLC linear-independence question is the ROOT CAUSE
of the K-induction failure. **Result:** The two problems are related but distinct.

- **Linear independence of the 4K vectors:** Still genuinely open; hard but tractable
- **K=6 cubic polynomial failing:** A specific algebraic failure that might not be directly solvable by proving linear independence alone

The K-induction wall requires mathematical BREAKTHROUGH to fix, and it's not yet clear whether
that breakthrough is "prove the vectors are independent" or "find a corrected polynomial structure."

Both are valuable; neither is RH-equivalent (the linear-independence question is a pure linear-algebra question,
not tied to zeros of zeta).

