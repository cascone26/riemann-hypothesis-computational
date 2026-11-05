# Next Steps Analysis — All Lanes (2026-10-23)

**Purpose:** Explicitly scope the next concrete move for each frontier lane, ranked by tractability and impact. Updated after K-induction/vPlusLC unification finding.

---

## Lane Rankings (by tractability ⟹ impact ratio)

### 1. D1/heig — Literature Hunt (IMMEDIATE, 2-4h)
**Current state:** b_1..b_4 proven unconditionally via Hamburger+Favard+interval arithmetic. Extension to all n blocked by one missing classical theorem.

**Next step:** Search Simon (OPUC vol 4), Szegö, Gautschi, Denisov-Rakhmanov, Stahl-Totik for theorem: "smooth measure density-of-states satisfying Weyl law ⟹ all Jacobi coefficients b_n > 0."

**Win condition:** Find the theorem statement + verify applicability to ξ function.
**Failure condition:** After 4h deep search, conclude it's an open problem (= no shortcut exists, D1 requires RH-scale work).
**Impact if found:** D1 CLOSES entirely (unconditional route to RH).
**Impact if not found:** Honest re-framing of D1 as "unproven lemma needed," feeds into research-scale work.

**Agent in flight:** abd8fbda61cf1e6a6 (expected return ~20:35 CDT).

---

### 2. K-induction — Null-Space Characterization (IN PROGRESS, 1-2h)
**Current state:** K=6 Vandermonde determinant = 0 (verified dps=1000); hsum constraint uses those vectors; linear dependence breaks cubic/quartic.

**Next step:** Compute the explicit null-space vector(s) of the K=6 vPlusLC/vMinusLC matrix. Characterize structure:
- Which z_i blocks contribute most?
- Are there sign patterns, conjugate symmetries, repeating structure?
- Does the null vector reveal hidden algebraic structure?

**Win condition:** Understand the structural reason WHY the determinant collapses at K=6 (not just the fact that it does).
**Failure condition:** Null vector is "generic" with no obvious structure (suggests K=6 collapse is accident, not design).
**Impact if structure found:** May suggest K-dependent reformulation of cubic (accommodate the dependence instead of fighting it).
**Impact if generic:** Confirms K≥6 is an algebraic dead end for this vector family; pivot to non-polynomial approaches.

**Agent in flight:** a16be8eae93ea56ec (expected return ~20:35 CDT).

---

### 3. CCM Simplicity — Proof Path Survey (MODERATE, 4-6h)
**Current state:** Numerical evidence strong (gaps 6e-12 to 5e-7, always nonzero). Analytic proofs all fail (P-F: mixed signs, perturbation: gap ratio > 1).

**Next step:** Systematically survey what's REQUIRED for simplicity at the functional-analysis level:
- Can Toeplitz structure + functional equation constraints + symmetry properties prove gap > 0?
- Is there a recent (2020s) paper on generic/Baire-category simplicity for operators related to tau?
- Can we extract a lower bound on the gap from Gershgorin circle / spectral bounds applied to tau's structure?

**Win condition:** Find a proof path (even if challenging) that applies to our tau.
**Failure condition:** Confirm that simplicity for this operator is genuinely research-scale (no known shortcut).
**Impact if found:** CCM Simplicity solved ⟹ all three CCM conditions proven (evenness+approximation depend on it).
**Impact if not found:** Honest characterization of CCM as "spectral mechanism verified numerically, analytic closure blocked on one non-obvious gap."

**Next agent launch (when null-space returns):** Assign one agent to simplicity survey.

---

### 4. vPlusLC Independence — Characterize K=6-K=7 Boundary (MODERATE, 3-4h)
**Current state:** K≤5 proven independent (Vandermonde det nonzero). K=6 selective collapse (2/3 configs zero). K=7 universal collapse (all configs zero).

**Next step:** Once K=6 null vector is computed:
1. Compare the null vector structure to the hsum constraint in CrossK6Key.lean.
2. Determine whether the null vector's pattern appears in the polynomial identity's monomials.
3. Test whether similar null structure exists at K=7 (universal collapse suggests common cause).

**Win condition:** Identify a SINGLE structural principle explaining both K=6 (selective) and K=7 (universal) collapse.
**Failure condition:** K=6 and K=7 have unrelated degeneracy mechanisms (suggests independent phenomena).
**Impact if unified:** Strengthens hypothesis of K-bounded algebraic limit, may guide future non-polynomial approaches.
**Impact if separate:** Clarifies that multiple independent K-dependent barriers exist (not one phase transition).

**Dependency:** Awaits K-induction null-space result (agent a16be8eae93ea56ec).

---

### 5. K-induction Reformulation — Polynomial Redesign (HIGH EFFORT, open-ended)
**Current state:** Cubic/quartic proven to fail at K=6 due to hsum null-space change. Both tested on Olympus, both fail identically (Lean `ring` tactic).

**Next step:** Given the K=6 null structure (from agent task #2), explore:
1. Can the polynomial identity be rewritten in a form that's agnostic to the null-space dimension?
2. Are there K-dependent variants of the cubic/quartic that "sidestep" the K=6 wall (adapt structure to the changing dependence)?
3. Is there a genuinely different polynomial family (not degree-3/4) that avoids the K=6 collapse?

**Win condition:** A new polynomial form that compiles clean at K=6+ on Olympus.
**Failure condition:** Verify that the cubic/quartic collapse is FUNDAMENTAL to this approach (not just technical).
**Impact if won:** K-induction lane revives; can potentially prove RH-relevant identities for K≥6.
**Impact if lost:** K-induction formally K-bounded at K≤5; close the lane and accept this as a boundary.

**Estimated effort:** This is research-scale and open-ended; 8-16h minimum, possibly weeks.

---

## Recommended Flow (Priority Order)

1. **Wait for D1 literature search (abd8fbda61cf1e6a6) to return.** If theorem found: D1 CLOSED, major win. If not found: reframe as independent lemma.
2. **Wait for K-induction null structure (a16be8eae93ea56ec) to return.** Characterize the K=6 degeneracy pattern.
3. **Launch vPlusLC boundary characterization** (low-effort cross-lane synthesis): compare null vectors to hsum, test K=7.
4. **Launch CCM simplicity survey** (moderate-effort research): map out proof paths for gap > 0.
5. **Decide on K-induction reformulation:** based on null structure, determine whether redesign is tractable or the lane is K-bounded.

---

## Metrics for Success This Session

- **D1:** Theorem found (huge win) OR formally named as open problem (honest close).
- **K-induction:** Null structure computed AND compared to hsum constraint (cross-lane synthesis).
- **vPlusLC:** K=6/7 boundary mechanism unified (single principle OR separate phenomena confirmed).
- **CCM:** Simplicity proof path mapped (even if not closed).

**Overall:** Every lane has a clear next concrete step; no "stuck" or "where next" language.

---

## Files to Update When Agents Return

- LEDGER.md: Append both agent findings with S-numbers (likely S2715, S2716)
- STATUS.md: Add UPDATE 43 summarizing results
- SESSION_K6_UNIFICATION_2026-10-23.md: Add "K=6 Null-Space Structure" section with agent findings
- FRONTIER_STATUS_2026-10-23.md: Update lane status based on new information

**Anti-loop check:** None of the tasks above are redundant with prior sessions (all are downstream analysis of NEW findings from this session).
