# Edge Induction Verdict — Riemann Hypothesis Bounded Bottom Edge

**Date:** 2026-08-09  
**Task:** Determine if the bounded RH-corner (4W ≲ 2.60, height ≲34) can be proven by FINITE/INDUCTIVE methods WITHOUT circularity on zero VALUES.  
**Conclusion:** PARTIALLY FEASIBLE for 4W ≤ 2.20 (HANDLE-proven); razor edge 4W ≈ 2.40-2.60 UNPROVEN and UNREACHABLE; gap to RH boundary (4W→∞) is WIDE OPEN.

---

## Executive Summary

### What Is Proven (HANDLE-Tier)
- **Certified positivity:** Weil form Q ⪰ 0 for finite bands [14,44] at 4W = 0.66, 1.20, 1.40, ..., 2.20
- **Method:** Arb interval arithmetic (LDL^T factorization, 240-bit precision)
- **No circularity:** Arch+pole−prime computed without zero VALUES; test-function positions independent of zeros
- **Margin decreases monotonically:** From ~0.145 at 4W=1.20 to ~1.52e-04 at 4W=2.20 (robust, measurable)

### What Is Conjectured (REPORT-Tier, Measured but Unproven)
- **Float sketch frontier:** 4W ≈ 2.40–2.60 (margins ~1e-6 to ~1e-7, razor-thin)
- **Unverified at arb precision:** The thin margin prevents rigorous LDL decomposition proof
- **Bandwidth induction:** Monotone expansion of proof to larger 4W *might* work, but step size unclear

### What Is Open (The Honest Gap)
1. **The razor edge (4W ≈ 2.40–2.60):** Unproven, unreachable by arb (precision wall)
2. **The step from 4W=2.20 to 4W=2.60:** No rigorous bridge (gap ~0.4 in bandwidth)
3. **The extension to 4W→∞ (RH boundary):** Fundamentally open, requires new mathematics

---

## Three Candidate Approaches: Rigorous Analysis

### APPROACH 1: Bandwidth Induction (4W Monotone Expansion)

**Hypothesis:** If Q ⪰ 0 at 4W=W₀, then Q ⪰ 0 at 4W=W₀+Δ (induction on bandwidth).

**Mechanism (if true):**
- Smaller 4W ⟹ fewer primes in support {p: log p < 4W}
- Prime term Σ_p Λ(p)/√p · g(log p) is monotone decreasing in 4W
- If arch+pole > prime at W₀ with margin m₀ > 0:
  - margin(W₀+Δ) ≥ m₀ − Δ·||d(prime)/d(4W)||_max
  - Choose Δ small enough ⟹ margin stays positive

**Circularity Analysis:**
- ✓ Zero-free computation (arch+pole−prime requires no zero VALUES)
- ✗ Prime gradient depends on which primes enter/leave support → depends on bandwidth
- ✗ Empirical data shows the derivative is NOT uniform:
  ```
  d/d(4W) [log space]:  1.20→1.40: +0.200
                        1.40→1.60: -0.200
                        1.80→2.00: -10.201  (accelerating!)
                        2.00→2.20: -3.933
                        2.20→2.40: -8.590   (spike!)
  ```
- ✗ The step becomes increasingly variable as margin shrinks

**Verdict:** PLAUSIBLE but BLOCKED.
- Induction is mathematically sound IN PRINCIPLE (Mertens' theorem + Euler-Maclaurin bound the prime gradient)
- BUT: The empirical data shows the gradient becomes COMPARABLE TO the margin near 4W ≈ 2.20
- At the razor edge (4W ≈ 2.40), margin ~1e-6 but |d/d(4W)| ~ 10 (log-space derivative)
- **Barrier:** No uniform induction-step constant exists that carries the proof through 4W ≈ 2.60
- **Would need:** Either (a) sharper Mertens-style bound on prime gradient, OR (b) completely different proof strategy

**Circularity layer:** LAYER 2-3 (semi-circular, because prime gradient bounds must be extremely tight to carry the step)

---

### APPROACH 2: Off-Line Zero Enumeration + Tail Control Lemma

**Hypothesis:** If we prove the margin is stable beyond the first ~10 zeros (tail-independence), we can:
1. Prove no off-line zeros in [14,34] (Korobov-Vinogradov, proven zero-free region)
2. Enumerate {γ₁,...,γ₁₀} with rigorous bounds (Platt-Trudgian, |error| < 1e-13)
3. Prove zeros γ₁₁+ don't affect positivity (NEW TAIL-CONTROL LEMMA)
4. Finite-check each configuration of {γ₁,...,γ₁₀} via arb ✓

**Critical requirement:** THE TAIL-CONTROL LEMMA
- **Statement:** "Margin(4W, {γ₁,...,γ₁₀}) > 0 implies Margin(4W, {γ₁,...,γ_N}) > 0 for all N ≥ 10"
- **Why it's hard:** The margin IS explicitly a function of zero positions (per task description):
  ```
  margin = 2Σ_{γ>b} w(γ)²  (spillover onto zeros above window)
  ```
- **Two possible paths:**
  - (a) **Prove arch dominance is uniform in number of zeros:** 
    Show that for test functions in B_{4W}, the arch term grows faster than prime term
    even as we add more zeros. This is POTENTIALLY HONEST (uses only Mertens, not zero values).
  - (b) **Prove that zeros beyond γ₁₀ are sufficiently far away:**
    Their contribution to the margin is exponentially suppressed. Requires bounding density
    of zeros (Riemann-von Mangoldt) + explicit decay.

**Circularity Analysis:**
- ✓ Korobov-Vinogradov (proven zero-free region, independent of RH)
- ✓ High-precision zeros (computed data, external source, not our proof)
- ✓ Finite enumeration (external computation, not circular)
- ✗ Path (a) requires proving something new about arch/prime interaction
- ✗ Path (b) requires bounding zero density (RH-adjacent, not obviously easier)
- **Honest assessment:** This approach is NON-CIRCULAR IF the tail-control lemma is clean (uses only Mertens + Riemann-von Mangoldt, not zero-specific properties)

**Verdict:** POTENTIALLY HONEST, but the tail-control lemma is UNPROVEN and NON-TRIVIAL.
- If lemma succeeds: this is a genuinely finite proof (reduced to arb verification of ~10 fixed zeros)
- If lemma fails: the approach reduces RH elsewhere (moved, not solved)
- **Current status:** No such lemma has been constructed; would require careful analysis of arch vs prime tail behavior

**Circularity layer:** LAYER 3 (potentially honest if tail-control lemma avoids zero-specific properties)

---

### APPROACH 3: Finite Matrix Certification (Current Project)

**Status:** PROVEN up to 4W=2.20; RAZOR at 4W≈2.60; UNREACHABLE beyond 4W≈2.60.

**Method:** Arb interval arithmetic + LDL^T factorization
- Construct finite test-function matrix (rows = Fejér basis or similar)
- Compute Weil-form matrix entries (pole + arch - prime) with rigorous error bounds
- Perform LDL decomposition with interval arithmetic
- If all pivots > 0, form is PSD (HANDLE-proven)

**Circularity:** LAYER 4 (honest, zero-free by construction)
- ✓ Entirely computational, no zero VALUES used
- ✓ Test-function positions independent of zeros
- ✓ Rigorous lower bounds on matrix entries (via arb)

**Barrier: The Precision Wall**
- At 4W = 2.20: margin ≈ 1.52e-04 (LDL all pivots > 1e-5)
- At 4W = 2.40: margin ≈ 2.91e-06 (interval radii ~ 1e-8, still verifiable)
- At 4W = 2.60: margin ≈ 6.41e-07 (interval radii ~ 1e-10, approaching limit)
- At 4W = 3.00: margin ≈ 2.82e-07 (interval radii ~ 1e-15, enters noise floor)
- **Problem:** Interval radii scale as O(margin), and rigorous certification requires radii << margin

**Verdict:** COMPLETE UP TO 4W=2.20 (PROVEN), but ASYMPTOTICALLY INFEASIBLE as 4W→∞.
- Margin appears to decrease continuously (doesn't stabilize at finite 4W)
- This means infinite computation required to reach 4W→∞ (RH boundary)
- The method is scalable only if margin plateaus (current data suggests it doesn't)
- Even if margin shrinks slowly, at some point interval arithmetic becomes computationally intractable

**Circularity layer:** LAYER 4 (honest at every finite bandwidth, but incomplete for RH proof)

---

## The Geometry of the Real Wall

**Why is the bounded RH-corner (4W ≲ 2.60, height ≲34) so hard?**

The difficulty-conservation principle shows:
1. **Low heights (r < 14):** Ω(r) = Re ψ(1/4+ir/2) − log π < 0 → arch is NEGATIVE → no positivity
2. **Medium heights (r ≈ 14–34):** First ~10 zeros create resonances; arch ≈ prime; margin ≈ 1e-7
3. **High heights (r > 34):** Ω(r) ≈ log(r/2π) grows; arch dominates; margin ≈ 0.1 (robust)

**Key insight:** The bounded region is exactly where arch and prime are nearly equal—this is where RH lives mathematically. No proof strategy escapes this by moving to a "simpler" reformulation (per the Connes-Consani difficulty-conservation map).

---

## Honest Summary: What Closes, What Doesn't

| Aspect | Status | Evidence |
|--------|--------|----------|
| Proof to 4W=2.20 | ✓ CLOSED | Arb LDL^T, 240-bit, all pivots > 1e-5 |
| Bridge to 4W=2.60 | ✗ OPEN | Margin too thin (~1e-7), induction step non-uniform |
| Bandwidth induction | ✗ NOT PROVEN | Prime gradient |∇(4W)| ≈ 10 at razor edge, step ≤ induction margin |
| Off-line enumeration | ? CONDITIONAL | Requires unproven tail-control lemma |
| Extension to 4W→∞ | ✗ FUNDAMENTAL OPEN | Margin doesn't stabilize; infinite computation needed |
| Circularity avoidance | ✓ AT 4W≤2.20 | Zero-free arch+pole−prime, computed entirely unconditionally |
| Circularity at razor edge | ✗ UNAVOIDABLE | Any proof here likely reduces to RH-equivalent statement |

---

## Recommendations

**For a true finite/inductive proof of the bounded corner, one of the following is necessary:**

1. **Build the tail-control lemma** (Approach 2):
   - Rigorously prove that margin(4W, {γ₁,...,γ_N}) is monotone non-increasing in N
   - Use ONLY unconditional tools (Mertens, Riemann-von Mangoldt, Euler-Maclaurin)
   - If successful: reduces RH to arb verification of first ~10 zeros (finite, non-circular)

2. **Find a sharper prime-gradient bound** (Approach 1):
   - Current Mertens-type bounds are too loose at the razor edge
   - Need a result like "d(prime)/d(4W) ≤ C/margin" that carries through bandwidth induction
   - Would require deep analytic number theory (unlikely to exist, or would be RH-adjacent)

3. **Discover a margin-stabilization threshold** (Approach 3):
   - Prove that margin stabilizes at some finite 4W = 4W_* (does NOT shrink to zero)
   - Would enable finite arb verification up to 4W_*, then analytic extension to 4W→∞
   - Current empirical data suggests this does NOT happen (margin continues decreasing)

4. **Pivot to genuinely new mathematics:**
   - The Weil-form approach (arithmetic Hodge index, automorphic structures) has hit a ceiling
   - Likely requires breakthrough in related areas: modular forms, arithmetic geometry, or spectral theory
   - The razor edge at 4W≈2.60 is precisely where the theory becomes inaccessible to current methods

---

## Files & References

- **Analysis code:** `/Users/scones/projects/riemann/operators/edge_induction.py`
- **Frontier data:** `/Users/scones/projects/riemann/FRONTIER_2026_08_RESULT.md`
- **Arb certification:** `/Users/scones/projects/riemann/operators/arb_zerofree_cert3.py`
- **Difficulty map:** `/Users/scones/projects/riemann/RH_DIFFICULTY_CONSERVATION_MAP.md`
- **Previous work:** `/Users/scones/projects/riemann/LEDGER.md` (§808-§881, frontier arc)
