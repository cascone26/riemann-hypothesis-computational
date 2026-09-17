# Next Moves: Hodge-Signature Rigidity Investigation (Post S2660)

**Status:** Rigidity conjecture well-motivated by classical theory; 4 obstructions identified; toy computation verified. 
**Next layer:** Deep technical analysis of each obstruction.

---

## Move 1: Formal Berkovich Argument

**Goal:** Extend semicontinuity theorems from classical families to formal q→1 limits.

**Action:**
- Study Huber's adic spaces or Berkovich's non-Archimedean analytic geometry
- Key question: If X → S is a family of formal schemes (where S = Spec(Z[[q-1]]) or similar), do rank semicontinuity + Hodge-number constancy theorems apply?
- Look for: Riemann-Roch on formal families, monodromy on adic spaces, Lefschetz on non-Archimedean varieties

**Papers to consult:**
- Huber, "Étale cohomology of rigid analytic varieties and adic spaces" (1996)
- Berkovich, "Spectral theory and analytic geometry over non-Archimedean fields" (1990)
- Abramovich-Karu, "Weak semistable reduction in characteristic 0" (2000) [discusses specialization in families with bad reduction]

**Expected outcome:**
- **Best case:** Theorem stating that rank/signature rigidity holds in formal limits
- **Worst case:** Precise statement of where classical theorems break, opening a new direction

---

## Move 2: Archimedean Tail-Order Analysis in Kurokawa Limit

**Goal:** Understand whether archimedean place introduces new signature information as q→1.

**Action:**
- Read Kurokawa's absolute tensor product construction (multiply zeta functions at all places, including archimedean)
- Read Connes-Consani's F1 framework + how they handle the archimedean place
- Compute: What IS the signature of the archimedean contribution in the limit?
- Does it combine additively with finite-place signature, or are there mixing terms?

**Key insight from S2660:** The finite places give (1, ρ-1). The question is whether the archimedean place (which becomes central as q→1) has its own signature structure that could modify the total.

**Papers:**
- Kurokawa, "Multiple zeta functions and analytic continuation," *Proc. Japan Acad.*, 1992+
- Connes-Consani, "The tensor product in the category of F1-modules" + surrounding papers (2010-2015)
- Deninger, "Motivic L-functions and regularized determinants," *Motives and L-functions*, 2013
- Groskin 2607.02828 on archimedean tail orders (already in project)

**Expected outcome:**
- Either: The archimedean signature is well-defined and combines nicely (rigidity preserved)
- Or: Precise formula for how new archimedean structure modifies the total (guides next step)

---

## Move 3: Explicit Rank Computation for Curves as q→1

**Goal:** Verify that rank(NS(C_q × C_q)) actually STAYS CONSTANT (doesn't drop) as q→1.

**Action:**
- Take a concrete curve C over Z (e.g., genus 2 hyperelliptic, Cremona model)
- For a sequence of primes q = 2, 3, 5, 7, ..., 97, compute rank(NS(C_q × C_q))
- Plot rank vs q; does it stabilize or decrease?
- If rank varies: characterize the pattern (does it drop at special primes? does it depend on reduction type?)

**Concrete test case: y² = x⁵ + x + 1 (genus 2)**
- Good reduction away from 2,5 (possibly others)
- Picard group of Jacobian Jac(C) is well-studied
- NS(C × C) related to endomorphisms of Jac

**Code:**
- Use SageMath or Magma to compute Picard number (or bounds via L-function of curves)
- LMFDB queries: curve hypergeometric L-functions, Picard ranks

**Expected outcome:**
- **If rank is constant:** Verifies Voisin's prediction, supports rigidity
- **If rank varies:** Identifies the precise primes where jump occurs, constrains the form of archimedean contribution

---

## Move 4: Stability of τ's Near-Kernel in the Limit

**Goal:** Check whether the 23-dim near-kernel of τ (in the ζ case) remains stable as we approach a q→1 analog.

**Action:**
- The ζ operator τ = A+B−P has:
  - Spectral gap ≈ 1.87e-4
  - Kernel dimension 23 (near-kernel: everything with eigenvalue |λ| < gap)
  - This is specific to the Weil form, not generic to Hodge-index signatures
  
- Question: Is this kernel an artifact of the finite-prime structure, or does it persist in a q→1 limit?
  
- Approach: Build a family τ_q parameterized by q, compute kernel dimension for q = 2, 3, 5, ..., check trend as q→1

**Expected outcome:**
- **Most likely:** Kernel is ζ-specific, not a general Hodge phenomenon. Stays small (≤ 23) or drops as q→1.
- **If kernel grows:** Signals that signature becomes ill-defined in the limit (obstruction 4 confirmed)

---

## Move 5: Synthesis—Toward a Handle-Tier Proof

**Goal:** Once obstructions 1–4 are analyzed, synthesize into either:
- (A) A rigorous proof that signature is rigid in formal q→1 limit, OR
- (B) A precise characterization of the obstruction(s)

**Path if (A):** 
- Combine Berkovich semicontinuity + Kurokawa-limit archimedean analysis + rank stabilization
- Write a formal theorem: "For C/Z with good reduction, rank(NS(C_q × C_q)) and hence signature remain constant as q→1 in the formal Kurokawa framework"
- Corollary: F1-geometry SHOULD preserve (1,1) signature, giving a theoretical gateway to RH

**Path if (B):**
- Document exactly where/why rigidity breaks
- Propose a modified conjecture that accounts for the obstruction (e.g., "signature is rigid modulo archimedean contribution" or "rank stabilizes except at finitely many q")
- Opens a new research direction

---

## Time Estimates & Complexity

| Move | Complexity | Time | Blocker |
|------|-----------|------|---------|
| 1 (Berkovich) | High | 2–4 hours reading + 1–2 hours synthesis | Literature access |
| 2 (Archimedean) | Very high | 4–6 hours deep Kurokawa + 1–2 hours computation | Need Kurokawa/Connes-Consani original papers |
| 3 (Rank computation) | Medium | 2–3 hours SageMath coding + data pull | LMFDB access, SageMath setup |
| 4 (Kernel stability) | Medium | 2–3 hours implementing family τ_q | Existing τ code (available) |
| 5 (Synthesis) | Medium | 2–3 hours once moves 1–4 done | Completion of moves 1–4 |

---

## Honest Risk Assessment

**Success condition:** Any of the following resolves the obstruction question:
- ✓ Rigorous Berkovich proof that rank/signature are rigid
- ✓ Explicit archimedean formula showing signature is preserved
- ✓ Rank computation confirming stabilization across many primes
- ✓ Even a partial result (e.g., "rank stabilizes except for finitely many q") that sharpens the picture

**Failure modes (less likely, but possible):**
- ✗ Rank(NS) genuinely drops as q→1 in a way that breaks rigidity
- ✗ Archimedean place introduces a competing signature that modifies (1, ρ-1)
- ✗ Formal limit framework is too degenerate for classical theorems to apply

**If any failure mode occurs:** We'll have a precise obstruction, which is itself valuable progress. The door isn't "closed," it's "here's exactly what's blocking it."

---

## Recommendation

**Best bang-for-buck next:** Start with Move 3 (explicit rank computation) because:
1. It's concrete and immediate (SageMath code, LMFDB queries)
2. If rank stabilizes, it's strong evidence that rigidity holds (motivates deeper work on moves 1–2)
3. If rank drops, it identifies the precise primes where Picard rank changes, which constrains the archimedean contribution (guides move 2)
4. Result is dispositive: either clear support for rigidity, or a hard constraint

**Then:** If Move 3 is positive, pursue Move 2 (archimedean analysis) to understand the deeper mechanism. If Move 3 is negative, use the rank-drop pattern to guide a modified conjecture.

**Move 1 (Berkovich)** and **Move 4 (kernel)** are supporting research that strengthen whichever of 2 or 3 succeeds.

