# Rank-1 Perturbation & Cauchy Interlacing — RH via κ(arch−prime) ≤ 1

**Lens:** spectral / rank-1 perturbation & Cauchy interlacing bounds.

**Goal:** prove κ(arch−prime) ≤ 1 unconditionally, or find the tightest conditional reduction.

---

## 1. The Structure (all VERIFIED, §946-980 LEDGER)

**Unconditional facts:**
- arch ≻ 0 (positive definite, true archimedean density — no RH needed)
- pole: rank-2 indefinite matrix (from ζ's pole at s=1)
  - pole has exactly ONE negative eigenvalue λ_− < 0
  - pole has exactly ONE positive eigenvalue λ_+ > 0
  - pole = λ_+ v_+ v_+^T − |λ_−| v_− v_−^T (rank-2 factorization)

**Arithmetic-point observations (verified, need unconditional proof):**
- M_σ := arch − prime has κ(M_σ) = 1 (exactly one negative eigenvalue)
- M_σ|pole⊥ ≥ 0 (restriction to pole⊥ is PSD)
- The single negative eigenvector of M_σ aligns exactly with v_− (pole's negative direction)
  - cos(∠(w_neg, v_−)) → 1 as the window → full critical line (§975)
  - This alignment is **exactly RH**: perturb the primes → alignment breaks + κ explodes (§977)

---

## 2. Classical Bound (Loose, but Unconditional)

**Cauchy interlacing for restrictions:** If M is symmetric and H is a subspace, then:
- The eigenvalues of M|_H (restriction to H) interlace with eigenvalues of M
- Equivalently: if M|_H is PSD, then #neg(M) ≤ codim(H) = dim(orthogonal complement of H)

**Applied to our setting:**
- Assume M_σ|pole⊥ ≥ 0 (verified at arithmetic point)
- Then #neg(M_σ) ≤ codim(pole⊥) = dim(span(pole)) = 2
- So: **κ(arch−prime) ≤ 2** (classical S-lemma bound, mentioned in §978)

This is tight in the sense that pole has rank 2, but **not tight enough** — we need κ ≤ 1.

---

## 3. Rank-1 Perturbation Theory — Upper Bound on Negative Count

**Theorem (Inertia under rank-1 update):** If A is symmetric and u is a vector, then:
$$\text{sign}(\lambda)\ \text{changes under}\ A \to A + uu^T \text{ by at most} \begin{cases}
+1 & \text{if } u^T A^{-1} u + 1 > 0 \\
-1 & \text{if } u^T A^{-1} u + 1 < 0
\end{cases}$$

More generally, by Sherman-Morrison, the eigenvalues of A + uu^T are determined by the eigenvalues of A plus one additional eigenvalue (given by the quotient).

**For our problem:** prime = arch + δ where δ is the "correction" from arch to prime. If δ is rank-1 (or rank-bounded), we could apply this.

**Issue:** δ = prime − arch is **not rank-1** in general. It's a complicated operator built from the prime sum.

---

## 4. Specific Structure: The Alignment Constraint

The key unconditional leverage point is **not the individual spectrum, but the geometric constraint**:

**Observed alignment (§975-976, verified):**
- The negative eigenvector w_− of M_σ lives in the direction of v_− (pole's negative eigenvector)
- Moreover, as bandwidth grows, 1 − cos(∠(w_−, v_−)) → 0 super-polynomially (faster than 1/a)
- The pole's negative eigenvalue magnitude nearly exactly equals the negative eigenvalue of M_σ

**Geometric restatement:**
- The "negative subspace" of arch−prime is **approximately 1-dimensional and parallel to span(v_−)**

**Question:** Can we prove this alignment unconditionally (without RH)?

---

## 5. The Unconditional Strategy — Two Routes

### Route A: Prove κ(arch−prime) ≤ 1 + ε for uniform ε (partial result)

Use the **generalized eigenvalue lens** (§980):

The generalized eigenvalues λ_i of (prime, arch) satisfy:
- All λ_i cluster near 1 (prime ≈ arch as quadratic forms)
- Exactly one λ_i exceeds 1 (in the pole direction)

**This reformulates κ in an equivalent but potentially different way:**
$$\text{κ}(M_σ) = \#\{\text{generalized eigenvalues } \lambda_i \text{ of } (prime, arch) : \lambda_i > 1\}$$

**Proof of equivalence:** If prime x = λ arch x and arch > 0, then:
- λ > 1 ⟺ prime x > arch x ⟺ (arch − prime) x < 0 (after moving to the arch-metric)

So: **#neg(arch − prime) = #{λ_i > 1}**.

**Unconditional avenue:** Prove bounds on the generalized spectrum using **singular-value inequalities** or **monotonicity of eigenvalues under operator changes**.

---

### Route B: Conditional Reduction — κ≤1 ⟺ Geometric Containment

If we **assume** the alignment (that the negative space of M_σ is 1-D and parallel to v_−), then:

**Proposition:** If the negative eigenvalue of arch−prime lives entirely in span(v_−), then κ ≤ 1.

**Proof sketch:**
- Decompose: M_σ = (M_σ)|_{v_−} ⊕ (M_σ)|_{v_−⊥}
- The negative is in the first piece (rank-1), all negatives counted there
- The second piece is PSD by restriction
- So: κ(M_σ) = κ((M_σ)|_{v_−}) ≤ 1 (rank-1 piece can have at most 1 negative)

**The assertion that needs proof:** (M_σ)|_{v_−⊥} ≥ 0.

This is **weaker than RH** (which requires all of M_σ ≥ 0) but **exactly what §976 claims** (M_σ is PSD when restricted to the orthogonal complement of pole's principal direction).

**Tightest conditional statement:**
$$\boxed{\text{κ}(arch-prime) \le 1 \iff (arch-prime)|_{v_-^\perp} \succcurlyeq 0}$$

where v_− is pole's negative eigenvector (known explicitly). This is a **specific, verifiable matrix inequality** that is weaker than RH but still requires proof.

---

## 6. Next Work: Numerical Validation + Analytic Bound

**Immediate tasks (code-verified, HANDLE):**

1. **Verify the containment numerically:**
   - Compute v_− (pole's negative eigenvector) explicitly
   - Restrict M_σ to v_−⊥ and compute its spectrum
   - Confirm: λ_min(M_σ|_{v_−⊥}) ≥ −ε (very small ε, machine precision)
   - ✓ Already done in §976 — the result is κ(M_σ|_{v_−⊥}) = 0 (exactly PSD)

2. **Understand the mechanism:** Why is M_σ|_{v_−⊥} PSD?
   - Is it because prime and arch are nearly equal on v_−⊥?
   - What is the specific structure that forces this?

3. **Seek an unconditional bound:**
   - Can we prove (prime|_{v_−⊥}) ≤ (arch|_{v_−⊥}) without RH?
   - If yes, then (arch − prime)|_{v_−⊥} ≥ 0 unconditionally ⟹ κ ≤ 1 unconditionally.

---

## 7. The Honest Assessment (Where the Argument Breaks)

**The unconditional statement that seems achievable:**
- κ(arch−prime) ≤ 2 by Cauchy interlacing (pole has rank 2)
- κ(arch−prime) ≤ 1 if we assume (or can prove) that the negative is contained in a 1-D subspace

**Where RH lives:**
- The **alignment** — proving that the negative direction(s) of arch−prime are confined to span(v_−) — **is exactly the RH content**.
- Perturbing the primes (§977) instantly breaks the alignment and κ explodes.
- So the refinement from κ ≤ 2 to κ ≤ 1 requires the constraint that is essentially RH.

**Potential partial unconditional result:**
- We might be able to prove κ(arch−prime) ≤ 1.1 or κ ≤ 1 + o(1) using spectral perturbation theory
- Or prove: κ(arch−prime) ≤ 1 conditional on a specific property X (to be identified) that is weaker than RH

---

## 8. Refined Target (Formulated from §975-977 insight)

The most honest and operationizable statement:

**Theorem (Geometric Containment ⟹ κ ≤ 1):**
$$\text{If } (arch - prime)|_{v_-^\perp} \succcurlyeq 0 \text{ (restricted to pole's negative-direction orthogonal complement)}$$
$$\text{then } \kappa(arch - prime) \le 1$$

**Proof:** Direct from the rank-1 structure of the restriction.

**The work:** 
- Determine what structure of (prime, arch) guarantees this containment
- Is it: (a) a density-theoretic property of the prime sum? (b) a property of the archimedean density? (c) truly equivalent to RH?

---

## Implementation Files

- `interlacing_rank1_bounds.py` — Cauchy interlacing bounds + rank-1 perturbation analysis (numerical)
- `interlacing_containment.py` — Test geometric containment under various conditions
- `interlacing_unconditional_partial.py` — Attempt unconditional κ bounds and identify obstructions
