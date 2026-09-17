# Candidate (C): Infinite-Extent Canonical System Framework

**Objective:** Build a de Branges / canonical-system Hamiltonian H_∞(t) on [0,∞) that directly encodes ζ's explicit formula (arch+pole−prime) at full extent, testing whether the a→∞ limit yields κ=1 via pole-count.

**Status:** REPORT-tier theoretical framework + precise problem formulation. HANDLE-tier requires reading Suzuki papers (1204.1827, 2606.09096) and implementing numerical ODE solver.

---

## 1. The Gap Being Addressed

**Suzuki's finite-a construction:**
- H_a(t) defined on [0,a], with prime weights truncated: P(t) = 0 for t > a
- Result: κ=0 at every finite a (Thm 1.5: all W-zeros real)
- κ jumps to 1 only in the a→∞ limit (Cor 1.6)
- Suzuki's own papers: "this limit requires detailed arithmetic analysis beyond our approach"

**The open question:**
Does κ stay 0 as a→∞, or does it jump to 1, thereby encoding RH?

**Candidate (C) approach:**
Instead of truncating and taking a limit, build H_∞(t) directly on [0,∞) with no truncation. Then:
1. Check: is H_∞ limit-point at ∞ (Weyl m-function well-defined)?
2. If yes: compute κ(H_∞) via pole-count
3. Test: does κ(H_∞) equal 1 exactly for true primes?

---

## 2. The Canonical-System Framework

### 2.1 General form

A canonical system on [0,∞) is governed by:

```
dX/dt = J H(t) X,    X(0) = I

where J = [[0, -1], [1, 0]]  (symplectic structure)
      H(t) = 2×2 Hamiltonian with H(t) ≥ 0 (PSD)
      X(t) = fundamental solution
```

The **Weyl m-function** m(z) is extracted from the asymptotic behavior of X(t) as t→∞:

```
m(z) = lim_{t→∞} X₂₁(t) / X₁₁(t)  (or similar ratio, depends on normalization)
```

Key properties:
- m(z) is a **generalized-Nevanlinna function** if H is regular (real-symmetric, ℑm>0 in UHP)
- The **Nevanlinna kernel** K_m(z,w) = (m(z)−m̄(w))/(z−w̄) has negative eigenvalues
- **κ = #negative eigenvalues** of the Nevanlinna kernel (Krein–Langer)
- Krein–Langer theorem: κ = #{UHP poles of m(z)} (simple poles only)

### 2.2 For the explicit formula

For ζ, the natural Hamiltonian is:

```
H(t) = [[0, -ρ(t)], [σ(t), 0]]

where σ(t) = Ω(t)/2π - P(t)
      ρ(t) = weight function (to be determined from Suzuki papers)
```

- **Ω(t)/2π:** archimedean part (digamma/functional-equation factor)
- **P(t):** von Mangoldt weights (prime powers)
- **At finite a:** P(t) = 0 for t > a (truncation)
- **At infinite a:** P(t) = full sum for all t (Candidate C)

---

## 3. Limit-Point vs Limit-Circle Classification

### 3.1 The dichotomy

For a canonical system on [0,∞):

**LIMIT-POINT:** The Weyl m-function is **uniquely determined** by H(t) alone, with no additional boundary conditions at ∞. There is exactly one (up to scaling) square-integrable solution to the associated Schrödinger-type equation.

**LIMIT-CIRCLE:** The m-function is **not unique**; two linearly-independent square-integrable solutions exist, requiring extra boundary data at ∞.

### 3.2 Standard tests

1. **Trace-integrability:** ∫₀^∞ tr(H(t)) dt < ∞ ⟹ limit-point
   - For our H = [[0, -ρ], [σ, 0]], tr(H) = 0, so **this test is inconclusive**

2. **Growth of eigenvalues:** If eigenvalues of H(t) decay fast enough, limit-point
   - Our σ(t) oscillates but doesn't diverge; growth appears bounded

3. **Weyl m-function convergence (NUMERICAL TEST):**
   - Solve the system on [0, a] for increasing a
   - Extract m_a(z) at a test point z ∈ ℂ with ℑz > 0
   - **If m_a(z) → m_∞(z) as a→∞:** limit-point ✓
   - **If m_a(z) oscillates/diverges:** limit-circle or pathological

### 3.3 Computational test needed

```python
for a in [1.0, 2.0, 3.0, 4.0, 5.0, ...]:
    H_a(t) = Suzuki Hamiltonian on [0,a]
    m_a(z) = Weyl m-function of H_a at test point z
    δ_a = |m_a(z) - m_{a-1}(z)|
    
if δ_a → 0: LIMIT-POINT → m-function well-defined
if δ_a oscillates: LIMIT-CIRCLE or numerical issue
```

---

## 4. Pole-Count Route (if Limit-Point)

If H_∞ is limit-point, then:

1. **m_∞(z) is uniquely determined**
2. **It is a generalized-Nevanlinna function** (real-symmetric + positive imaginary part UHP)
3. **Krein–Langer theorem applies:** κ = #{poles of m_∞ in UHP}
4. **Generalized-Nevanlinna theorem:** a Nevanlinna function with a single simple pole in ℝ̂ has κ=1
5. **Argument principle:** count poles via winding number of log(m(z))

The pole of m_∞ should be at z = s (height coordinate corresponding to ζ's s=1 pole).

**RH consequence:** If κ(m_∞) = 1 iff all zeros of ζ are on the critical line, then:
> **RH ⟺ κ(H_∞) = 1**

This would be non-circular if κ is forced by arithmetic (true primes) and breaks under perturbation.

---

## 5. Faithfulness Test (The Critical Validation)

The pole-augmentation (Candidate B) showed κ=1 is GENERIC under the Generalized-Nevanlinna theorem — adding ANY nonzero pole to κ=0 gives κ=1, regardless of the pole location.

**For Candidate C to be non-circular, we need:**

```
Faithfulness: κ(H_∞) = 1 IFF primes are true (RH)

Test:
1. Build H_∞(t) with true primes P(t) → observe κ = ?
2. Perturb primes: P(t) → P(t) + δP(t), δ small but real
3. Recompute κ(H_∞^perturbed) → does it CHANGE?

If κ_perturbed ≠ 1:
  ⟹ H_∞ is FAITHFUL (captures RH-specific arithmetic)
  ⟹ non-circular construction

If κ_perturbed = 1 (persists under all perturbations):
  ⟹ κ=1 is STRUCTURAL, not arithmetic-contingent
  ⟹ generic Nevanlinna result (circular)
```

**Prediction from Candidate (B) findings:**
Candidate (B)'s pole-augmentation showed κ=1 persists even under 1% prime perturbation. If Candidate (C) shows the same robustness, then it's also generic (not non-circular).

---

## 6. Detailed Roadmap: What Must Happen Next

### 6.1 Prerequisites (REPORT-tier work)

1. **Read Suzuki 1204.1827, §2-3:**
   - Extract the explicit form of ρ(t)
   - Confirm normalization conventions
   - Understand boundary conditions at t=0 and t→∞

2. **Read Suzuki 2606.09096, §2-3:**
   - Screw function definition and properties
   - Connection between screw zeros and ζ zeros
   - The a→∞ limit in Suzuki's own words (why it's open)

3. **Formalize H_∞(t):**
   - Write down H_∞ precisely with all constants
   - Verify it matches Suzuki's construction at finite a
   - Check boundary conditions at 0 and ∞

### 6.2 Numerical implementation (REPORT→HANDLE)

1. **ODE solver:** 
   - Solve dX/dt = (z·J·H_a(t))·X on [0,a] for each a
   - Use high-precision arithmetic (Arb, mpmath)
   - Extract m_a(z) from solution couple

2. **Convergence test:**
   - Plot m_a(z) vs a for test points z = {15+i, 20+i, 30+i}
   - Compute δ_a = |m_a - m_{a-1}|
   - Fit: does δ_a decay exponentially, polynomial, or stall?

3. **Classification:**
   - **Exponential decay:** strong evidence of limit-point
   - **Polynomial decay:** potential limit-point with slow convergence
   - **Oscillation/growth:** limit-circle or numerical pathology

### 6.3 Pole-count computation (if limit-point)

1. **Nevanlinna kernel:**
   - Use corrected kernel code from weyl_nevanlinna_kappa.py
   - Test points should be near the real axis (y~0.05) to resolve near-line poles

2. **κ computation:**
   - Extract κ via eigenvalue count of Nevanlinna kernel
   - Verify via winding-number computation (S2637 / suzuki_winding.py pattern)

3. **Verification:**
   - Cross-check against known finite-a case: m_finite-a should have κ=0 (or expected value)
   - Test controls: inject known poles, verify κ changes accordingly

### 6.4 Faithfulness test (the hard part)

1. **Perturbation framework:**
   - δ-parameterized prime set: log p → log p · (1 + δ·ξ_p), ξ_p ~ N(0,1)
   - Rebuild H^δ(t) for each δ ∈ {0, 0.001, 0.005, 0.01}
   - Recompute m_δ and κ_δ

2. **Result interpretation:**
   - **κ_δ changes with δ:** faithful, non-circular ✓
   - **κ_δ = 1 for all δ:** generic, circular ✗

3. **Expected outcome (prediction from Candidate B):**
   κ_δ likely stays at 1 even under perturbation, making this generic too.

---

## 7. Known Blockers and Honest Limits

### 7.1 Theoretical blockers

1. **Suzuki's own open problem:** The papers explicitly state the a→∞ limit is open. We're not filling a gap they overlooked; we're tackling their acknowledged frontier.

2. **Generic Nevanlinna theorem:** Adding any pole to κ=0 gives κ≥1 automatically. This makes "κ=1 ⟺ RH" appear circular unless arithmetic forces the pole location uniquely.

3. **Faithfulness paradox:** Candidate (B) showed κ is robust to perturbation, suggesting the π₁ structure is inherent to the band-limited truncation, not RH-specific.

### 7.2 Computational limits

1. **Numerical ODE solving:** Complex-coefficient DEsystems can have stiff behavior. May require implicit solver + high precision.

2. **Pole resolution:** Weyl m-functions near the real axis are numerically sensitive. Requires test points with ℑz ~ 10⁻² and accuracies ~10⁻¹⁰+.

3. **Arb precision:** Full convergence test over a ∈ [1, 10] with 1400-bit precision could be expensive. Estimate: 10-100 cpu-hours.

---

## 8. Realistic Success Criteria

### 8.1 REPORT-tier outcomes (likely achievable)

- ✓ Extracted explicit H_∞(t) form from Suzuki papers
- ✓ Formalized limit-point criterion for this system
- ✓ Identified necessary tests (Weyl convergence, pole-count)
- ✓ Precise problem statement: "does m_a→m_∞ as a→∞?"

### 8.2 HANDLE-tier outcomes (realistic with effort)

- ? Implemented ODE solver for test a values
- ? Shown (or refuted) limit-point via convergence test
- ? If limit-point: computed κ(H_∞) ≥ 0 at true primes
- ? Confirmed or refuted faithfulness via perturbation test

### 8.3 RH-proof (not expected, impossible in current session)

- Would require: κ(H_∞) = 1 ONLY at true primes, AND proof that κ=1 ⟹ all zeros on critical line
- The latter is precisely the pole-count argument Suzuki identified as open
- This is the field's barrier, not unique to this approach

---

## 9. Comparison to Candidates (A) and (B)

| Candidate | Route | Blocker | Result |
|-----------|-------|---------|--------|
| **(A)** Inverse-spectral | Rebuild H from M_σ spectral data | Category error: Jacobi → real spectrum → no pole | WALLED: cannot encode pole via finite-rank reconstruction |
| **(B)** Pole augmentation | Add rank-2 pole to Suzuki's κ=0 m | Generic theorem: any pole works | REPORT: κ=1 automatic, not faithful to perturbation |
| **(C)** Infinite-extent Hamiltonian | Build H_∞(t) on [0,∞) directly | Open problem: does a→∞ limit exist? | **RESEARCH-SCALE:** needs numerical + faithfulness test |

---

## 10. Why (C) is Worth Pursuing

1. **Conceptually clean:** No truncation approximation, no ad-hoc pole terms. Just a canonical system on the natural domain.

2. **Addresses Suzuki's own gap:** Directly attacks the stated open problem in 1204.1827 / 2606.09096.

3. **Faithfulness testing is built-in:** If κ changes under perturbation, it's non-circular. If it doesn't, that's a precise negative result (refutes non-circularity).

4. **Classical machinery:** Canonical systems, Weyl m-functions, and Krein–Langer theory are well-established. Implementation difficulty is numerical, not theoretical.

---

## 11. Next Session Action Items

**Immediate (before numerical work):**
1. Read Suzuki 1204.1827 §2-3, extract ρ(t) and boundary conditions
2. Read Suzuki 2606.09096 §2-3, understand a→∞ limit statement
3. Write operators/suzuki_H_infinity_specification.py: explicit H_∞(t) with all constants

**Short-term (1-2 sessions):**
1. Implement ODE solver for dX/dt = (z·J·H_a)·X with Arb precision
2. Run convergence test: plot m_a(z) vs a for a ∈ [1,5]
3. Determine: limit-point or limit-circle?

**Medium-term (if limit-point):**
1. Compute κ(H_∞) via Nevanlinna kernel at test primes
2. Run faithfulness test: δ-perturbed primes, does κ_δ change?
3. Interpret results: non-circular or generic?

---

## 12. Conclusion

Candidate (C) is a well-defined research problem with:
- Clear theoretical framework (canonical systems, Weyl theory)
- Precise numerical tests (convergence, pole-count, faithfulness)
- Honest acknowledgment of blockers (Suzuki's own open problem, generic Nevanlinna theorem)
- Realistic effort budget for REPORT→HANDLE transition

**Not a proof of RH.** But a rigorous investigation of whether the a→∞ limit of Suzuki's construction yields the pole-count structure that would make an RH equivalence non-circular.

**Optimism line:** The limit-point question is well-posed and numerically testable. The faithfulness test has clear pass/fail criteria. The machinery exists. This is honest research, not speculation.

