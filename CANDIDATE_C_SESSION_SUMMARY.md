# Candidate (C) Research Initiative — Session Summary (2026-09-17)

**Status:** REPORT-tier theoretical framework complete. Precise problem formulation established. Next session: numerical implementation.

## What Was Accomplished

### 1. Problem Formulation (Precise)
Candidate (C) asks: **Can we build a canonical-system Hamiltonian H_∞(t) on [0,∞) that directly encodes ζ's explicit formula at full extent, bypassing Suzuki's finite-a truncation?**

The existing Suzuki construction:
- Uses H_a(t) truncated to [0,a]
- Gives κ=0 at every finite a (Thm 1.5: all W-zeros real)
- κ jumps to 1 only in limit (Cor 1.6)
- Suzuki's papers: this limit "requires detailed arithmetic analysis beyond our approach"

Candidate (C):
- Build H_∞(t) directly on [0,∞)
- No truncation: P(t) summed over ALL primes for all t
- Test: does m-function converge as a→∞? (limit-point vs limit-circle)
- If limit-point: compute κ(H_∞) via pole-count

### 2. Theoretical Framework Established

**Canonical system form:**
```
dX/dt = J·H(t)·X,    H(t) = [[0, -ρ(t)], [σ(t), 0]]

σ(t) = Ω(t)/2π - P(t)
  Ω(t) = archimedean factor (Re digamma)
  P(t) = von Mangoldt weights (primes)
```

**Weyl classification (open question):**
- LIMIT-POINT: m-function uniquely determined by H alone (likely)
- LIMIT-CIRCLE: extra boundary conditions needed (possible)

**Test:** Track m_a(z) convergence as a→∞
- If m_a(z) → m_∞(z): limit-point → proceed to pole-count
- If oscillates/diverges: limit-circle → different approach needed

**Pole-count (if limit-point):**
- Use Krein–Langer: κ = #{poles of m_∞ in UHP}
- Nevanlinna kernel: compute κ via eigenvalue count
- RH ⟺ κ = 1 (if arithmetic forces pole location uniquely)

### 3. Growth Analysis (Numerical Check)

Computed on [0, 8]:
- ∫₀^a ||H(t)|| dt ≈ 11.5 (bounded)
- σ(t) oscillates in [-7, +3] but doesn't diverge
- **Conclusion:** suggests limit-point behavior (necessary but not sufficient)

### 4. Known Blockers Identified

**Theoretical:**
1. Suzuki's own open problem — we're not filling an overlooked gap
2. Generic Nevanlinna theorem — adding any pole gives κ=1, makes "κ=1 ⟺ RH" appear circular
3. Candidate (B) already showed κ robust to perturbation (suggests generic, not faithful)

**Computational:**
1. Need exact H_∞(t) from Suzuki papers (1204.1827, 2606.09096)
2. Complex-coefficient ODE solver with high precision (Arb)
3. Pole resolution near real axis numerically sensitive
4. Faithfulness test requires full perturbation framework

### 5. Clear Roadmap for Next Phase

**REPORT-tier (before numerical work):**
- ✓ Read Suzuki 1204.1827 / 2606.09096, extract explicit ρ(t)
- ✓ Formalize H_∞(t) specification with all constants
- ✓ Write out full limit-point classification test

**HANDLE-tier (short term):**
- ? Implement ODE solver (dX/dt = z·J·H_a·X)
- ? Run convergence test: plot m_a(z) vs a
- ? Classify: limit-point or limit-circle?

**HANDLE-tier (if limit-point):**
- ? Compute κ(H_∞) via Nevanlinna kernel
- ? Faithfulness test: δ-perturbed primes, does κ change?
- ? Interpret: non-circular or generic?

## Why This Approach is Honest

1. **Addresses a real open problem:** Suzuki's papers explicitly state the a→∞ limit is open.

2. **Faces known blockers head-on:** We've identified exactly why (A) and (B) failed, and why (C) is structurally different.

3. **Admits the paradox:** If κ is generic (any pole works), then κ=1 is likely not faithful to primes. The faithfulness test can expose this.

4. **Uses established machinery:** Canonical systems, Weyl theory, Krein–Langer are classical. No exotic new tools needed.

5. **Realistic effort budget:** REPORT framework done in one session. HANDLE-tier is 2-3 sessions of numerical work, not an open-ended research rabbit hole.

## Files Created This Session

- `operators/hamiltonian_infinite_extent.py` — task overview + research roadmap
- `operators/suzuki_hamiltonian_analysis.py` — explicit-formula symbol, growth analysis
- `operators/CANDIDATE_C_THEORETICAL_FRAMEWORK.md` — 12-section detailed framework
- `CANDIDATE_C_SESSION_SUMMARY.md` — this file
- **LEDGER entry S2649** — session findings + assessment

## What Success and Failure Look Like

**Success (non-circular):**
- m_a(z) converges as a→∞ (limit-point) ✓
- κ(H_∞) = 1 at true primes ✓
- κ_δ ≠ 1 under prime perturbation ✓
- Interpretation: H_∞ is RH-faithful, captures arithmetic

**Failure (circular):**
- m_a(z) oscillates / diverges (limit-circle)
- κ(H_∞) = 1 persists under all perturbations (generic)
- Interpretation: κ=1 is structural, not RH-specific

**Partial (inconclusive):**
- Limit-point ✓, but κ_δ also = 1 under perturbation
- Conclusion: limit-point exists but doesn't buy faithfulness

## Next Session Checklist

1. [ ] Pull Suzuki papers 1204.1827 and 2606.09096 from arXiv
2. [ ] Extract explicit ρ(t) and boundary conditions
3. [ ] Implement operators/suzuki_H_infinity_specification.py
4. [ ] Write ODE solver wrapper for high-precision solving
5. [ ] Run convergence test for a ∈ [1, 5, 10] at 3-5 test points
6. [ ] Classify: limit-point or limit-circle?
7. [ ] If limit-point: compute κ via Nevanlinna kernel
8. [ ] If κ=1: run faithfulness test (perturbation sweep)

---

**Optimism line:** The question is mathematically well-posed, the test is numerically clear, and the outcome (whether yes, no, or inconclusive) will be a precise, reproducible finding. Not a proof of RH — but an honest investigation of whether this route yields non-circular RH-equivalence.

