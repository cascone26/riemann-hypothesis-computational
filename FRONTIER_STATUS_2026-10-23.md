# RH Proof Frontiers — Status & Blockers (2026-10-23)

| Lane | Primary Blocker | Nature | Status | Next Concrete Step | Est. Effort |
|------|-----------------|--------|--------|-------------------|-------------|
| **D1/heig** | Classical theorem: Weyl law ⟹ b_n>0 ∀n | Lit search | OPEN (bounded) | Search Simon/Szegö/Gautschi monographs | 2-4h |
| **D1/heig** | b_n>0 unconditional extension (n≥5) | RH-equivalent | RESEARCH | Treat as independent lemma; prove Hankel ⪰0 ∀k | High (≈RH) |
| **K-induction** | vPlusLC/vMinusLC linear independence at K≥6 | Algebraic phase boundary | BLOCKED | Compute K=6 null-space vector, characterize structure | 1-2h |
| **K-induction** | K-dependent cubic/quartic reformulation | Polynomial reconstruction | OPEN | Can the polynomial structure adapt to K-varying null space? | High |
| **vPlusLC** | K=6/7 Vandermonde determinant = 0 (exact) | Genuine singularity (dps=1000) | VERIFIED | Understand WHY the determinant collapses at K≥6 | Moderate |
| **vPlusLC** | General-K induction from the zero determinant | Structural characterization | OPEN | Prove independence ∀K<6, find obstruction ∀K≥6 | High |
| **CCM operator** | Simplicity (gap nonzero, isolated eigenvalue) | Analytical proof | BLOCKED | P-F fails (mixed signs), perturbation fails (gap ratio >1) | High (open) |
| **CCM operator** | Evenness (eigenvector γ-invariant) | Conditional | PROVEN | Depends on simplicity; numerically verified | — |
| **CCM operator** | Approximation (k_λ tracks ξ zeros) | Analytic bound | NUMERICALLY VERIFIED | Prove convergence ∀λ,N unconditionally; currently RH-equivalent | High (≈RH) |
| **Spectral (Angle 3)** | H(t) extraction from ξ | De Branges machinery | HARD | Would equal RH proof in another language; not a shortcut | Very High |

---

## Key Patterns

1. **RH-equivalence tier:** D1 (unconditional b_n>0 ∀n), CCM (all three conditions unconditional), Spectral (H⪰0).
   - These are not bugs or dead ends — they're CORRECT statements that proving them = proving RH.
   - Real value: different algebraic language, structural insight, potential new lemmas.

2. **Algebraic bottlenecks tier:** K-induction (K≥6 linear dependence), vPlusLC (same determinant zero).
   - Both hit genuine **phase transitions**, not technique limitations.
   - K-induction was top priority Sept 2026; now understood to be K-bounded, max K=5.

3. **Literature/computation gaps tier:** D1 (classical OPs theorem), CCM (simplicity proof path).
   - Circumscribed, bounded tasks with clear win conditions.

---

## Recommended Priority Order

### Immediate (high yield, bounded effort)
1. **D1 Angle 2 literature search** — If classical theorem found: D1 CLOSES. If not: reframe as open lemma.
2. **K=6 null-space characterization** — Reveals the algebraic structure underlying the phase boundary.

### High-value open questions (research-scale)
3. **K-induction reformulation** — Can polynomial structure adapt to K-dependent null space? (Architecture question, not proof.)
4. **CCM simplicity** — Only one missing condition; recent papers (Chatzakou-Helffer) are relevant but not decisive.

### Unlikely to yield (no shortcut)
5. Spectral/canonical-systems H(t) extraction — Valid, but equivalent to full RH proof.
6. Moment-based D1 extension (Hankel ⪰0 ∀k) — Exponential decay floor at k≥5 prevents numerics; analytically ≡ RH.

---

## Session Agents in Flight

- **a16be8eae93ea56ec:** K=6 null-space vector computation (Vandermonde decomposition, structure analysis)
- **abd8fbda61cf1e6a6:** D1 literature hunt (Weyl law ⟹ b_n>0 classical theorem search)

Expected completion: ~20-30 min from launch (2026-10-23, ~20:20 CDT).
