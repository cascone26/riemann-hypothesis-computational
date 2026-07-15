# Track 2 (Connes finite spectral model) — rigorous status, 2026-05-19

## What is now ESTABLISHED (numerically, non-circularly)

Define, at λ²=14, L = log 14, basis V_n(x)=L^{-1/2}e^{2πinx/L} for |n|≤N, the
finite Weil matrix

    τ_N[n,m] = (1/2π) ∫_ℝ Ω(r) ê_n(r) conj(ê_m(r)) dr      (archimedean)
             + W_{0,2}(n,m)                                  (boundary / ζ pole)
             − Σ_{k≤λ²} Λ(k) k^{-1/2} q_{n,m}(log k)         (primes)

with Ω(r) = Re ψ(¼ + i r/2) − log π   (digamma symbol; **zero-free**),
ê_n the sinc transform of V_n.  **No zeta zero enters the definition.**

Verified (`connes_decisive_v3.py`, `connes_v3_noncirc_audit.py`,
`connes_v3_control.py`):

1. **Positivity.** λ_min(τ_N) collapses from v1's −3.6 to machine zero
   (≈ −1e-13; kernel eigenvalue ≈ −1e-14) for N = 16…64.  τ_N is
   positive-semidefinite with a single closing near-kernel direction ξ.

2. **Non-circular reproduction.** The secular function
   ξ̂(z) = (2/√L) sin(zL/2) Σ ξ_j/(z−2πj/L) has real zeros equal to the
   Riemann zeros γ₁…γ₁₈ to ≤ 1e-11 at N=48 — including γ₁₃…γ₁₈
   (59.347, 60.832, 65.113, 67.080, 69.546, 72.067) which appear in **no
   source file** and were checked against an independent `mpmath.zetazero`
   computation.  Genuine prediction, not a planted constant.  Resolution
   ceiling tracks the basis bandwidth s_N=2πN/L (expected).

3. **Critical-line specificity.** Replace Ω's argument ¼ by a≠¼ (i.e.
   pretend zeros lie on Re s = 2a).  Then λ_min jumps to −0.01…−0.60
   (9–11 orders of magnitude negative) and zero accuracy degrades from
   1e-11 to ~1e-3.  Only a=¼ — the true critical line — yields a PSD τ_N
   and exact zeros.  The positivity is therefore an informative RH signal,
   not a normalization artifact.

## The exact reduction (this is a theorem, Weil/Bombieri/Connes)

τ_N is the Gram matrix of the Weil quadratic form W(f,f)=Σ_ρ|f̂(ρ)|²
restricted to span{V_n:|n|≤N}.  Hence:

    [ τ_N ⪰ 0 for all N, with the spaces exhausting the admissible class
      as N→∞ and λ ranging suitably ]   ⟺   RH.

Self-adjointness of D_log^{λ,N} is automatic (real symmetric τ ⇒ real
spectrum); it is **not** the content.  The content is the *uniform,
unconditional* positivity λ_min(τ_N) ≥ 0.

## Where the irreducible gap is — stated honestly

Proving λ_min(τ_N) ≥ 0 for all N from the arithmetic side (A+B−P, no
zeros) **is** proving RH; it is not a step before RH.  Every reformulation
collapses to the same wall:

- de Branges route: ξ̂'s zeros are provably real ⟺ its structure function
  E is Hermite–Biehler ⟺ the reproducing kernel (= τ_N) is PSD.  Same wall.
- Li/Keiper route: RH ⟺ λ_n>0 ∀n, with λ_n having a zero-free
  Bombieri–Lagarias (Λ + digamma) form.  λ_n>0 ∀n is the same positivity
  on a different basis.

No computation over finitely many N can close this; the literature
(Connes, Bombieri, Li, Lagarias) has not closed it either.  Numerically
re-confirming an equivalence is worthless as proof.

## Net of this session

Track 2 moved from "broken, circular (e_hat·e_hatᵀ feeds zeros in),
eps→−∞, 2π/L lattice" to **a correct, non-circular, RH-equivalent finite
construction with sharp critical-line-specific positivity evidence and
verified prediction of unseen zeros.**  This does not prove RH.  It does
give a single, clean, verified analytic target — *uniform* positivity of
the explicit-formula Gram τ_N — and rules out the prior circular path.

## Only genuinely forward sub-problem

Attack the **uniform lower bound**: find an N-independent decomposition
τ_N = Φ_N* M Φ_N with M a fixed positive operator (the digamma symbol's
Hadamard/de Branges kernel) and Φ_N the explicit-formula evaluation map,
OR an unconditional Hermite–Biehler proof for ξ̂.  Both are RH-hard; any
"proof" of either that uses zero locations is circular and must be
rejected.  This is mathematics, not more scripts.
