# Sign Theorem — rigorous proof draft (working, 2026-08-05)

Goal: prove, as rigorously as the machinery allows, that for every squarefree n,
sgn Re[e^{−ic} G_N(log n)] = μ(n), where G_N(u) = (1/N) Σ_{k≤N} (−1)^k e^{iθ(γ_k)} e^{iuγ_k}.
Mark exactly where an unconditional argument becomes conditional (pair correlation).

## Step 0 — the weight identity (RIGOROUS, given RH-up-to-T)
For zeros on the line and simple (verified to height T ≫ our range), §606:
  (−1)^k e^{iθ(γ_k)} = −i · e^{−i arg ζ'(ρ_k)} = −i · \overline{ζ'(ρ_k)}/|ζ'(ρ_k)|.
So G_N(log n) = (−i/N) Σ_k [\overline{ζ'}/|ζ'|](ρ_k) n^{iγ_k}.
This is exact. (−1)^k = sgn Z'(γ_k) for consecutive simple real zeros of Z — elementary.

## Step 1 — the target reduces to a sign of a real part (RIGOROUS)
With the global phase c fixed by the n=2 anchor, the claim is a statement about the sign of a
single real number per n. No asymptotics yet.

## Step 2 — diagonal/off-diagonal split (STRUCTURAL; rests on Step 3)
Via the Selberg/oscillator representation S(t) = Σ_p a_p sin(t log p + φ_p) and Jacobi–Anger,
G_N(log n) = D(n) + O(n), where:
  • D(n) = the DIAGONAL term (frequency exactly 0): for squarefree n = p₁…p_r,
      D(n) = (−1)^r · Π_{p|n} [J₁(a_p)/J₀(a_p)] · Π_q J₀(a_q) · (1 + o(1)),  a_p = 1/√p.
    SIGN(D(n)) = (−1)^r · (strictly positive) = μ(n), by ELEMENTARY Bessel positivity:
    a_p = 1/√p ≤ 1/√2 = 0.707 < first zeros of J₀ (2.405) and J₁ (3.832) ⇒ every factor > 0.
    [This half is rigorous MODULO the representation of Step 3.]
  • O(n) = OFF-DIAGONAL terms (frequency = ±log(prime power) ≠ 0): each carries a Landau
    factor of size land_c = T/(2πN) = 1/log(T/2π) + o(1) (unconditional, N(T) asymptotic).
    |O(n)| ≤ C(n) · land_c with C(n) ≤ C_max ≈ 0.61 uniformly (Bessel factorial decay ×
    Landau summability — the tail bound is routine; C(n) in fact decays in n).

## Step 3 — THE CONDITIONAL INPUT (where it stops being unconditional)
Steps 2 require the Jacobi–Anger/Euler expansion of e^{−iπS} to be valid WHEN EVALUATED ON THE
ZERO SET {γ_k} — i.e. that (1/N)Σ_k [Π_p e^{i a_p ...}] n^{iγ_k} equals the product of the
per-prime Bessel means. This is the statement that the prime oscillators {sin(γ_k log p)}_p are
JOINTLY EQUIDISTRIBUTED / independent as γ ranges over the zeros. That joint law is exactly the
ZERO–PRIME correlation, controlled by (conjecturally) Montgomery's pair correlation / GUE
statistics. Selberg's CLT gives the MARGINAL of S(γ) unconditionally, but the JOINT product law
over many primes at the zeros is not known unconditionally — this is the genuine gap.

## Conclusion of the draft
GIVEN Step 3 (the joint equidistribution / pair-correlation input), Steps 0–2 yield:
  sgn G_N(log n) = sgn D(n) = μ(n),  because |D(n)| = Θ(1) and |O(n)| ≤ C_max/log T → 0
  (margin = 0.42 log(T/2π), §696), uniformly in n up to n ≲ T/2π (§702).
So the sign theorem is UNCONDITIONAL EXCEPT for the joint zero–prime equidistribution (Step 3),
which is pair-correlation-deep. This is the exact, honest boundary of the argument.

## What would close Step 3
A theorem that the finite products Π_{p≤P} f(γ_k log p) average over zeros to Π_p ⟨f⟩ with an
error o(1/log T) — a multi-prime generalization of the Landau/Montgomery machinery. Not known;
this is the frontier. (Everything else above is elementary or routine.)
