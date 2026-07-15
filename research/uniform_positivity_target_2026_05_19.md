# The sharpened analytic target — Track 2, 2026-05-19

Purpose: reduce "uniform positivity of τ_N" to the single unconditional
inequality whose proof *is* RH, stated as sharply as it can be, so the
target is unambiguous and no effort is wasted re-deriving the equivalence.

## Setup (fixed, non-circular)

For test functions g in the admissible class (even, Paley–Wiener, the
finite span span{V_n:|n|≤N} as N→∞ and L=log λ²→∞ exhausts it), the
Weil explicit formula is the *identity*

    W(g) := Σ_ρ ĝ(ρ)              (sum over nontrivial zeros ρ)
          = A(g)  +  B(g)  −  P(g)

with, all three terms zero-free / zero-independent:

  A(g) = (1/2π) ∫_ℝ Ω(r) |ĝ(½+ir)|² dr ,  Ω(r)=Re ψ(¼+ir/2) − log π
  B(g) = pole term at s=0,1   (the W_{0,2} block; rank-bounded)
  P(g) = Σ_{k≥2} Λ(k) k^{−1/2} ( ĝ(log k) + ĝ(−log k) )

τ_N is the Gram matrix of W restricted to span{V_n}. Weil/Bombieri:

    W(g) ≥ 0 for all admissible g   ⟺   RH.

## The single target inequality

All of Track 2 now reduces to exactly this, with **no zero ever entering**:

    (★)   A(g) + B(g)  ≥  P(g)      for every admissible real even g,
                                     with a margin uniform in the
                                     bandwidth (i.e. in N, L).

(★) is not a step toward RH. (★) **is** RH, rewritten so that the only
objects appearing are: the digamma symbol Ω (explicit, zero-free,
Ω(r) = log(r/2)+O(1) for large r, Ω(0) ≈ −5.37 — *sign-indefinite*), the
von Mangoldt weights, and the rank-2 pole block. The difficulty is now
isolated to one place and one place only:

  A(g) is a signed quadratic form. Its symbol Ω is negative on a
  bounded r-interval around 0 and grows like log r outside it. So A is
  **not** coercive by itself; positivity of W is a non-trivial
  cancellation in which the log-growth of Ω must dominate the prime sum
  P after subtracting the bounded negative core and adding the pole
  block. Any factorization τ_N = Φ_N* M Φ_N with M ⪰ 0 fixed must
  account for the indefinite core of Ω — that is the entire content.

## Why every reformulation lands on (★)

- de Branges: build E from ξ; H(E) reproducing kernel = limit of τ_N;
  E Hermite–Biehler ⟺ kernel PSD ⟺ (★). Equivalent, not weaker.
- Li/Keiper: λ_n = Σ_ρ[1−(1−1/ρ)^n] ≥ 0 ∀n. The Bombieri–Lagarias
  arithmetic form splits λ_n into a digamma (archimedean) part and a
  Λ(k)-prime part; λ_n ≥ 0 ∀n is (★) diagonalized in the Li basis
  instead of the V_n basis. Same form, same wall.
- Jensen/Pólya (Track 3): hyperbolicity of all Jensen polynomials of ξ
  ⟺ ξ ∈ Laguerre–Pólya ⟺ (★). Griffin–Ono–Rolen–Zagier prove it for
  each fixed degree d only asymptotically in the shift — partial, and
  still the same wall for the full statement.

These are genuinely the same inequality in three coordinate systems.
There is no reduction of (★) to anything strictly easier in the
literature; producing one would itself be the proof.

## What is and isn't progress from here

NOT progress (explicitly rejected, circular or vacuous):
  - any numerical re-confirmation over finitely many N or finitely many λ_n;
  - any "factorization" or "Hermite–Biehler proof" that uses zero
    locations, ξ's known zero-free regions tied to RH, or assumes (★)
    on a subspace;
  - re-deriving the equivalence in a fourth basis.

Genuine progress would be exactly one of:
  (i)  an unconditional, bandwidth-uniform coercivity bound
       A(g) + B(g) − P(g) ≥ c·‖g‖² with c independent of N — i.e. a
       direct positive-definite lower bound on τ_N's smallest eigenvalue
       proven from Ω's analytic structure (its Hadamard product over the
       trivial zeros) plus an unconditional prime-sum upper bound; or
  (ii) an unconditional Hermite–Biehler proof for the secular function
       ξ̂ that does not invoke the location of any ρ.

Both are open problems of exactly RH's difficulty. That is the honest
state: Track 2's *mechanism* is fully built and verified non-circularly;
the *bound* (★) is RH and is not closable by computation or by
re-coordinatizing.

## Disposition

Track 2 numerical/structural phase: CLOSED, verified, documented
(track2_status_2026_05_19.md, connes_decisive_v3.* , control + audit).
The remaining object is (★). Continuing to attack (★) by reformulation
is provably going in circles. The only non-circular avenues are (i)/(ii)
above, which are research-frontier analysis, not scripting — and any
apparent "proof" of them must be checked against the circularity tests
in this note before it is believed.
