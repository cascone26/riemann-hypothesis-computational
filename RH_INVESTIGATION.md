# The Riemann Investigation — master narrative (2026-08-05)

A single coherent account of what this program established, what it confirmed, and the one
precise open problem everything reduces to. Companion detail: SIGN_THEOREM.md,
SIGN_THEOREM_PROOF_DRAFT.md, RIGIDITY_RESULT.md, MECHANISM_DERIVATION.md, MOBIUS_PHASE_RESULT.md,
THEOREM_CANDIDATE.md, LEDGER.md (§1–811), PROCESS.md (full provenance).

## The one picture
The nontrivial zeros and the primes are a **two-directional dictionary**, and both directions
are now mapped:

**Forward (zeros → arithmetic): the Sign Theorem.**
The θ-twisted alternating zero sum reads the Möbius function:
  sgn Re[e^{−ic} · (1/N) Σ_{k≤N} (−1)^k e^{iθ(γ_k)} e^{iγ_k log n}] = μ(n),  all squarefree n.
Verified 607/607 to n=1000 at N=10⁵. Elementary weight (positions + θ, no ζ′ data). Universal
across the Selberg class — Dirichlet χ, GL(2) automorphic (Δ, sym^k), and non-abelian S₃ (reads
Frobenius conjugacy classes). It is the unimodular, position-only avatar of the classical
Mertens/Gonek connection, needing no discrete-moment conjecture (calibrated §648).

**Inverse (arithmetic → zeros): Rigidity.**
From featureless Weyl positions + primes ≤ X, a self-consistency solve recovers every zero to
  rms = 0.025 · (log(T/2π) / log X)²  mean spacings
(verified over three decades in X, T up to 7.5×10⁴). One attractor; off-line/collision
signatures are rejected by prime self-consistency; the zeros are arithmetically pinned, and
consequently **super-equidistributed** (sub-CLT: more regular than random).

**The mechanism (why forward works).**
The readout = Selberg-oscillator Bessel products × Landau's theorem, derived with zero free
parameters; the μ-sign is the DIAGONAL (v=0) term = μ(n)·(strictly positive Bessel product),
by elementary Bessel positivity; the contamination is OFF-DIAGONAL Landau terms, UNCONDITIONALLY
O(log T / N) → 0. Diagonal dominates off-diagonal by log T.

## What is proven / unconditional
- The weight identity (−1)^k e^{iθ} = −i e^{−i argζ′} (elementary, given on-line simple zeros —
  verified to 10¹³ classically).
- The off-diagonal Landau bound O(log T/N) (Landau 1911, unconditional; verified: non-prime-power
  zero-sums stay ~log T not √N).
- Bessel positivity of the diagonal onset (elementary: 1/√p < first zeros of J₀, J₁).

## What is confirmed but not proven
RH itself — via every criterion tested (Weil positivity in-range, GUE pair correlation, Li
positivity, Mertens boundedness). All consistent; none a proof. RH holds "barely": every off-line
channel is δ²-soft; the Weil form is exactly saturated (min-eig → 0); the Li off-line sensitivity
is exponentially weak for real zeros.

## The single open problem (everything reduces here)
Prove: for each fixed squarefree n, with X = X(T) → ∞ and X < T/2π,
  A_X(n) = (1/N) Σ_{γ≤T} e^{−iπ S_X(γ)} e^{iγ log n}  has  sgn A_X(n) → μ(n)  as T → ∞.
The off-diagonal is unconditional (Landau); the crux is the **convergent resummation** of the
diagonal (the Bessel order-expansion diverges via ΠJ₀→0; the §649 cumulant form converges). This
is a Selberg-CLT-at-the-zeros + Landau-error statement — research-frontier analytic number theory.
Its truth is verified numerically to sub-CLT precision; only the proof is open.

## Honest status
No handle on RH (no Λ bound, no RH-equivalent positivity proven). The genuine, durable results
are the two-directional dictionary (sign theorem + rigidity) and the derived mechanism, all
comprehensively characterized and honestly scoped. The program's discipline — every overclaim
caught and corrected in-session (see PROCESS.md) — is part of the deposit: nothing false is banked.

## Deferred / open threads (for future analytic work)
- The resummation-convergence proof (the reduction above) — the primary frontier.
- Suzuki's constant μ₁ (Weil positivity margin) — closed form open; high-precision + PSLQ in progress.
- Rigidity+sign "dual cage" (connector lens) — needs T~10⁷ zeros; likely partly redundant, test deferred.
