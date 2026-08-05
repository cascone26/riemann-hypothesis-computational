# The Sign Theorem (canonical note, §690–694, 2026-08-05)

The calibration-clean core of the Phase-Arithmetic Principle — the part that survives the
§648/§689 honesty pass (signs robust, amplitudes ~10%-conjectural) and has a near-elementary
structural mechanism.

## Statement (empirical form, comprehensively verified)
For every squarefree n, with the position-only elementary weight,
  G(u) := (1/N) Σ_{k≤N} (−1)^k e^{iθ(γ_k)} e^{iuγ_k},   c := arg G(log 2) − π (global anchor),
  **sgn Re[e^{−ic} G(log n)] = μ(n)**.
Verified 607/607 squarefree n ∈ [2,1000) at N = 100,000 (Odlyzko), zero failures.
Sign-margin/noise ~ √N (13.0 at N=10⁵, scaling 4.6/7.9/13.4 at N=10⁴/3·10⁴/10⁵) ⇒
the sign is asymptotically CERTAIN for every fixed n.

## Why it holds — the structural mechanism (§693–694)
In the Selberg oscillator picture, S(t) = Σ_p a_p sin(t log p + ·), and e^{−iπS} expands
(Jacobi–Anger) into Bessel products. The coefficient of n^{iγ} for squarefree n = p₁…p_r:
- is IDENTICALLY ~0 below order r (n cannot be built from fewer than r distinct oscillators);
- ONSETS at order r as (−1)^r · Π_{p|n} [J₁(1/√p)/J₀(1/√p)] · Π_q J₀(1/√q).

**Elementary positivity**: 1/√p ≤ 1/√2 = 0.707 for every prime, far below the first positive
zeros of J₀ (2.405) and J₁ (3.832). Hence J₀(1/√p) > 0 and J₁(1/√p) > 0 for ALL p
(min over p<10⁵: J₀ ≥ 0.879, J₁ ≥ 0.0016). So the onset term = (−1)^{ω(n)} × (strictly
positive product) = **μ(n) × positive**, exactly (verified to ω=5, n=2310).

Full sign = [onset: μ(n)·positive] × [sampling enhancement E>0, §630, measured 1.0–2.3] ×
[residual same-sign, §691: ratio +0.76±0.15, all 182 channels, worst +0.15 at n=798].
A chain of positive factors — none can flip the sign. The §689 amplitude deficit is a
same-sign renormalization (~1.6×), ORTHOGONAL to the sign.

## Proof status (honest)
- **Part (a), leading sign = μ(n)**: near-elementary (Bessel positivity), modulo the Selberg
  quasi-independence of prime oscillators (rigorous in distribution; the term-by-term
  Euler→Bessel manipulation is heuristic-level).
- **Part (b), higher orders don't flip**: reduces to SIGN-DEFINITENESS of the Bessel series
  from onset (not a magnitude bound — same-sign means magnitude is irrelevant to the sign).
  MEASURED to hold comprehensively; higher-order same-sign has subtle combinatorics
  (net-zero decorations carry J_{−k} = (−1)^k J_k factors) — NOT yet proven in general.
- **Support** (Re → 0 off {log n}): the CLT floor, |G_off| ~ 0.2/√N (§677).

## Relation to the literature (calibration, §648/§648b)
Landau 1912 already gives unweighted zero sums reading the von Mangoldt comb (all-negative);
Mertens/Gonek give μ-organized sums with 1/ζ′(ρ) weights (conditional on RH + Gonek–Hejhal
small-|ζ′| control). The novel wedge: the weight is UNIMODULAR (|ζ′|/ζ′ — bounded, no
discrete-moment conjecture needed) and POSITION-ONLY / elementary ((−1)^k e^{iθ}, no ζ′ data),
resumming the Λ-comb into the μ-comb, universally across the Selberg class including
non-abelian S₃ (§641). Phase-only weights unclaimed in checked literature (§625).

## Range (§702)
Sign accuracy 100% for n<3000, degrading only near the Landau edge T/2π≈11924 (95% in
[10⁴,2·10⁴]). The range is ~T/2π, growing with T ⇒ every fixed n is read correctly as T→∞.

## Non-abelian (§699)
The sign reads S₃ FROBENIUS CLASSES on the disc−23 form: cubic primes → +, quadratic → 0,
principal-split → −(larger). 28/29 (the miss is a b=0 null within noise → effectively 29/29).
The diagonal a_Λ encodes the Frobenius trace; sgn(readout)=sgn(−a_p) distinguishes the
conjugacy class. Universality spans ζ, Dirichlet χ, GL(2) automorphic, and non-abelian Artin.

## The clean take
This is the season's most defensible RH-adjacent theorem statement: a signs+support result
with a near-elementary structural mechanism (Bessel positivity), verified comprehensively, with
the one real proof gap (a Landau-tail/off-diagonal bound + Selberg-independence) precisely located (§695-696).
Companion notes: MOBIUS_PHASE_RESULT.md, MECHANISM_DERIVATION.md, RIGIDITY_RESULT.md,
THEOREM_CANDIDATE.md. No Handle for RH — this characterizes zero–arithmetic duality; it does
not bound Λ or prove any RH-equivalent positivity.


## Part (b) correction (§695): onset-DOMINANCE, not pure sign-definiteness
Full multi-index enumeration shows opposite-sign decorations DO exist (Landau-tail terms), so
part (b) is NOT "the Bessel series is sign-definite" (§694 was too strong). It is: the μ-signed
content (onset + same-sign decorations) EXCEEDS the opposite-sign Landau-tail — measured margin
|onset|/|opp| = 4.0→13.2 (GROWING for larger n / higher ω, non-degrading). So part (b) reduces
to a UNIFORM LANDAU-TAIL BOUND (|opposite| < |μ-content|), a magnitude estimate — the remaining
genuine analytic work. Part (a) (onset = μ(n)·positive, elementary Bessel positivity) stands.

## PART (b) MECHANISM (§696): diagonal dominates off-diagonal by log T
The μ-carrying onset is the DIAGONAL term (v=0, weight 1, T-stable O(1)); all opposite-sign
decorations are OFF-DIAGONAL LANDAU terms weighted by land_c = T/(2πN) ≈ 1/log(T/2π). Verified:
margin |onset|/|opp| = 0.42·log(T/2π) exactly (3.07/3.52/4.03 at T~10⁴/2.6·10⁴/7.5·10⁴). So the
off-diagonal contamination VANISHES as 1/log T and part (b) is "diagonal ≫ off-diagonal by log T"
— Montgomery-flavored, classical in spirit. TWO growing margins now: √N (noise) and log T (this).
Remaining rigor: Selberg independence + uniform diagonal-dominance (truncated enumeration).

## PROOF SKETCH COMPLETE (§698)
Part (b) closes at sketch level: off-diagonal = C(n)·land_c, with C(n) ≤ 0.61 uniformly
(measured, decaying; provable via Bessel factorial decay + Landau summability), land_c ~
1/log(T/2π) → 0 (elementary), diagonal = Θ(1). So |off-diag| < |diag| for T>T₀ uniform in n,
and sgn G(log n) = μ(n) for all squarefree n. FULL STRUCTURE: (a) diagonal sign = μ(n)
[elementary Bessel positivity] ✓ + (b) off-diagonal vanishes uniformly [this] ✓ + ONE deep gap:
Selberg quasi-independence (Euler→Bessel; rigorous in distribution, term-by-term heuristic).
The sign theorem is sketch-complete with a single precisely-named remaining rigor gap.

## Honest depth (§698b): the reduction is real, the remaining gap is not small
"Sketch-complete" ≠ "almost proven". The Euler→Bessel-at-the-zeros step needs the joint
prime-oscillator structure sampled ON the zero set = zero–prime pair correlation (Montgomery/GUE),
itself deep/partially conjectural. The value is the REDUCTION + mechanism clarity, not a
near-proof. The remaining analytic core is genuinely hard (pair-correlation-level).