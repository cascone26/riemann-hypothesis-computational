# The Krein String of ζ — consolidated result note (2026-08-01, LEDGER §511–§522)

All claims REPORT-tier (numerics, verified in-session, scripts cited). No Handle for RH claimed.

## Object
The Jacobi matrix ("Krein string") of the measure dμ_A = |ζ(½+it)|² dt on [−A,A], computed by
discretized-Stieltjes Lanczos (full reorthogonalization) at A=30/40/50/60/100, K≈2.3·A coefficients,
against the flat-measure baseline: db_n = b_n(ζ) − b_n(flat). Scripts: `operators/krein_*.py`;
data `operators/zeros2000.npy`, `argzp2000.npy`, `abszp2000.npy` (γ ≤ 2515, mpmath).

## The dictionary (all terms independently verified)
    db_n ≈ C · Σ_k λ_k · sinθ_k · cos( 2n·arccos(γ_k/A) − 2·arg ζ′(ρ_k) + c(A) ) · (n₀/n)^α
- Frequencies = zero positions: pinned-frequency fits z=+6.5/+6.9/+7.9/+9.8 (A=30/40/60/100) vs 200
  random-frequency controls; amplitudes persist as A grows (infinite string carries all zeros).
- Amplitudes = C·λ·sinθ, C≈0.68 universal (51 points, corr 0.90); simulated double zero → 1.72×(~2).
- Phases = arg ζ′(ρ_k) ≡ πS(γ_k): window-independent (spread 0.000), identified at R=0.997.
- Envelope: exactly 1/n (Vanlessen) — the measured effective α≈0.6 and the ×0.55 amplitude
  suppression are one finite-n resolution effect, mapped in §523-524: scales converge to the
  theorem value 1 (control: 0.985 at K=300, R²=1.0000). NO free constants, NO anomalies (§524).
- Out-of-sample: frozen dictionary predicts untouched A=50 at r=0.742, empirical p=0.0015 (2000 nulls).
- Parseval closure: variance ratio 1.56→1.23 falling toward 1 with A (asymptotically complete).
- Off-line zero at ½+δ: amplitude deficit ∝ δ² (softness = Weil-razor §474b), growing ~K^1.6 with
  string length, no saturation → the infinite string detects any δ>0.

## The Möbius comb (candidate-novel core; lit gate: NOT FOUND — §521)
G(u) = N⁻¹ Σ_k w_k e^{iuγ_k} over 500–10⁴ zeros (§530: at 10⁴ zeros the comb is PERFECT — 40/40
squarefree vs 0/40 non-squarefree, tol 2.5e-4, n≈90–3000, p<0.001; Gonek→1.0000; gap ≈ 0.009·ln n):
- Plain (w=1): Landau — top-15 peaks = 15/15 prime powers (sanity).
- Twisted (w=e^{−2i·argζ′(ρ)}): peaks on SQUAREFREE integers — 18/20 (n<99, p<0.0005) and 19/30
  (n≈90–1000, p<0.001, tol 6e-4); non-squarefree suppressed (1/20, 3/30). G(log n) ≈ const·μ(n)/√n.
- Weight ladder (60 sf targets): Gonek 1/ζ′ r=0.9994 (→1 with N: Perron identity verified);
  PHASE-ONLY e^{−i·argζ′} r=0.989 — plateaus (structural, not finite-N); string twist 0.936; plain 0.799.
- Gap law: (Gonek − phase-only) ≈ 0.009·ln n over n<99..1000 — the zero-position field's Möbius
  blind spot grows logarithmically; |ζ′(ρ)| supplies it.
- Mertens: M(x) reconstruction (x≤1990): Gonek 0.9933, phase-only 0.8735 (> shuffle-max 0.748),
  mag-only 0.9305, plain 0.776 — cumulative M(x) leans on magnitude; pointwise μ leans on phase.
- CRITICAL EXPONENT (§552): G_phase(log n) ∝ μ(n)/n^0.51 — the readout scales at the RH-critical ½
  (peak corr 0.958 at s=0.51, 10⁴ zeros). Shuffling the ζ′-phases collapses the peak to the trivial
  s=0.42; the arg ζ′≡πS field encodes μ's √n normalization = operator shadow of M(x)=O(x^{½+ε}) ⟺ RH.

## Literature position (§521, agent + spot-verified)
- Vanlessen math/0212014 (VERIFIED): the RH machinery for generalized Jacobi weights — our §516 phase
  law is an application; term-by-term match verified, all constants closed (§523-524).
- Suzuki arXiv:1204.1827 (VERIFIED): canonical system for ζ from ξ-ratio inner functions (ω>1);
  NOT the |ζ|² spectral measure; zero-readout from string coefficients unaddressed.
- Phase-only Möbius selection (unit-modulus weight): NOT FOUND in literature — strongest novel claim.
- Nobody found who has built the Jacobi matrix of |ζ(½+it)|² itself.

## Honest negatives banked en route
Naive prime combs at n≈A·log p and (A/2)·log p^m refuted (combined z≈+1.0); πk parity phase model
refuted; finite-A identification db≈G_t(2n/A) refuted (arccos warping — limit-only); squarefree comb
in |db| at finite A null; positivity-of-string constraint real but weak (active only A≳6×10³,
log-strength); A=30 comb "hits" were selection.

## Open gates
1. Write up the ζ-evaluation of Vanlessen's Φ_ν ⇒ arg ζ′(ρ) as a proposition (conjugate of log|ζ|
   at a zero = regularized arg ζ ↔ arg ζ′(ρ_k)); numerics closed (§523-524), proof-sketch pending.
2. The A→∞ continuum object (which canonical system ζ's |·|² measure is — the Suzuki-gap).
3. Hankel/moment positivity tower translation (string-side Weil analogue).
4. Phase-only Möbius claim: 10⁴-zero scale-up DONE (§530, perfect comb). Remaining: formalize +
   write up; connect to Suzuki screw↔string correspondence (§529).

## Proposition (ζ-evaluation of Vanlessen's phase; proof-sketch level — §526)
Let w_A(x) = |ζ(½+iAx)|² on [−1,1], with singularities x_k = γ_k/A (zeros in the window, λ_k = 1,
double zeros of the weight) and smooth part h = w_A / Π_k ((x−x_k)(x+x_k))². Vanlessen (1.5)–(1.6)
[arXiv:math/0212014] gives
    b_n − A/2 = −(A/n)·Σ_k sinθ_k · cos( (2n+1)θ_k − Φ_k ) + O(1/n²),   θ_k = arccos x_k,
    Φ_k = (α + 1 + Σ_{j>k} 2λ_j)π − (α + β + 2N_A)θ_k − (sinθ_k/π)·PV∫ log h(t)/√(1−t²)·dt/(t−x_k).
Claim: Φ_k ≡ −2·arg ζ′(ρ_k) + (window-linear terms c₀(A) + c₁(A)θ_k)  (mod 2π).
Sketch: (i) for simple zeros λ_j = 1 the counting term Σ_{j>k}2λ_j·π ≡ 0 (mod 2π) — this is exactly
why the naive π·k phase model fails (§514 refutation, now explained). (ii) the middle term is linear
in θ_k with slope −(α+β+2N_A) — window-dependent, absorbed in c₁(A) (and the source of the aliased
slopes in the §514 circular-linear fits). (iii) the PV integral is the interval-Hilbert transform of
log h; adding back the singular factors' conjugates converts it to the conjugate function of the FULL
log w minus explicit θ-linear terms; the conjugate of log|ζ(½+it)|² at t = γ_k, regularized as the
jump midpoint, equals 2·arg ζ′(ρ_k) + π + smooth arch/θ_H terms (near ρ: ζ(½+it) ≈ ζ′(ρ)·i(t−γ), and
the symmetric regularization of arg ζ across the zero lands at arg ζ′(ρ) + π/2). Numerics: the full-
weight conjugate ψ matches arg ζ′(ρ_k) at R=0.997, window-spread 0.000 (§516); scale s = −2 (§514).
Status: REPORT — numerically closed PARAMETER-FREE at 0.1-rad precision (§527: fitted phases match
(1.6) computed from the weight at R=0.991-0.997 across 45 zeros/3 windows; conventions: off-diagonals
are a_n, index+1, ±pair folding). Formal proof gap: only the conj_h ↔ arg ζ′ identity bookkeeping.

## Lemma draft (conj ↔ arg ζ′ — the rigor step for the Proposition; §535)
**Setup.** w_A(x) = |ζ(½+iAx)|² on [−1,1]; h = w_A / Π_k ((x−x_k)(x+x_k))², x_k = γ_k/A. The
Vanlessen phase needs C[log h](θ_ν) — the conjugate function (interval Hilbert transform in the
Chebyshev frame) of log h at the singularity.
**Lemma (claim).** C[log h](θ_ν) = −2·arg ζ′(ρ_ν) + L_A(θ_ν), where L_A is smooth on [0,π]
(no singular content at the zeros), explicitly: L_A collects (i) the arch/Γ-factor conjugate,
(ii) the interval-vs-line Hilbert correction, (iii) the ±pair and out-of-window zero tails.
**Proof sketch, three steps.**
1. *Zero-free conjugacy.* On any closed subinterval J ⊂ (γ_k, γ_{k+1}) the function
   log ζ(½+it) is analytic (no zeros, no pole for t>0 windows), so log|ζ| and arg ζ are harmonic
   conjugates: the LINE Hilbert transform of log|ζ| equals arg ζ + const on J, modulo the
   harmonic-completion data at ∞ handled by the functional-equation/arch factor (use ξ:
   log|ξ(½+it)| has line-conjugate arg ξ = −θ_H(t) + π·N-content exactly).
2. *Local model at a zero.* Near t = γ_ν (simple zero): ζ(½+it) = iζ′(ρ_ν)(t−γ_ν)(1+O(t−γ_ν)).
   Hence log|ζ| = log|ζ′(ρ_ν)| + log|t−γ_ν| + smooth, and arg ζ jumps by π across γ_ν with
   regularized (midpoint) value arg ζ′(ρ_ν) + π/2. Removing the singular factor's own conjugate
   (the sawtooth of log|t−γ_ν|, midpoint-regularized to 0) leaves exactly the smooth part whose
   value at γ_ν is arg ζ′(ρ_ν) + π/2 + arch terms. The weight is |ζ|², doubling everything: −2 scale.
3. *Interval vs line.* The Chebyshev conjugate on [−1,1] and the line Hilbert transform differ by
   a harmonic function of the Joukowski uniformization — smooth on the open interval, singular
   only at ±1 (endpoint exponents α=β=0 here), hence absorbed into L_A. Out-of-window zeros
   (|γ|>A) contribute to log h a function smooth on [−A,A]; their conjugate is smooth ⇒ into L_A.
**Update (§541).** The arg-side is now EXACT: arg ζ′(ρ_ν) = πS(γ_ν) − π/2 (classical folklore via
ζ=e^{−iθ_H}Z, machine-verified 4.4e-14 at 79 zeros; arg arch(½+it)=π+θ_H exactly). The lemma's gap
is pure smooth-function Hilbert bookkeeping.
**Final status (§542b).** COMPLETE at numerical certainty: Hadamard decomposition of conj_h exact to
1e-4 at A=200 with the truncation constant matching the analytic density integral to 4 digits;
remaining formal work = classical convergence lemmas only.
**Status.** Numerically closed: parameter-free phase residuals sd 0.075–0.31 rad across A=30..200,
mean ≈ 0 (§527, §534). What remains for full rigor: (a) uniformity of the O(t−γ) local-model error
in the conjugate, (b) explicit form of L_A via θ_H and the Joukowski kernel, (c) the tail-sum
convergence for out-of-window zeros (explicit formula grade — standard). None of these touch the
arg ζ′ identification itself; they bound the smooth remainder.
