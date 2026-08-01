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
G(u) = N⁻¹ Σ_k w_k e^{iuγ_k} over 500–2000 zeros:
- Plain (w=1): Landau — top-15 peaks = 15/15 prime powers (sanity).
- Twisted (w=e^{−2i·argζ′(ρ)}): peaks on SQUAREFREE integers — 18/20 (n<99, p<0.0005) and 19/30
  (n≈90–1000, p<0.001, tol 6e-4); non-squarefree suppressed (1/20, 3/30). G(log n) ≈ const·μ(n)/√n.
- Weight ladder (60 sf targets): Gonek 1/ζ′ r=0.9994 (→1 with N: Perron identity verified);
  PHASE-ONLY e^{−i·argζ′} r=0.989 — plateaus (structural, not finite-N); string twist 0.936; plain 0.799.
- Gap law: (Gonek − phase-only) ≈ 0.009·ln n over n<99..1000 — the zero-position field's Möbius
  blind spot grows logarithmically; |ζ′(ρ)| supplies it.
- Mertens: M(x) reconstruction (x≤1990): Gonek 0.9933, phase-only 0.8735 (> shuffle-max 0.748),
  mag-only 0.9305, plain 0.776 — cumulative M(x) leans on magnitude; pointwise μ leans on phase.

## Literature position (§521, agent + spot-verified)
- Vanlessen math/0212014 (VERIFIED): the RH machinery for generalized Jacobi weights — our §516 phase
  law is an application; the α=0.6 anomaly is the open delta vs its 1/n order term.
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
1. Analytic derivation of the dictionary from Vanlessen-type RH asymptotics for the ζ weight —
   must produce the arg ζ′ phase and explain α≈0.6 (collective effect?).
2. The A→∞ continuum object (which canonical system ζ's |·|² measure is — the Suzuki-gap).
3. Hankel/moment positivity tower translation (string-side Weil analogue).
4. Phase-only Möbius claim: extend zeros (10⁴+), formalize as a statement about
   Σ e^{−i·argζ′(ρ)} x^{iγ} vs Σ μ(n)-combs, then write up.
