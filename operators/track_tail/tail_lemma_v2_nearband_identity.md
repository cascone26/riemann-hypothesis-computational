# Tail Lemma v2 — the near-band cancellation identity (exact), and what remains to prove

**Status (2026-07-24, roundtable session):** the cancellation Rie identified as the crux of the
modes>K₀ tail lemma is now DERIVED IN CLOSED FORM and VERIFIED numerically against the raw kernel
(rel. err 1e-6…1e-4, quadrature-limited). It is an exact identity, not an asymptotic estimate.
Supersedes the numerical content of `tail_delta_k0_lemma.md` (broken draft — see DEPRECATED marker
there; its C_k ≈ 2.2k off-diagonal model is wrong in form). The theorem statement + Schur skeleton
of the old draft (its Part 0) remain valid and are not repeated here.

Verification artifacts: `~/tmp/nearband_identity_check2.py` (to be moved here), measured answer key
`results/tail_rowsum_mode51_out.txt`, matrix `results/Q51_a0549_nl4000_float.npy`.

## Setup

Even-sector Weil block at window parameter a (tightest a = 0.549), b = π/(2a), odd modes k, l.
From `rigor/k0_hp_nl24000.py` (verbatim kernel): for k ≠ l,

    Q[k,l] = (1/A) · [ a_l·J(lb) + a_k·J(kb) ],      J(ω) = ∫₀^{2a} H(t)·sin(ωt) dt

with

    a_l = 2 k l² b / (k²−l²),    a_k = −2 k² l b / (k²−l²)          (resonant, O(k²/|k−l|) near band)
    H(t) = G(t) + w_p·(t−log 2)·𝟙[t>log 2]
    G(t) = −4e^{t/2} − 4e^{−t/2} + 8 − (Ω₀/2)t − Φ₀/4 + Σ_{n≥0} ¼(n+¼)^{−2} e^{−(2n+½)t}

(G is the summed archimedean/Lerch kernel; w_p = log2/√2 the p=2 comb weight. H ∈ W^{1,1}(0,2a);
H′ has an integrable log singularity at t=0, H′(t) = **+½·log t** + O(1) (→ −∞; measured
H′ − ½log t ∈ [1.19, 1.36] on the resolved range — the sign was stated backwards in the first
version of this doc, corrected same day), and a jump w_p at t = log 2.)

## The two structural facts

1. **Commensurability.** 2ab = π, so at the D-term frequencies ω = mb (m odd):
   cos(2a·mb) = cos(mπ) = −1 and sin(2a·mb) = 0. Boundary terms of integration by parts
   evaluate at integer multiples of π — this is exact, by construction of the mode basis.
2. **Coefficient reciprocity.** a_l/(lb) = 2kl/(k²−l²) = −a_k/(kb), exactly.

## The identity

Integrate J by parts once (valid, H ∈ W^{1,1}):

    J(ω) = [H(0) − H(2a)cos(2aω)]/ω + (1/ω)∫₀^{2a} H′(t)cos(ωt) dt
         = [H(0) + H(2a)]/ω + (1/ω)∫₀^{2a} H′(t)cos(ωt) dt        (ω = mb, m odd)

By fact 2 the boundary contributions cancel **exactly**:

    a_l·[H(0)+H(2a)]/(lb) + a_k·[H(0)+H(2a)]/(kb) = 0.

What survives, using cos(lbt) − cos(kbt) = 2 sin((k+l)bt/2)·sin((k−l)bt/2):

    ┌──────────────────────────────────────────────────────────────────────────────┐
    │  Q[k,l] = (4kl / (A(k²−l²))) · ∫₀^{2a} H′(t)·sin((k+l)bt/2)·sin((k−l)bt/2) dt │
    └──────────────────────────────────────────────────────────────────────────────┘

**Exact** (no asymptotics), all quantities explicit. Numerically verified against the direct
computation (Q51 matrix, NL=4000): pairs (15,13),(15,17),(27,29),(39,37),(39,41),(39,7),(51,47),
(51,7),(15,49),(3,51) — rel. err 1.2e-6 … 1.5e-4, limited by trapezoid quadrature of the
log-singular H′, not by the identity.

## Why this kills the old draft's O(k) blowup

The draft bounded |a_l|, |a_k| separately (→ C_k ∝ k, row sums ~73, self-contradicting its own
claimed 0.04). In the identity the resonant coefficient 4kl/(k²−l²) = 4kl/((k−l)(k+l)) is paired
with sin((k−l)bt/2), which (i) vanishes linearly at t = 0, taming the log singularity of H′, and
(ii) is bounded by min(1, |k−l|bt/2), supplying the |k−l| that cancels the resonant denominator.
The fast factor sin((k+l)bt/2) supplies decay in (k+l) via one further IBP / van der Corput.
Heuristic count: (kl/((k−l)(k+l))) · O(|k−l|/(k+l)) = O(kl/(k+l)²) = O(1) uniformly — matching
the measured envelope sup_{k,l} |Q[k,l]|·|k−l| = 1.15 and bounded nearest-neighbor entries.

## Remaining proof obligations (this is the actual analysis to write)

1. **Near-band bound.** From the identity, prove |Q[k,l]| ≤ C₁/|k−l| uniformly in k, with explicit
   C₁: split H′ = (−½ log t) + (bounded, BV part incl. the w_p jump); on the log piece use the
   linear vanishing of sin((k−l)bt/2); on the BV piece IBP the fast factor. Measured target:
   any C₁ ≥ 1.2 is consistent (measured sup = 1.15); aim C₁ ≤ 2.
2. **Infinite row-sum bound.** Σ 1/|k−l| alone log-diverges over the infinite tail; the far band
   needs the extra (k+l)-decay from obligation 1's fast-factor IBP (measured: far entries run
   ~13× below the C₁/|k−l| envelope at distance 38). Target: uniform row sum ≤ S with S < the
   diagonal lower bound. Measured row sums (window 51, non-edge rows k∈[15,49]): ≤ 1.07.
3. **Diagonal lower bound — STRUCTURE NOW DERIVED AND VERIFIED (same-day update).** Exact
   diagonal identity, same family as the off-diagonal one (IBP once; boundary vanishes since
   sin(kπ)=0 and the (2a−t) factor kills t=2a):

       Q[k,k] = −(kb/a) · ∫₀^{2a} H′(t)·(2a−t)·sin(kbt) dt

   Verified vs Q51: rel. err 6.4e-8 … 4.0e-6 (`diag_identity_check.py`,
   `results/diag_identity_check_out.txt`). Splitting H′ = ½log t + φ (φ ∈ BV, jump w_p at log 2)
   gives the two-term law, fit over odd k ∈ [15,51] with max residual **0.0148**:

       Q[k,k] = log(kb) − 1.8527 + 0.3608·cos(kb·log 2 + φ₀) + O(1/k)

   The oscillation amplitude is ANALYTIC and matches to 4 digits: w_p·(2a − log 2)/a = 0.3612
   (the jump at log 2 IBPs to a boundary term at t = log 2 with exactly this weight). Implied
   floor: log(25b) + c₀ − amp = **2.057** for k ≥ 25 (measured min 2.1425 — the analytic floor
   is the worst-phase envelope; actual odd k never hit worst phase in-window). Remaining to
   rigorize: (i) analytic value of c₀ (fitted −1.8527; comes from γ_E + the φ-boundary terms
   of the ½log t and BV pieces — elementary), (ii) explicit constant in the O(1/k) error (TV
   bounds on φ). Both are bookkeeping on the verified identity.

   **SAME-DAY UPDATE 2 — c₀ IS CLOSED FORM: c₀ = −log 2π.** Derivation chain (verified,
   `diag_closed_form_check.py` / `results/diag_closed_form_check_out.txt`):
   - Lerch asymptotic Σ_{n≥0} xⁿ/(n+¼) = −log(1−x) − ψ(¼) − γ_E + o(1) (x→1⁻) gives
     H′(t) = ½log t + φ(0⁺) + O(t·log t) with **φ(0⁺) = ½log 2 + ½ψ(¼) + ½γ_E − Ω₀/2
     = ½(γ_E + log 2π)** — the ψ(¼) in Ω₀ = ψ(¼) − log π cancels exactly. Numeric: 1.2075
     predicted vs 1.190 measured at t = 0.01 (O(t) consistent).
   - Term A (log piece), using ∫₀^X log t·sin(ωt)dt = (1/ω)[−γ_E − log ω + Ci(ωX) − cos(ωX)log X]
     at ωX = kπ: contributes log(kb) + γ_E − Si(kπ)/(kπ) with Si(kπ) = π/2 + 1/(kπ) + O(1/k³)
     (asymptotic, per Ein's scoping note — same applies here).
   - Term B (φ piece): boundary −2φ(0⁺); the jump w_p at log 2 contributes the oscillation
     −(w_p(2a−log 2)/a)·cos(kb·log 2) with amplitude 0.361439.
   - γ_E cancels between A and B:  **Q[k,k] = log(kb/2π) − 0.361439·cos(kb·log 2) − 1/(2k) + O_osc(1/k)**.
   Measured residual·k ∈ [−0.77, −0.09] over odd k ∈ [15,51] — bounded, oscillating around the
   −1/2 predicted mean.

   **SAME-DAY UPDATE 3 — C_d IS CLOSED FORM AT 1/k ORDER, BOTH PARITIES.** The 1/k coefficient:
       R(k)·k = −1/2 − (2w_p/π)·sin(kb·log 2) + O(1/k)
   Sources: −1/2 from −Si(kπ)/(kπ) (Si(kπ) = π/2 ± 1/(kπ), sign by parity — the parity flip is
   O(1/k²) only); the oscillation from the SECOND IBP boundary at log 2 — g′ = φ′(2a−t) − φ
   jumps by −w_p through the −φ term, and ab = π/2 exactly gives coefficient w_p/(ab) = 2w_p/π
   = 0.31203. Envelope: **|R(k)·k| ≤ 1/2 + 2w_p/π + ε_k = 0.81203 + ε_k.**
   Verified BOTH parities (`cd_both_parities_check.py`, `results/cd_both_parities_out.txt`):
   even k = 16/26/40/50 match the model at 1.4e-3–1.4e-2 (k=26: measured −0.8012 vs model
   −0.8005 — envelope nearly attained, so 0.812 is sharp); odd k match within the 1/k² band.
   The SAME leading law log(kb/2π) − 0.3614·cos(kb·log 2) holds for even-k modes (odd-parity
   sector) — independently confirms the sector-independence of the diagonal law at 1/k² precision.
   REMAINING (last line item): the ε_k = O(1/k) bookkeeping constant via one more IBP
   (φ_s ∈ C^∞[0,2a] by the Mellin expansion of the Lerch sum, so V₃ = ∫|g‴_s| is finite);
   measured model-residuals ≤ 0.014 at k ≥ 25 say ε_k is small. **Working constant: C_d = 0.9
   for k ≥ 25, both parities** (margin 0.088 over the near-attained envelope).
   Floor formula for the M* optimizer, analytic at EVERY cut, BOTH sectors:
       **d₀(k) ≥ log(kb/2π) − 0.361439 − 0.9/k    (k ≥ 25, odd or even modes)**
   Structural note (Report-tier, but exact): log(kb/2π) = 2π·N′(kb) — the diagonal of the Weil
   form in this basis equals 2π × the Riemann–von Mangoldt zero-counting density at height kb.
   The archimedean diagonal IS the local zero density; the Γ-factor of the explicit formula
   seen as a matrix entry.
4b. **WARNING for obligation 2′ (norm split) — measured this session, do not inherit the wrong
   inequality.** (a) The claim "modulated Hilbert-type kernels have norm strictly below π" is
   FALSE for sum-modulation in complex form: e^{iθ(k+l)}/(k+l) = diag(e^{iθk})·Cauchy·diag(e^{iθl})
   is a unitary conjugation — norm exactly equal to the Cauchy section (measured 0.8965 = 0.8965,
   window odd 25..3023). The REAL cosine version does drop (measured 0.4507) via Re-part
   interference, but no "linear in θ" law is established — prove it or measure it, don't cite it.
   (b) The β-kernel's exact decomposition is
       2[k·cos(θl+φ) − l·cos(θk+φ)]/(k²−l²)
         = [cos(θl+φ)+cos(θk+φ)]/(k+l)  −  2·sin(θ(k+l)/2+φ)·sin(θ(l−k)/2)/(k−l)
   and the DIFFERENCE-channel (second) piece has in-window norm **1.4869** — larger than the
   Cauchy piece (0.897) and dominant in the norm budget. Crude bound: it is a diagonal-modulated
   Toeplitz with symbol the [−θ,θ] indicator (norm ≤ π) — but crude-π on this piece plus the
   others will not clear d₀; the small measured total far-field norm (0.185, Ein) must come from
   small amplitudes γ_eff/β_eff extracted from ĥ, and THOSE are what obligation 2′ must pin down.
4. **Compose:** δ_tail ≥ d₀ − S > 0 for all odd k > K₀ = 23 (measured floor at K₀=23: 1.22),
   then feed δ_tail into the (valid) Schur skeleton of the old draft with the NL=24000/Arb
   low-block numbers (`results/k0_hp_nl24000_out.txt`, K₀=23: λmin 6.077e-8, residual 6.40e-12).

## Answer key (all measured, on disk)

| quantity | value | source |
|---|---|---|
| envelope sup |Q[k,l]|·|k−l|, k∈[15,51] | 1.15 | Q51_a0549_nl4000_float.npy |
| row sums, non-edge k∈[15,49] | ≤ 1.07 | tail_rowsum_mode51_out.txt |
| min Gershgorin margin, window 51 | 0.94 (k=15) | same |
| margin floor for K₀=23 tail (k≥25) | 1.22 (k=29) | same |
| near-band (|k−l|≤4) share of row sum | 0.46–0.69, median 0.51 (edge row k=51: 0.32, excluded) | this session |
| nearest-neighbor entries | bounded ≤ 0.30, oscillatory, NOT monotone in k (0.27@15, 0.13@27, 0.30@29, 0.05@51-edge) | this session |

Caution recorded for accountability: the earlier "nearest neighbors decrease in k (0.52→0.045)"
read used the k=51 EDGE ROW (half its near band lies outside the computed window) and pairwise
sums; mid-window rows oscillate up to ~0.59 (sum of both neighbors). Boundedness holds; monotone
decrease does not. Do not build the lemma to predict decay in k — build it to predict uniform
boundedness, which is what the identity gives and the data shows.
