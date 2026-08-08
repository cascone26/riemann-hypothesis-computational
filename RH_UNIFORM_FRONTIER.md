# The RH-Uniform Frontier — honest synthesis of Riemann8 §808-821

**One-paragraph state:** We have a *certified + Lean-verified* Weil-positivity theorem for a described
infinite class at one rung (§808), and a thorough, mostly-confirmed *characterization* of what the
uniform (all-widths) version requires. The uniform theorem is exactly Weil's criterion (⟺ RH) and
remains open — as it does for the whole field, including the leading spectral program (CCM 2025, §820).
What Riemann8 adds is a precise, measured map of *where* and *why* the difficulty sits, plus certified
sub-pieces. Six proposed scalar "mechanisms" were tested and refuted (below) — the durable results are
structural, not scalar.

## What is PROVEN / CERTIFIED (Handle-tier)
- **§808 the anchor ∀-theorem** (certified + Lean-kernel-verified, trusted axioms): Weil Q(ψ)≥0 for
  every even ψ∈PW_{0.66} at the open rung with weighted off-window tail T(ψ)≤1e-4 (certified Q≥0.02).
  First fully-certified ∀-positivity theorem for the Weil form. All constants arb-certified on one
  self-consistent basis; closing inequality in arb AND Lean.
- **§796 frame constant above t***: Σ_{γ≤800}ψ̂(γ)²≥0.10·‖ψ‖² certified on the above-t* class (ZERO
  side — conditional on the zeros used; see caveat).
- Prior certified apparatus: §790 zero-free ladder to 4W=2.0, §791 strip lemma, §792 Lean Gram certs.

## What is CHARACTERIZED (Report-tier, tested + explained)
- **§812** difficulty LOCALIZES to a boundary shell of concentration modes (m_E·K W-stable for the
  core; margin crashes at the boundary).
- **§813/§813b** the crash is CROSS-COUPLING (Q positive per-mode), and the coupling is ARITHMETIC:
  the archimedean form is a PD floor, the PRIME form is INDEFINITE, its positive eigenvalues eat the
  floor on boundary-coupled directions.
- **§816** the archimedean floor is W-INVARIANT (~0.153, CV 5% over 10 widths), EXPLAINED by the fixed
  window (Ω over [14,44] is W-independent). The prime obstruction GROWS with support.
- **§818/§821** the SHARP mechanism, Beurling-derived + confirmed: Q=2Σ_γφ̂(γ)² is a SAMPLING form;
  t*=2πe^{4W} is the Nyquist threshold; zeros oversample above t* (frame-positive) and undersample
  below (the §794 corner); the Weil margin collapses EXACTLY as the window's zeros fall below t*
  (measured: 2/8 below→0.093, 8/8 below→0). Q stays ≥0 (→0⁺) because the window zeros are REAL; an
  off-line zero → negative (§780 detector).

## §832 THE KNIFE-EDGE, cleanly isolated (new — the sharpest finite RH-shadow so far)
The Weil margin on the window [14,44] is CERTIFIED robust (≥0.10, arb) while sampling FEWER modes than
the window has zeros, and thins to 0⁺ EXACTLY as N → #window-zeros (=8), with the Gram matrix
WELL-CONDITIONED throughout (cond(G)≤4 at the collapse — so it is NOT a conditioning artifact, unlike the
§829 over-capacity collapse). The margin goes to 0⁺, never negative — precisely RH's prediction (marginal
positivity at critical sampling; an off-line zero drives it negative). Certified ladder: m(4W)=0.151,
0.145, 0.103, 5.6e-5 at 4W=1.30/1.60/1.90/2.20 (N≈0.78·N_cap). The finite razor is set by the WINDOW's
zero-count, not the frame capacity — the clean finite-dimensional shadow of RH. See LEDGER §832.

## The HONEST open core (= RH)
Prove **arch(ψ)+pole(ψ)−prime(ψ) ≥ 0 for all band-limited ψ, uniformly in width** — the PRIME side,
with NO assumption on zeros. This is Weil's criterion. The margin is razor-thin (→0) on the sub-t*
corner directions; Beurling density (§821) explains WHY (undersampling) and WHERE (below t*); but the
proof does not get to assume the zeros, so the illuminating zero-side split (§821) is NOT an
unconditional reduction (§821-CAVEAT). The unconditional target is monolithic prime-side positivity.

## §821 status, precisely (after self-caveat, then §829 refinement)
CONFIRMED and durable: the margin-tracks-t* prediction (mechanistic understanding of the collapse).
WITHDRAWN: the "unconditionally reduces RH to the finite corner" framing — the above-t* frame half is
zero-side/conditional. Use §821 as understanding, not as a proof reduction.

**§829 REFINEMENT (clean capacity-matched probe, cert_beurling_clean.py).** The "margin collapses to
0 as the window's zeros fall below t*" numbers (§821: 2/8 below→0.093, 8/8→0) were CONDITIONING-
CONFLATED — measured on dense/near-capacity spans where the t* sweep and the frame-capacity knee move
together. With mode count held FIXED below capacity and G-conditioning reported every row, sweeping 4W
so t* crosses the whole window does NOT collapse the margin: it stays robustly positive (0.15–0.19).
The dramatic collapses coincide EXACTLY with λmin(G)→0 at the frame capacity N_cap=30·2W/π (critical
sampling), independent of t*. So §821 is right about WHERE the razor lives (critical sampling) but the
collapse-to-zero was the capacity degeneracy, not the pure t* crossing. Honest core, sharpened: below
capacity margin is robustly positive across the crossing [float, solid]; at critical sampling the margin
is instrument-undecidable (float=noise ~1e-11; arb radii swallowed it in §823) — resolving whether the
critically-sampled margin stays ≥0 uniformly IS the hard core = RH.

## Connection to the field (§820, CCM 2511.22755)
The leading spectral program's key basis = the PROLATE spheroidal functions = our concentration
eigenspace. Their conjectured μ_λ→0 = our measured Q_floor→0 (we measure →0⁺, slightly stronger).
Their stated remaining obstacle (prolate→Hermite convergence, Lemma 7.2) lives in the same basis.
Neither the arithmetic-surface route (barrier-map B18) nor the CCM spectral route closes the core;
both bottom out where our measurements say the difficulty is. Our apparatus can numerically probe
CCM's open lemma (a concrete, unused bridge).

## The six refuted scalar mechanisms (DO NOT re-propose)
m_E·D=const (§810b), saturation (§811), peak+decline (§811-corr), eigenvector-misalignment (§814/§815),
prime∝prime-sum (§817), "proof-strategy reduction" (§821-caveat). Each proposed, tested, refuted/
tempered. LESSON: Weil positivity is irreducibly multi-dimensional; it does not reduce to a scalar
bound. Progress requires structural/theoretical work, not more scalar float scans.

## Sharpest next targets (for a fresh attack)
1. **Prime-side spectral inequality** (the real prize): bound the prime form's positive eigenvalues by
   the archimedean floor STRUCTURALLY (not scalar) — the intersection-form/Hodge-index object is what
   would do it. Engage de Branges / BRS sampling theory (§793, §818 reframe) for the frame lower bound
   handling the zeros' small gaps (Red §788).
2. **Certify the characterization**: turn §816 (arch-floor W-invariance) and §818/§821 (t*-margin) into
   arb-certified structural facts.
3. **CCM bridge**: numerically test their Lemma 7.2 in our apparatus; connect our Q_floor→0⁺ rate to
   their μ_λ→0 conjecture — a concrete contribution to the leading program's open problem.
4. **Bank the anchor theorem as a writeup** (§808 is publishable-shaped: a certified+Lean-verified
   Weil-positivity ∀-theorem with an honest prior-work section).

*Provenance: LEDGER §808-821 + caveats, PROCESS arc, ANCHOR_THEOREM.md, cert scripts. Riemann8,
2026-08-07. The over-reaches are kept in-ledger deliberately (process provenance) — they show the
propose→test→refute methodology and prevent re-walking dead ends.*
