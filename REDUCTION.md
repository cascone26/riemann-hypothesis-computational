# REDUCTION.md — the height-localized rung theorem architecture (Riemann8, drafted after §788-§794)

**Scope, stated first (Red's §788 constraint, honored in the first paragraph):** this architecture
targets Weil positivity on a **height-localized** restricted class — Fejér probes (quartic decay,
NOT Gaussian; s = center spacing, Gauss's §791-addendum correction), centers ≤ C_max, Fourier
support 4W — using only verified zeros ≤ T₀. It is NOT the full band-limited rung statement: the
frame inequality over finitely many zeros is FALSE on any class of dimension > #zeros (Red's
§788 attack, upheld). The full-class statement needs the uniform-in-W machinery (Grothendieck's
dévissage) plus something genuinely new for the unbounded-height completion.

## The object

Q(ψ) = pole(ψ) + arch(ψ) − prime(ψ) for ψ = Σ vⱼ F_{cⱼ}, F Fejér pairs, cⱼ ∈ [c_min, C_max],
support parameter 4W. Guinand–Weil: Q(ψ) = 2Σ_{γ>0} ψ̂(γ)² under RH; unconditionally, zeros
contribute ψ̂(t−iσ)-terms per §791.

## The three pieces + the corner formula

1. **Corner (below t\* = 2πe^{4W}).** The undersampled region: zero density (1/2π)ln(t/2π) <
   class capacity 4W/2π exactly for t < t\* (Gödel §788, now QUANTITATIVELY confirmed §794-add).
   **Corner dimension formula:** dim ≈ ∫_window (4W/2π − (1/2π)ln(t/2π))₊ dt — the Beurling
   deficiency integral (measured fit: monotone, slope right, constants within ~½; float=Report).
   The corner is GEOMETRIC (density), not arithmetic (§794-add: prime-node reading refuted).
   Handled by: FINITE certificates — the existing zero-free interval-LDLᵀ machinery (§788, §790),
   λ_min-grade via certify_lambda (§791-add-3; pivots certify sign, not margin).

2. **Frame inequality (window above t\*, below T₀).** Σ_{t* < γ ≤ T₀} ψ̂(γ)² ≥ c(W,T₀)‖ψ_above‖²
   on the corner-complement. Above t\* the zeros oversample the class capacity — a genuine frame
   bound is plausible with explicit constants. STATUS: OPEN, the real analytic work. BRS
   (2005.02996) is interpolation NOT frames (§793 downgrade) — their Dirichlet-kernel estimates
   are scaffolding only. The needed statement is finite-dimensional in practice (class dim below
   any fixed T₀ is finite) — candidate route: certified λ_min of the zeros-only sampling Gram on
   the above-t\* subclass (the §794 probe machinery, certified; its measured spectra are already
   comfortably positive above the tiny cluster).

3. **Strip tail (above T₀).** §791 (PROVEN + certified + adversarially verified): arbitrary
   configurations cost ≤ cosh⁴(Wσ)·(on-line envelope), two-sided. Hypothesis set FOUR items
   (§791-add-3). Adversaries only — this piece can never supply positive mass (Red's constraint;
   using unverified zeros as friends re-imports the RH-fragment §788 excised).

**Composition:** Q(ψ) ≥ [corner cert margin on ψ_corner] + [frame constant on ψ_above]
− [cross-terms] − [strip tail S(W,T₀)] — with the cross-term control between the corner projection
and its complement the second open analytic item (the §476g 8-probe near-singularity lives here).

## Hypothesis-set discipline (the §791 pattern, mandatory)

Every certified statement carries its complete enumerated hypothesis set:
(1) zero-list completeness to the cited height (counting-formula check, Gauss),
(2) on-line-ness to that height (Platt–Trudgian),
(3) the classical counting bound (cite Backlund/Trudgian with constants at paper time),
(4) ordinate enclosures (currently mpmath dps-30 balls at 1e-9; hardening path = interval-Newton
    re-enclosure on ξ, queued §791-add-3).
Zero-free-route statements (§788/§790/§792 pattern) carry NONE of these — prefer that route
wherever it reaches.

## What is honestly missing (the two open analytic items)

- The frame constant c(W,T₀) above t\* with explicit constants (piece 2).
- Cross-term control corner↔complement in the composition.
Neither is RH-complete on its face (both are statements about verified zeros and finite geometry);
that is exactly what makes this the right next theorem-shaped target. The full-class,
all-heights statement remains walled behind the field's density crux (Suzuki Cor 1.6 / Kim R8) —
this architecture deliberately does not pretend to cross it.

*Provenance: §788-§794 + addenda, RUNG_STRIP_LEMMA.md, roundtable Riemann8 round 1-2
(Gauss/Gödel/Grothendieck/Euler/Red), blind-prism session. 2026-08-06.*
