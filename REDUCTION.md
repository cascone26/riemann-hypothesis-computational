# REDUCTION.md — the height-localized rung theorem architecture (Riemann8, drafted after §788-§794)

**Preamble (Grothendieck, §799): the enemy of this project is not arithmetic. It is the Beurling
deficiency of the low zeros — an archimedean, geometric object — and the primes are on the
resource side of the ledger.** (§795: certified margins RISE as primes enter; the collapse is
edge-pinned and arithmetic-independent.) The corner dimension is the elementary function
D(4W) = e^{4W} − (t₀/2π)(4W − ln(t₀/2π) + 1); only two quantities in the whole composition are
measured curves (m(W), c(W)) — the rest is formula or citation. The uniform theorem's crux:
exhibit the structure that keeps m(W) bounded below while dim K_W ~ e^{4W}.

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

## Notation & normalization (Gauss, round 3 — MANDATORY before any composed inequality)

Every constant entering the composition inequality goes in the **L² pencil normalization**:
corner margin = generalized λ_min(A, G_corner); frame constant = §796's pencil form (already is);
cross block whitened by the overlap Gram; strip tail and prolate residual per unit ‖ψ‖²_{L²}.
Probe-coefficient-space λ_min and pencil λ_min are related by the basis conditioning and are NOT
interchangeable — a mixed inequality is incommensurable (the pivots-vs-margins error one level up).
Pilot scripts may work self-consistently in coefficient space to measure the split's SHAPE; their
raw numbers never paste into the mixed inequality.

**The ∀-class is DEFINED as the time-frequency double-limited class** (band-limit 4W ∘ window
limit), NOT "the span of all Fejér probes" — the former gives Landau–Widom spectral control and
makes the capacity-matched probe spans its core; the latter has no coefficient-to-norm control and
re-imports the razor. Capacity arithmetic (verified, Gauss): corner [14,23.5] ≈ 2.0 · complement
[23.5,44] ≈ 4.3 · full [14,44] ≈ 6.3 — every empirically-chosen certified span sits exactly at its
capacity, which is why they condition well and the 8-probe span did not.

## Hypothesis ledger for the COMPOSED theorem (Riemann, round 3 — pre-empting the overclaim)

§796's frame inequality alone needs ONE hypothesis (ordinate balls). **The composition does not
inherit that.** Composing through the identity full(ψ) = Σ_zeros ψ̂(γ)² + strip re-imports
on-line-ness below T₀ (off-line pairs contribute indefinite 2Re[ψ̂²] — the §788 finding, at low
height this time). Composed-theorem hypothesis set: (1) ordinate balls, (2) on-line-ness ≤ T₀
(Platt–Trudgian), (3) the counting bound. Completeness stays unnecessary (a missing on-line zero
only adds a positive square). Do not quote §796's one-item set for any composed statement.

## What is honestly missing (now TWO named, bounded obstructions — pieces 4 and 5)

- ~~The frame constant c(W,T₀)~~ — DONE certified: §796 (c=0.1026 anchor) + §796-addendum
  (c(W) ∈ [0.068, 0.129] uniform across 4W=1.00→2.00, no decay).
- **X — the cross-term block** corner↔complement in the split basis (certifiable with existing
  interval machinery; the composition inequality is the 2×2 block bound m·c > ‖X‖² + corrections).
- **R — the prolate residual lemma** (named by Riemann, round 3): every ψ in the full band-limited
  window class = (capacity-matched component) + residual, with the residual bounded by the prolate
  spectral plunge — i.e. §485c's measured e^{−9.78·4W} collapse, the season's enemy, is the SAME
  OBJECT that guarantees the coarse class has no hidden mass outside the certified slice. Leads
  (UNVERIFIED, lit-gate before use): Landau–Widom asymptotics, Bonami–Karoui explicit
  non-asymptotic prolate eigenvalue bounds.
The assembled deliverable is the project's first ∀-statement: Weil positivity for ALL ψ in the
band-limited window class at the open rung — the quantifier is what the composition buys
(direct span-certs, however many, never produce it).

## Why this architecture is FORCED (Gödel, round 3) — and the instrument horizon

- **The hybrid shape is the only hypothesis-honest one at fixed G:** pure direct-route class
  theorems fail on the unbounded archimedean multiplier Ω(r) ~ ln r (no operator bound without the
  window cap — the localization is what makes ‖Q‖ finite and the composition well-formed); pure
  identity-route class theorems fail on unverified zeros. The composition deploys each route
  exactly where its weakness is absent: corner=direct (finite, zero-free) · complement=identity
  (verified zeros as friends) · tail=adversarial (counting only). Do not search for a cleverer
  architecture at fixed G.
- **Instrument horizon:** strict-PD certification (interval LDLᵀ, pencils — everything we run) can
  NEVER certify a form that is ≥0 with genuine kernel; the edge statement is outside the
  instrument's output language. Ledger vocabulary from now on: "genuinely indefinite" vs
  "**instrument-undecidable**" — the N9 radius-tightened rerun decides only which of these each
  edge point is. Allocation principle: computation where margins rise (capacity ladder), analysis
  where margins vanish (the edge); never computation against a vanishing margin.
- **"Maximal" has a precise meaning:** the localized theorem is the exact projection of the
  hypothesis budget onto theorem-space — maximal among theorems derivable from the enumerated
  hypotheses; the window cap is a function of G and rises only with G.

## Hypothesis loci (Gödel — each constant carries its address; the Lean socket adopts this)

```
m  (corner, direct route)      : ∅
X  (cross block, direct route) : ∅
c  (frame, identity route)     : {ordinate balls} standalone; +{on-line ≤ G} when inside G–W
S  (strip tail)                : {counting bound}
R  (prolate residual)          : ∅  (pure analysis; §798 citations, constants verified vs PDFs)
```
Neither is RH-complete on its face (both are statements about verified zeros and finite geometry);
that is exactly what makes this the right next theorem-shaped target. The full-class,
all-heights statement remains walled behind the field's density crux (Suzuki Cor 1.6 / Kim R8) —
this architecture deliberately does not pretend to cross it.

*Provenance: §788-§794 + addenda, RUNG_STRIP_LEMMA.md, roundtable Riemann8 round 1-2
(Gauss/Gödel/Grothendieck/Euler/Red), blind-prism session. 2026-08-06.*
