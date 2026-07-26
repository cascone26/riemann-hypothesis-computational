# ROUNDTABLE GROUNDING — state as of 2026-07-26 evening (LEDGER through §342)

Paste-ready banked-results block. Everything here is DEPOSITED — restating it scores nothing;
build on it or attack it.

## Banked from Roundtable #4 (overnight 07-26, §334) — do not re-derive
- Every RH-equivalent toured (Weil ⟺ GNS ⟺ reflection-positivity ⟺ Li ⟺ κ_N) is a biconditional provable in a
  weak base theory — CONSERVATIVE paraphrases, no middle term. A real middle term must be non-conservative +
  constructible-without-zeros + manifestly positive (3-clause gate).
- The complementary series is UNITARY ⟹ every positivity frame (W≥0, reflection positivity, definite pairing)
  is blind across the tempered edge λ=¼. The true threshold is an INTEGRABILITY/decay condition (L², no slack).
- Kubota kill: on the modular surface the archimedean scattering determinant is φ(s)=ξ(2s−1)/ξ(2s) — "index at
  ∞" = compute ζ. Circular. General law: every coordinate placing the wall at the CENTER (Re=½) came back as ζ;
  only EDGE (Re=1) coordinates ever delivered theorems.
- Retired for cause: ordinal-of-RH (Specker), III₁ invariant (universal ⟹ discriminates nothing),
  cocycle/coboundary + flow-of-weights (asymptotic), index-at-∞, archimedean motive/Frobenius as primary target,
  positivity as governing frame.
- GATED (t01 prefilter + conflation catch): the table's own "follow" prescriptions — sym⁵/sym⁶ Maass
  functoriality, Guth–Maynard density extension, Sato–Tate compactness-as-signature — all REDUCE-TO-WALL for
  our goal. Also: Kubota's φ=ξ/ξ is the CONTINUOUS (Eisenstein) spectrum; Selberg's λ₁≥¼ is the DISCRETE
  cuspidal spectrum. Selberg's conjecture is RH-analogous, NOT RH-implying — no known bridge to RH for ζ.

## NEW since that table (same day, §335-342) — the day's measured results
All on the validated instrument operators/l2_edge_sweep.py (analytic side pole+arch−prime vs zero side must
close; residual 1e-10 = its own correctness certificate; reproduces every banked §286-287 number).

1. **The dense-regime collapse is a CRITICAL-DENSITY SAMPLING phenomenon, not a sign phenomenon (§335).**
   The Weil form on Gaussian probes stays exactly PSD at all densities; its margin plunges super-exponentially
   to EXACT rank deficiency precisely when probe count exceeds the Riemann–von Mangoldt zero count of the band
   (Nyquist in the zeros), independent of Gaussian width; sign-blind controls (overlap Gram, arch-only,
   smooth-Weyl-density Gram) show NO threshold. The entire collapse is carried by the DISCRETENESS of the zero
   set. Bent primes (Λ(2)×1.05) stay invisible in the sparse regime and go negative EXACTLY at the edge —
   "marginality/no-cushion" measured: arithmetic is only tested at critical density.
2. **Exact fluctuation identity:** zeroGram − (arch-density Gram) = pole − prime, EXACTLY (explicit formula).
   The fluctuation form IS the prime side; its norm is O(1) but it DECOUPLES from the minimizer direction in
   the mid regime (norm bounds certify nothing; directional control is everything).
3. **Literature status (§336, scanned + spot-checked):** foundations all published — Borichev–Gröchenig–
   Lyubarskii 2009 (Gaussian-Gabor frame bounds near critical density), Landau–Widom (prolate plunge),
   Bondarenko–Radchenko–Seip + Kulikov (zeros as interpolation set at exactly RvM density), Groskin
   arXiv:2607.02828 (truncated Weil machinery), Kurasov–Sarnak (zeros NOT a Fourier quasicrystal). The
   SYNTHESIS — "Weil positivity margin = frame bound collapsing at RvM critical density; rank = zero count" —
   appears unpublished. Do not re-derive the components.
4. **Certificate chains (§337-339):** rigorous zero-free lower-bound chains on λ_min(A−P), A=arch+pole,
   P=prime (both zero-free). 2-block gap bound + adaptive exact prime core reaches the WEYL RANK LIMIT on every
   band tested ([14,44],[50,80],[100,130]) and every width tested (w=1.5..4), tracking true margins at
   95-99.9%. The directional block NEVER grows (k=3 always suffices); only the exact core grows.
5. **Bragg/RS structure (§337-338):** at each probe density ONE prime power is near-resonant on the minimizer
   (log n* ≈ 2π/s; others suppressed 10-1000x by equidistribution); at critical density the resonant scale is
   n* = γ/2π — the Riemann–Siegel crossover. The prime form's diagonal cost on the minimizer is tiny
   throughout (|p11|≤0.09, negative/helping in the sparse regime); the obstruction is off-block coupling.
6. **Arb-rigorous (§340):** B_full PD at the plunge floor s=3.5 AND s=3.4 (band [14,44], w=3) by interval
   LDL^T — theorem-grade for these configs, margins 1.9e-10 / 9.4e-12.
7. **Fully a-priori chain (§341):** replacing the arch integral by an explicit piecewise-Riemann MINORANT
   (rigorous by Ω-monotonicity, erf-explicit) makes every chain ingredient an elementary closed formula (no
   integration, no zeros). Certifies to s=4.2 (83% of true margin from formulas alone); loss (~3e-4) localized
   ([0,7] flat floor + partition residue) with mechanical fixes to reach ~3.9.
8. **THE WALL, MEASURED (§342):** the chain is UNIVERSAL (reaches the rank limit at every w and band — no
   structural saturation anywhere). The uniformity wall appears as a COST LAW: the exact core needed at the
   edge grows ~e^{c/w} (max element 7→11→17→29→113 as w: 4→1.5) — as the test class completes (w→0), the core
   → ALL primes. "Prove RH on this front" = bound the core-growth cost as the class completes. Difficulty
   conservation as a measured growth law, not a mood.

## Open questions this table could actually move
- **Q1 (the core-growth attack):** can the e^{c/w} exact-core growth be beaten a-priori? The Bragg
  equidistribution suppression (10-1000x per off-resonant prime) is measured but uncharged in any bound —
  is there a resummation/ambiguity-function argument that controls the INFINITE core using the oscillation,
  turning the growing core into a convergent structure? This is the sharpest open target on the front.
- **Q2 (the frame-theory transplant):** Seip-type density theorems are UNIFORM over all configurations at
  once — the kind of uniformity matrix analysis can't reach. Can Gabor/Fock density theory prove a statement
  uniform over the test class where our per-configuration chains can't? What exactly would "the zeros form a
  maximal-density real sampling set for the Gaussian family" assert — RH-equivalent, weaker, or a genuine
  middle object? Run it through the 3-clause gate.
- **Q3 (finite ingredients):** the two named analysis pieces for full a-priori-ification at fixed w:
  localized-weight Toeplitz floor for the arch matrix; finite-n λ₄(A) lower bound (Bogoya–Böttcher–Grudsky).
- **Q4 (the RS endpoint):** the resonant scale at criticality is n*=γ/2π — the same crossover as
  Riemann–Siegel / Berry–Keating. Is there a duality argument (short-orbits ↔ long-orbits) that handles the
  critical resonance term the way RS handles the main sum — i.e., the last enemy at the edge has a known dual
  representation?
- **Q5 (honest independence check):** is the sampling identification a genuinely different coordinate, or
  will it, pushed to uniformity, re-encode as the same wall (per the edge-vs-center law)? What's the
  conservativity status of "zeros = maximal-density sampling set"?

## Standing rules for the table (unchanged)
LP/Jensen fixed-degree = closed wall. Explicit formula = banked trace identity (scores nothing). Criteria
certify, objects prove. 3-clause gate + constructible-without-zeros + manifest positivity. Check t01
route-prefilter + this doc before claiming novelty. Tools live in ~/projects/riemann/operators/ (l2_edge_sweep,
directional_cert, resonant_cert, band_sweep_cert, kr_frontier, arb_dense_cert, apriori_chain, arch_floor,
w_sweep_chain) and ~/projects/riemann/arsenal/.
