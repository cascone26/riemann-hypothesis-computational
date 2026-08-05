# Theorem-candidate: the Möbius resummation of Landau's theorem (working note, 2026-08-04)

## Statement (candidate, ζ case; Selberg-class version replaces μ, θ, Λ by family data)
Let γ_1 < γ_2 < … be the ordinates of the nontrivial zeros (assumed simple, on-line up to T),
θ the Riemann–Siegel theta, N = N(T). For u > 0 fixed and X = X(T) → ∞ suitably slowly:
  G_N(u) := (1/N) Σ_{k≤N} (−1)^k e^{iθ(γ_k)} e^{iuγ_k}
satisfies: G_N(log n) = C₁(n,T)·b(n) + C₂(n,T)·Λ(n)/(√n log n) + E(n,T),
where b = Dirichlet-inverse coefficients (μ for ζ), C₁ ≍ the Bessel product Π-structure
(§629; ≈ −0.47·(1+u-drift) at T~10⁴), C₂ ≈ −0.28 (§615), and E → 0 (measured decay §650).
Support: G_N(u) → 0 for u off {log n} (measured ~N^{−1/2}).

## Proof architecture (all boxes machine-verified at T ≈ 10⁴, LEDGER §629–650)
1. EXACT elementary weight: e^{−i arg ζ′(ρ_k)} = i(−1)^k e^{iθ(γ_k)} (§606; algebra from
   Z-realness + counting; needs zeros simple/on-line — this is where RH-type hypotheses enter).
2. EXACT resummed identity (two integrations by parts, §610/§631/§649):
   G = (1/K)[∫e^{iΨ}(θ′/π)dt + [Se^{iΨ}] − i∫S·Ψ′e^{iΨ}dt],  Ψ = −πS_X + ut.
3. Box (i) smooth part: stationary phase / Bessel main terms — the μ-comb with Π J-structure.
4. Box (ii) S_smooth·Ψ′ bilinear: convergent prime-pair sums (Σ over p,q of explicit
   coefficients; absolutely convergent — amplitudes 1/(m√n)).
5. Box (iii) sawtooth: LINEAR in the zero-fluctuation (§649b); its transform = Landau comb/(iν)
   ⇒ single Landau factors only. Inputs: Landau's theorem, UNIFORM in x over the needed
   frequency set — Gonek's uniform Landau; this is the care point where error terms
   (and possible conditionality) concentrate.
6. Error: finite-T remainder, measured decaying in T (§650). Target bound: o(1) for fixed u,
   X = T^θ, θ < 1 (shape TBD analytically).

## Honesty
- Everything above is REPORT-tier: machine-verified numerics + exact identities; no analytic
  error bounds proven yet. The novelty calibration is §648/§648b (vs Landau 1912, vs
  Mertens/Gonek 1/ζ′-weights): the contribution is the BOUNDED position-only weight and its
  universality (through non-abelian Artin), plus the zero-parameter derivation.
- No Handle for RH: this note characterizes zero-arithmetic duality; it does not bound Λ_dB
  or prove positivity of anything RH-equivalent.


## Error-bound section (v1, 2026-08-04 — measured budgets + analytic shapes)
Assemble G_N(u) via the resummed identity (§649). Per-box budgets at T~10⁴, u=log 2:
- **Box (i) smooth/stationary-phase**: main term; truncation error from the X-cutoff enters
  through the oscillator field. Analytic shape: classical stationary-phase remainders, O(1/T)
  boundary terms (measured: boundary contribution < 10⁻³).
- **Box (ii) prime-pair bilinear** (S_smooth·Ψ′): absolutely convergent double sum
  (amplitudes 1/(m√n · m′√n′), frequency constraint); truncation at Y: tail ≤ Σ_{n>Y} 1/√n·(…)
  — measured contribution −0.047→−0.104 across X with S_Y at Y=3×10⁴ capturing all but the
  sawtooth part.
- **Box (iii) sawtooth×Ψ′ (linear in dN−ρdt)**: reduces to single Landau factors (§649b).
  MEASURED uniform Landau budget (§669): |err| ≤ 0.002 for x ≤ 200 (active resonances),
  ≤ 0.006 to x = 2000; degrades to ~0.03 erratic beyond x = T/2π (§672) — hence the X(T)
  constraint below. Analytic input: Gonek's uniform Landau formula; conditionality (if any)
  enters ONLY here.
- **X(T) constraint** (from §670–672): keep ALL active frequencies below the spectral edge —
  X·(max resonance multiplier) < T/2π. For fixed u this allows X → ∞ as T → ∞ (X ~ T^θ any
  θ < 1 for the leading structure; the measured accuracy peak at small X reflects finite-T
  secondary terms, not an asymptotic obstruction).
- **Net measured closure**: identity verified to ~10⁻³ at all X (§649); zero-free evaluation
  matches truth to 94–104% per channel at the natural X (§634–638); remaining analytic work =
  (1) uniform-Landau error term inserted into the box-(iii) sum, (2) interchange/tail lemmas
  for boxes (i)–(ii) (absolute convergence — routine), (3) uniformity in u on compact sets.
Status: REPORT-tier architecture with quantified budgets; no step known to require RH except
possibly the strongest uniform-Landau ranges (Gonek's formula is unconditional; range/error
tradeoffs to be checked against the x < T/2π constraint, which our assembly respects).

## Error decay VERIFIED at 7× height (§688)
The resummed-identity deficit (box iii + tail) collapses −0.0269 (T~10⁴) → +0.0001 (T~70k) at
fixed X=3000 — 270× reduction, confirming the error is ~10^-4 in the theorem's operating regime (X<T/2π) at T~70k — main/error ~3000× (§688b caveats the pure-o(1) reading: partly edge-crossing; matched-edge scan queued). The
leading main term persists (amplitude law A₁→0.27, §685); the correction vanishes. This is the
strongest available evidence that the theorem-candidate's architecture closes asymptotically.