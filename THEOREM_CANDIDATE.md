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
