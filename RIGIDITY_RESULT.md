# Arithmetic Rigidity of the Zero Configuration (canonical note, §651–661, 2026-08-04)

## The result
The nontrivial zeros are, in an instrumented and quantified sense, **the unique configuration
consistent with the primes**. Solving the multiscale self-consistency system
  θ(t_k)/π + 1 + S_δ(t_k) = k − ½,   S_δ = Gaussian-smoothed prime series,  δ: 0.8 → 0.12
by damped Newton (full Jacobian θ′/π + S′_δ, floored, step-clipped) recovers the true zeros:

- **The law**: rms error = (1.04 ± 0.03)/(log X)² mean spacings, verified over three decades
  (X = 500 → 5·10⁵; 0.6% of a spacing at X = 5·10⁵; 5000 zeros simultaneously). [§652–653]
- **Local rigidity**: any single zero displaced ½ spacing is detected and repaired to floor;
  neighbors undisturbed. [§654]
- **No alternative solutions**: 15/15 defect configurations (pair-merge = the off-line
  double-ordinate signature, triple-merge, block-shift × 5 heights) repaired to floor by the
  multiscale schedule. Fine-scale-only solves DO stick (narrow basins, §656) — multiscale is
  the right functional. [§655–658]
- **Effectively global basin**: full recovery from iid corruption of ALL zeros at 4, 8, 16
  mean spacings (99.7% at floor from a start with no local structure); frays only at ~32 sp
  via isolated stragglers — no competing attractor. Edge scales with the anneal schedule
  (instrument-limited). [§659–660]
- **Universal**: Δ's 683 degree-2 zeros recover from Ramanujan λ(p) data the same way
  (0.054 sp at X=300; ζ-law shape holds, constant ~1.7×). [§661]

## Anatomy
- Optimal smoothing δ* ~ 1/log X (contraction limit — S′_δ must stay below θ′; derived §651,
  matches measured optimum). Full-Newton is second-order in δ* ⇒ the (log X)^{−2} law.
- The constant ≈1 is COLLECTIVE: the solve beats the pointwise S-residual bound by 4.4× via
  correlations between neighboring counting constraints (covariance-driven; the analytic
  derivation goes through the sawtooth-field covariance — same object as the mechanism note's
  box (iii)). [§657]

## Relation to RH (honest scope)
Within counting-function data, an off-line pair at height t₀ is indistinguishable from a
double ordinate — exactly the pair-merge defect, which the prime field rejects and repairs
(residual 3.4× detection scale). This is rejection by the actual primes at accessible heights
and resolution; the theorem-class statement (uniqueness of the fixed point over admissible
configurations, uniformly in T) is the standing analytic target. No Λ bound; No Handle for RH.

## Scripts
operators/inverse_v3.py, inverse_v3_xscan.py, rigidity_probe.py, exclusion_map.py,
exclusion_anneal.py, basin_radius.py; Δ: inline §661. Companion notes:
MOBIUS_PHASE_RESULT.md (forward direction), MECHANISM_DERIVATION.md, THEOREM_CANDIDATE.md.
