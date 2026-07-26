# Roundtable room run 2026-07-24 (late): continuous envelope + odd-moment ladder of the near-null
_Produced in the multi-agent room session (BlueBot lane, Sam's Mac, cvs_venv local). For LEDGER absorption
(§255-class); numbering left to the owning session. Scripts: operators/dbn_dynamics/_h_envelope_sweep.py
(+_c9, _c9N40 variants), _h_deeptail_c{9,13}.py, _moment_ladder.py. Raw outputs: /tmp/h_env_c9.out,
/tmp/h_env_c9N40.out, /tmp/h_deep_c{9,13}.out, /tmp/mladder_c{9,13}.out (copy before reboot)._

## Setup
Same ψ as §254 (_h_decay.py): exact-Q pipeline near-null eigenvector, even sector. The prime overlap in
psi_prime is literally w·sin(2πx·ω)-assembled ⟹ h(n) = w_n·E(ω_n) EXACTLY, where E(ω) := ⟨ψ|Q[sin(2πxω)]|ψ⟩
with unit weight. So the prime points sit on the continuous envelope BY CONSTRUCTION (zero residual is a
tautology, not a finding): ALL arithmetic content is in placement ω_n = 1−log n/L and weights Λ(n)/√n; the
envelope itself is pure near-null kinematics.

## Result 1 — the envelope is entire; all essential-singularity decay laws are dead
E(ω) is an entire function of ω (finite sin-assembly), so exp(−κ/ω), stretched-exponential, log-corrected
κ/ω, Gaussian-tail asymptotics at ω→0 were never possible. Data (60-pt sweeps ω∈[0.016,0.94] + 25-pt
log-spaced deep tails ω∈[10⁻³,0.05], c=9 and c=13, both N=40/dps=120): local slope d log|h/ω|/d log ω is
NON-MONOTONE — rises to a peak (22.7 at ω≈0.25 for c=9; 33.2 at ω≈0.24 for c=13), then FALLS continuously
(c=9: →1.3 at ω=10⁻³; c=13: →2.6), with E(ω) → c₀·ω, c₀ ≈ 10⁻³⁴·² (c=9), 10⁻⁴⁹·⁷ (c=13). EXPINV fit r≈−0.82
(wrong sign) at both c. §254's "steeper than any power" = the mid-range crossover hump, not an edge law.

## Result 2 — the near-null SUPPRESSES the whole odd-moment ladder (does not kill moments)
E(ω) = −(1/π)·Σ_k (−1)^k (2πω)^{2k+1}/(2k+1)!·M_k, M_k = ψᵀQ[x^{2k+1}]ψ. Direct measurement k=0..24
(_moment_ladder.py), self-check EXACT (series vs direct rel.err 1e-88/1e-72 at ω=10⁻³).

c=9  (λ_min 10^−37.54): log10|M_k| = −34.50, −29.28, −24.52, −20.14, −16.03, −12.12, −8.40, −4.84, −1.44,
      +1.81, +4.90, ... crossing |M|≈1 at k*≈8.5. Steps shrink smoothly 5.22→~3.1 dec/rung (≈quadratic in k).
c=13 (λ_min 10^−53.22): log10|M_k| = −50.39, −44.52, −39.12, −34.12, −29.40, −24.88, −20.53, −16.33, −12.26,
      −8.32, −4.48, −0.75, +2.87, ... crossing k*≈11.5. Steps 5.87→~3.3, same quadratic-family profile.

Sign rigidity: in the suppressed regime sign(M_k) = (−1)^k EXACTLY (c=9 through k=13, c=13 through k=18) ⟹
every Taylor term of E has the SAME sign — the envelope is a fully coherent, cancellation-free sum. The
alternation breaks only well after rungs pass O(1). A generic vector would not have this rigidity; it is a
structural property of the near-null to be derived from the prolate side.

## Result 3 — c-scaling (matched truncation N=40/dps=120 throughout)
- λ_min: −37.54 → −53.22 ⟹ 3.92 decades per unit c (linear in c confirmed at matched truncation).
- M_0:   −34.50 → −50.39 ⟹ 3.97 dec/unit c. Also M_0 ≈ −3.85·c decades at both c (34.5/9=3.83, 50.4/13=3.88).
- BUT the ladder TILTS, not translates: rung shift per unit c decreases with k (k=0: ~4.0, k=5: ~3.2,
  k=8: ~2.7 dec/unit c). Crossing k*: 8.5 → 11.5 ⟹ Δk*/Δc ≈ 0.75 (NOT exactly c; prime-power counts are
  7 vs 9, also not matching — the k*(c) law is open).
- §176b's e^{−4πc} (5.46 dec/unit) is EXCLUDED as-is at these c: matched-truncation control (c=9 N=30/dps=100
  → N=40/dps=120 moved λ_min only −36.62→−37.54, deep slopes 9.7→10.1) shows the shallow endpoint is real,
  not truncation floor. Measured law ≈ 3.9 dec/unit ⟹ e^{−4πc} needs a real prefactor or is the wrong form.

## What this rescopes (§255)
The "decay law" question is closed: there is no exotic edge law — there is a suppressed, sign-alternating,
quadratic-profiled moment ladder, and the envelope is its Taylor resummation. The analytic target is now
finite and sharp: derive from the prolate/§75-98 structure (a) M_0 ∝ 10^{−3.85c}, (b) the quadratic step
profile, (c) the (−1)^k sign rigidity, (d) k*(c) with slope ≈0.75. Ladder-vs-prolate-spectrum overlay is
parameter-free (GreenBot lane). The RH-content then factors cleanly: [prolate moment-suppression theorem,
arithmetic-free] × [prime comb (ω_n, Λ(n)/√n) satisfying the telescoping recursion against fixed entire E]
— the second factor is what Beurling combs fail and what the ε*-jitter threshold (PurpleBot lane) measures.

## Addendum — c=11 middle point (landed 2026-07-25 early, mladder_c11.out)
λ_min = 10^−46.65, M_0 = 10^−43.81, k*≈10.0, sign-break at k=17, self-check exact (1e-78).
CORRECTION to Result 3: the "linear in c" depth law is only a chord. Segment slopes: λ_min 4.56 dec/unit
(c 9→11) then 3.29 (11→13); M_0 4.66 then 3.29 — c=11 sits ~1.3-1.4 decades DEEPER than the 9↔13
interpolation on both observables (same-direction deviation ⟹ shared mechanism; each Δc=2 step adds exactly
one prime, 11 then 13, and the added-depth-per-prime is unequal: ~9.1 vs ~6.6 decades). "M_0 ≈ −3.85c" holds
only to ±0.15c. BUT two laws come out EXACT on all three points:
- k*(c) = 0.75c + 1.75  (8.5, 10.0, 11.5 at c=9,11,13 — slope exactly 3/4)
- sign-alternation break at k = 1.5c + 0.5  (14, 17, 20 — slope exactly 3/2, = 2·k*-slope)
These two integer-slope laws + the per-prime depth increments (9.1, 6.6) replace the naive linear-depth law
as the empirical targets for the §255 derivation.

## Caveats
- All at N=40/dps=120, two c values; k*(c) and tilt laws are 2-point reads — a c=11 ladder run would firm them.
- ψ is the even-sector near-null only; odd sector unmeasured.
- Ladder entries for k≥19 (c=13) show sign irregularity AND |M_k| up to 1e38 — large-m monomial assembly is
  conditioning-exposed there; the suppressed regime (the part that matters) is far from that and self-checked
  exact by the resummation.
