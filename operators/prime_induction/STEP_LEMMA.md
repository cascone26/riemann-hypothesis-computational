# Phase 1 — The Step Lemma (prime-by-prime induction on Weil positivity)
Drafted 2026-07-20 (Opus), riemann4 §76 execution. Continuous formulation from Suzuki's screw function
(operators/suzuki_screw_lambda_a.py = code-truth).

## Setup (continuous, Suzuki eq 1.3)
Fixed even function
  g(t) = −4(e^{t/2}+e^{−t/2}−2)                              [POLE — negative, exponential]
       − (|t|/2)·(ψ(¼)−log π)  = +2.686·|t|                  [ARCH CORE — positive, linear (Ω(0)=−5.372)]
       − ¼[Φ(1,2,¼) − e^{−|t|/2}Φ(e^{−2|t|},2,¼)]            [ARCH TAIL — bounded]
       + Σ_{n≤e^{|t|}} Λ(n)/√n·(|t|−log n)                    [PRIMES — each a positive ramp]
λ_a = inf_{v∈C_c^∞(−a,a)} Q(v)/‖v‖²,  Q(v)=∫∫ g(x−y) v'(x) v'(y) dx dy.  RH ⟺ λ_a ≥ 0 ∀a>0.

Crucially g is FIXED; the ONLY a-dependence is the support (−a,a). Since v'(x)v'(y) is probed at |x−y|≤2a,
prime n is INVISIBLE for 2a<log n and becomes active at **a_n = (log n)/2** (⟺ λ²=e^{2a}=n). At a=a_n the
ramp (|t|−log n)_+ meets the support only at the single point |t|=2a=log n, where the ramp = 0 → **prime n
enters with zero weight, continuously** (verified: q-kernel window (1−log n/L) and off-diag sin(2πm)→0 at
λ²=n, LEDGER cascade builder).

## Screw (Krein) decomposition — the sign of each piece
For even h with h(0)=0, ∫∫h(x−y)f(x)f(y) with f=v' mean-zero (∫v'=0 automatically): the sign is governed by
h''  (distributional). h is a POSITIVE screw increment iff h'' ≥ 0 (as a measure) [Krein–Langer].
- **Arch core** +2.686|t|:  (|t|)'' = 2δ(t) ≥ 0  → POSITIVE increment (weight 5.372·δ). Dominant near 0.
- **Prime n ramp** Λ(n)/√n·(|t|−log n)_+ :  ''  = Λ(n)/√n·(δ_{log n}+δ_{−log n}) ≥ 0 → POSITIVE increment,
  a point mass at ±log n of weight Λ(n)/√n. **Each prime ADDS positivity, localized at |t|=log n.**
- **Pole** −4(e^{t/2}+e^{−t/2}−2):  '' = −(e^{t/2}+e^{−t/2}) < 0 → NEGATIVE increment, exponential. THE obstruction.
- Arch tail: bounded, sub-dominant.

So **g'' = 5.372·δ₀ + Σ_n Λ(n)/√n(δ_{log n}+δ_{−log n}) − (e^{t/2}+e^{−t/2})** (+ tail). The positive part is a
DISCRETE atomic measure at {0}∪{±log n} (arch + primes); the negative part is the smooth exponential pole.
RH ⟺ the atomic positives dominate the exponential pole in the screw/Rayleigh sense on every (−a,a).

## The induction
- **Base case (PROVEN):** a < a_2 = (log2)/2 = 0.347. No primes. λ_a > 0 — Suzuki Thm 1.4
  (λ_a = log(1/a)+μ₁−log2π+ψ(2)−1+O(a), μ₁≈0.25 purely archimedean; code verifies primes empty here). ✓
- **Inductive step H(n)⟹H(n′):** if λ_a ≥ 0 for all a ≤ a_n, then λ_a ≥ 0 for a ∈ (a_n, a_{n′}]
  (n′ = next prime power).

## Why the step is NOT trivial (triviality trap, resolved)
Naive: "each prime adds a positive increment (h''≥0), so Q only increases, so positivity is preserved." FALSE —
because the DOMAIN grows simultaneously. Two competing effects as a increases through a_n:
  (i) domain growth: v' can now reach |x−y| up to 2a, accessing MORE of the negative pole −4e^{|t|/2}
      (whose magnitude ~4e^{a} at the new boundary) → pushes λ_a DOWN.
  (ii) prime-n activation: the new positive point mass Λ(n)/√n·δ_{±log n} at |t|=log n < 2a → pushes λ_a UP.
The step lemma is the statement that (ii) ≥ (i) at every threshold — a RACE, not a monotone gift.

## The hard kernel (honest — where RH concentrates)
Quantitatively: at the moving boundary |t|=2a the pole magnitude is 4e^{a}=4√(e^{2a})=4√(λ²)=4√n. The total
prime mass activated by then is Σ_{m≤n}Λ(m)/√m ~ 2√n (PNT). **Both scale as √n** — the prime mass and the pole
magnitude grow at the SAME rate (√n), and RH ⟺ they balance to keep λ_a≥0. This is the √x-strength PNT
cancellation (§66) re-localized: the step lemma at threshold n requires the increment of prime mass Λ(n)/√n
positioned at |t|=log n to offset the increment of pole-access from domain growth dt=da. The HARD KERNEL is:
  **does the positive point mass Λ(n)/√n at |t|=log n, weighted by |v_a'(x)v_a'(y)| on the diagonal band
  |x−y|=log n, dominate the pole's contribution on the newly-accessible corner |x−y|∈(2a−da, 2a)?**
This is a LOCAL (per-prime), EXPLICIT inequality on the minimizer v_a — a genuinely more structured object than
"global λ_a≥0". Acceptable Phase-3 outcomes: (a) prove it with a bounded residual, or (b) show it = global
positivity verbatim (door-closure).

## Phase-2 measurables (to characterize the race before Phase 3)
Using suzuki_screw_lambda_a.py (λ_a on the CONTINUOUS domain — the right object, unlike §53's fixed-loc τ):
1. λ_a on a FINE a-grid through a_2,a_3,a_4,a_5,a_7 (=½log{2,3,4,5,7}=0.347,0.549,0.693,0.805,0.973). Does λ_a
   DIP between thresholds and get RESCUED at each a_n (the relay), or stay pinned? (§53 said pinned for τ, but
   that was fixed-loc — RE-TEST on λ_a.)
2. dλ_a/da near each a_n: the (i)−(ii) race — sign just below vs just above a_n.
3. The minimizer v_a: does it concentrate its v' mass at the diagonal band |x−y|=log n of the entering prime?
   (would confirm the local-inequality picture).
NOTE Phase-2 caveat: suzuki_screw_lambda_a.py flags λ_a<0 at a≳1 as a DISCRETIZATION artifact (pole exp-growth
leak). So trust the a-grid only up to a≈0.9 (through prime 7). Beyond needs better discretization.

## Prior-art gate (Phase 0, running): MUST position vs Connes–Consani 2006.13771 semilocal positivity.
If their semilocal result already gives λ_a≥0 for every finite prime set unconditionally, the "induction" base+
steps are all their theorem and only the a→∞ limit is content (known story) → fallback. Novelty candidate: our
BANDWIDTH-COUPLED filtration (support 2a tied to prime cutoff e^{2a}), where prime n enters exactly as the
domain reaches it — is that their filtration or a genuinely different one?
