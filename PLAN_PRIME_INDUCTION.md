# PLAN: Prime-by-Prime Induction on Weil Positivity (drafted 2026-07-20, riemann4, on Fable; execute on Opus)

## The idea (ours — from §52/§53/§56 mechanism work)
Each prime power n enters the localized Weil form τ(λ²) at threshold **λ² = n with exactly zero weight**
(window factor 1 − log n/L, L = log λ², vanishes at n = λ²) and turns on **continuously**. §53 verified the
passage through λ²=3 is smooth; §52 showed the entering prime *takes over stabilizing* the directions the
archimedean term is losing. So decompose RH:

- **Base case (PROVEN):** prime-free regime λ² < 2 — Suzuki Thm 1.4, λ_a > 0 unconditionally.
- **Inductive step (THE TARGET):** if τ(λ²) ⪰ 0 just before prime power n enters, positivity survives
  through n's entry window [n, n′) (n′ = next prime power).

One infinite global statement → countably many local-in-prime steps, each: one prime, one window, one
continuous deformation, with the §52 exact-cancellation identities in hand.

**Honest framing upfront (§54):** there is NO margin buffer — the step must preserve an exact knife-edge
(min-eig pinned ≈ 0 by the Ξ-section null tower), not a cushion. The difficulty concentrates into the step;
the bet is that the step is *structured* enough (rank-additive prime entry, zero-weight onset) to crack.

## The sharp candidate invariant (formulate precisely in Phase 1)
Observed: min-eig(τ(λ²)) ≈ 0 for ALL λ² (§53 knife-edge pinned). If exactly true, then along the flow the
first-order condition is **⟨v₀ | dτ/dλ² | v₀⟩ = 0** where τ v₀ = 0 (null vector ⊥ the flow derivative).
- Prime entries preserve this automatically at onset (entering term has zero weight AND its coefficient's
  λ²-derivative acts only through the window — compute it).
- ⚠ TRIVIALITY TRAP to check: if min-eig ≡ 0 identically, its derivative is trivially 0 — the *content* is
  that min-eig stays 0 rather than dipping negative, i.e. the SECOND-order form on the null space + the rest
  of the spectrum staying ≥ 0. The invariant that matters may be the second-order object. Phase 1 must
  disentangle this before Phase 2 measures it.

## Phases

### Phase 0 — Rule-13 prior-art check (subagent, ~15-20 min) — GATE
Has anyone run induction-on-primes / add-primes-one-at-a-time for Weil positivity?
Search: "Weil positivity finitely many primes", "semilocal trace formula positivity", "Weil criterion finite
Euler product induction", Beurling add-a-prime monotonicity, "smallest eigenvalue monotone in primes".
**Known prior art to position against (MUST read): Connes–Consani arXiv:2006.13771 "Weil positivity and
trace formula" — proved SEMILOCAL positivity (archimedean + finitely many places).** Critical question: does
their result already give every finite stage (making "induction" vacuous and the limit the whole content —
the known story), or is our λ²-coupled family (prime set grows WITH bandwidth L = log λ², Suzuki A_a) a
genuinely different filtration where the per-step statement is NOT theirs? Our filtration ties test-function
support to the prime cutoff — that coupling is the novelty candidate. If the framing is dead in the
literature → fallback (see Gates).

### Phase 1 — Precise formulation (~30-45 min, analytic writing)
1. Fix the family: τ(λ²) = A + B − Σ_{n ≤ λ²} w_n(λ²) Q_n, w_n(λ²) = Λ(n) n^{-1/2} (1 − log n / log λ²)
   [verify exact form from operators/suzuki_screw_lambda_a.py — code is truth].
2. State H(n): τ(λ²) ⪰ 0 for all λ² ∈ [2, n]. Step: H(n) ⟹ H(n′).
3. Write the step as a perturbation lemma: on [n, n′], τ = τ_smooth(λ²) − w_n(λ²) Q_n with w_n(n) = 0.
   What controls a zero-crossing of the min eigenvalue: first-order (the derivative identity above), and if
   that vanishes, the second-order quadratic form on the null space N(λ²). Write BOTH explicitly.
4. Disentangle the triviality trap: state exactly which object must be shown ≥ 0 per step.
5. Deliverable: one page, the Step Lemma + its explicit-formula meaning.

### Phase 2 — Numerical probes (Olympus, dps ≥ 30, ~1-2h)
- **P2a:** min-eig(τ) and 2nd even eigenvalue on a FINE λ² grid through thresholds 3, 4, 5, 7, 8, 9, 11.
  Is the knife-edge pinned exactly (≡0 within noise) or structured dips/rises between thresholds?
- **P2b:** the handoff law, quantified: at each threshold compute ⟨v₀|dτ/dλ²|v₀⟩ just before/after entry;
  decompose dτ/dλ² into (arch-drift) + (window-drift of established primes) + (entering prime). Verify the
  entering term contributes 0 at onset and measure the first-order compensation (§52 quantified).
- **P2c:** null-vector continuity: does v₀(λ²) rotate smoothly through a threshold? Do NEW null directions
  appear at thresholds (new Ξ-sections, §51 cascade) or between them?
- Tools: extend operators/session51_59/{parity_blocks,threshold_expand}.py; archive results + scripts in
  operators/prime_induction/ (new dir); log runs in LEDGER.

### Phase 3 — The analytic step (the real math; first pass ~1-2h, then open)
- Prove what's provable in the Step Lemma: smoothness between thresholds; zero-weight onset; the rank/sign
  structure of Q_n on the descending directions (§52 says ⟨u|Q_n|u⟩ has the compensating sign — prove it for
  the explicit q-kernel).
- Isolate the HARD KERNEL: the piece that resists. Two acceptable outcomes:
  (a) proof sketch of the step with a clearly-bounded hard kernel (progress), or
  (b) proof that the step's kernel = global positivity verbatim (a §54-style door-closure — still valuable,
      record and pivot).
- Watch for: the second-order form on the null space being COMPUTABLE (explicit q-kernels + explicit null
  vectors = Ξ-sections) — if it is, the step may reduce to a verifiable inequality per prime. That would be
  the breakthrough shape.

### Phase 4 — Position vs Connes–Consani semilocal (~30 min)
Read 2006.13771 properly (membrane-pull the source). Write the comparison: their fixed-adelic-space finite-
place positivity vs our λ²-coupled filtration. Adjust claims; cite in LEDGER. Do NOT re-prove their result.

## Gates / fallbacks
- **Gate 0** (after Phase 0): framing already tried & dead → fallback to engaging Connes' prolate/Poisson
  convergence formulation directly (2602.04022), our §48 PW_λ in hand.
- **Gate 2** (after Phase 2): knife-edge NOT structurally maintained (real dips at high precision, not
  discretization) → the induction hypothesis is wrong as stated; reassess with the measured profile.
- **Gate 3** (after Phase 3): hard kernel = global positivity verbatim → record door-closure §-entry, pivot
  to geometric individual-isolation (§59) as the standing deep thread.

## Discipline
- Rule 13 before every new sub-attempt; verify-before-ship on every numeric; LEDGER § per result;
  PROCESS.md append per block; archive ALL scripts in repo (no /tmp-only tools); Olympus for compute;
  model: Opus for Phases 1/3 (real math), Sonnet fine for Phase 2 runs & Phase 0 subagent.
- ETAs: Phase 0-2 bounded (~2-3h total). Phase 3 open-ended — report the kernel honestly when isolated.
