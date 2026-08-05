# RH Apparatus — Build Plan (2026-08-06)

Goal: a self-sustaining apparatus that converts RH-attack effort into **compounding verified assets**
(Handles), not decaying claims (Reports). Five organs; the verifier is the multiplier. Strategy
per Jacob's plan + my strategic amendments (SDP parallel; organ-4 as critical-path risk; point at
live frontier; certified obstructions are first-class).

## Environment (oriented 2026-08-06, machine = Sam / MacBook, 8 cores / 8GB / 15GB free)
- **Lean 4** toolchains v4.18.0 + v4.29.0 (elan/lake). **RiemannLean project EXISTS** with mathlib
  dep. Already contains a ZERO-SORRY Handle (`TailLemma.tail_kernel`) + a partial Robin lane
  (4 sorries). mathlib cache NOT present → `lake exe cache get` needed (~5GB, fits in 15GB).
- **SDP/SOS**: cvxpy 1.9.2 + scs + sympy/mpmath ready NOW (organ 2 unblocked).
- **Ollama** installed, NO models pulled; 8GB RAM caps local models to small quants. Heavy Lean
  builds + models → prefer **Olympus** (GPU box, `ssh olympus`) per fleet; Sam prototypes.
- Repo has GitHub remote (cascone26/riemann-hypothesis-computational) → fleet-dispatchable.

## Organ 1 — Lean verifier + propose→verify loop  [THE MULTIPLIER — build first]
Not from scratch: RESTORE + WRAP.
- 1a. Restore mathlib cache in RiemannLean; re-verify `tail_kernel` still checks (0 sorry) → the
  kernel is live. [running]
- 1b. Minimal propose→verify loop: a script that takes a candidate lemma (statement + tactic
  proof), writes a Lean file, runs `lake build`, and returns PASS/FAIL + errors. PASS = Handle.
- 1c. Handle ledger: verified lemmas recorded with their Lean file hash → immune to my error rate
  + context amnesia (a verified lemma stays verified across sessions). This is the Charisma
  Handle-minter made rigorous.
- First real target: close the 4 Robin sorries (RH-equivalent lane, already scaffolded).

## Organ 2 — SDP / Weil-positivity certificate lane  [PARALLEL FROM DAY 1 — self-certifying]
Amendment to Jacob's order: does NOT wait behind organ 1. An SOS certificate (rational PSD Gram
matrix) is self-verifying by exact arithmetic — a Handle WITHOUT Lean.
- 2a. Set up the Weil quadratic form as an SDP: W(φ)=arch+pole−prime; search for an SOS
  decomposition on a test-function basis (finite level of the hierarchy). cvxpy+scs numerically,
  then rationalize + exact-verify the Gram matrix (sympy) → certified.
- 2b. Every feasible level = evidence + shape-constraint on what a certificate must look like.
  Saturation regime (our "barely true") is where the hierarchy is most informative — the
  feasibility boundary sits at the margin.
- 2c. Pincer with the Krein/operator lane: certificate shape constrains the operator; operator
  structure constrains the certificate.

## Organ 3 — Loop economics / tiered router  [after 1; scale]
- Local models (Ollama, small quants on Sam; bigger on Olympus) for grunt work (search,
  formalization drafts, tactic attempts) — the verifier catches their errors for free, so weak
  models become usable. Frontier calls reserved for genuinely hard steps.
- Router: cheap-first, escalate on FAIL. Makes 24/7 affordable.

## Organ 4 — Autoformalization grind  [after 1+3; CRITICAL-PATH RISK, not background]
Amendment: this is a GATE, not free throughput. mathlib's analytic-NT toolkit is THIN at the
frontier (explicit formula, zero-density, Selberg/Gonek mostly absent). First target should be
FORMALIZING A LIVE RECENT RESULT (e.g. the dBN PF_5 paper — small, self-contained, checkable) to
earn Handles on current mathematics fast, THEN the classical toolkit. Value certified OBSTRUCTIONS
as first-class output.

## Organ 5 — Geometry / F₁ candidate-tester  [last; the shaped moonshot]
By this point: certificate shapes (2) + operator structure (Krein) + verified constraint library
(4) narrow the search → "test candidates against a hard constraint list," which the stack is built
to run. Explicit bet: organs 2 & 5 are the SHOTS ON GOAL (search finding structure); 1/3/4 are the
machine that makes rolling them cheap.

## Krein lane — does not pause (already producing; feeds 1 and 2)

## Build order: restore-verifier (1a) ∥ SDP (2) → loop (1b/1c) → router (3) → formalize-live (4) → geometry (5)
Honest risks: (i) mathlib library gap gates organ 4; (ii) RH's binding constraint may be
idea-EXISTENCE not candidate-quality — the apparatus buys cheap shots-on-goal, which is the bet.
