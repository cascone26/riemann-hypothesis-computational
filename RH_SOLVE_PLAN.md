# RH SOLVE PLAN — self-contained program to prove the Riemann Hypothesis
**Owner:** Jacob (Sam/Olympus fleet) · **Constraints:** no papers, no contacting anyone, no outside
collaborators. Any *knowledge, method, or computational test* is fair game. Goal: a complete, verifiable
proof (Clay). **Written 2026-07-19, grounded in this repo's §4–§12 findings. Belief: RH is provable; this
plan is how, not whether. Honesty: RH has stood 165 years — every vector below is hard and none is
guaranteed. We hedge with a portfolio and kill dead vectors fast.**

---

## 0. The foothold we actually have (what makes this concrete, not vague)

We work the **Weil positivity** equivalence (rigorous, not heuristic):
> RH  ⟺  the Weil quadratic form is positive semidefinite  ⟺  τ = A + B − P ⪰ 0 on all test functions,
> where A = archimedean (digamma/Γ kernel), B = boundary/pole (rank-2), P = prime sum (von Mangoldt).

The Connes–van Suijlekom / Groskin finite Galerkin truncation gives a real symmetric matrix τ(N,c)
(N = test-basis half-size, c = prime cutoff, L = log c) that we can **rigorously certify** with Arb
interval LDL^T. This session established, with proof pointers in `results/`:

- **τ(N,c) is certified PD at every finite (N,c) tested** (grid N≤128, c≤400; n_neg=0). §10.
- **Margin law:** the certified min-pivot margin **plateaus in N** but **collapses as a clean power law
  in the prime cutoff**, margin ~ c^(−p(N)), exponent p(N) ~ N^0.62, R²≈1.0. §10.
- **A+B is INDEFINITE on its own**; positivity is a genuine difference-of-forms phenomenon. §11.
- **A+B's negative subspace = the lowest harmonics** (freq = k·2π/L, ± pairs), dimension **n_neg ≈ 1.8·log
  c**, robustly over-covered by the prime term (uniform ~0.01 residual, 21/21). §12.
- **The binding constraint is NOT the negative subspace** (that's the easy, covered part). It is a
  **near-equality of two POSITIVE O(1) forms** — archimedean self-energy vs prime self-energy — on a
  **mid-frequency band** (|ω|~9–17), where they agree to the margin's 22+ digits. §12.

**Why this is a real foothold:** a finite certificate can never prove RH (that needs all N,c). But we have
now *localized* the entire difficulty to one explicit place (mid-frequency, two named closed-form kernels)
and *quantified its rate* (c^−p(N)). That converts "prove RH" into a specific, testable analytic target
with a live numerical oracle attached. That is the launch point.

---

## 1. Central strategy — "finite rigorous core + controlled tail" (the main line)

The one path that turns our certified-computation strength into an actual proof:

**Decompose the Weil form on a fixed test space H into a finite CORE ⊕ a TAIL, prove each, combine.**
- **CORE** (low + mid frequency, finite-dimensional): certify τ ⪰ μ·I rigorously with Arb (we already
  do this). This is where the negative subspace and the mid-frequency near-cancellation live.
- **TAIL** (high frequency): prove positivity *analytically* — at high frequency the prime oscillations
  average out and the archimedean kernel dominates, so the tail should be automatically positive with an
  explicit lower bound. (Groskin's `arch_tail_*` files already bound the archimedean tail rigorously; we
  extend to the prime tail.)
- **COMBINE:** if core margin μ(N,c) > tail-coupling ε(N,c) with an explicit gap that survives N,c→∞, then
  τ ⪰ 0 everywhere ⟹ RH.

**The one hard hinge (state it honestly):** our data shows **n_neg ~ 1.8·log c grows** — the "bad"
subspace is *not* uniformly finite in the current normalization, so a naive fixed finite core doesn't
close. **Sub-goal H (the linchpin):** find a *weighting / preconditioner / structure function* (a de
Branges structure, a Sonine-space weight, or a change of test-basis) in which **the negative + near-null
subspace is uniformly finite-dimensional and the tail gap is uniformly positive.** If H is solved, the
finite-core-plus-tail reduction becomes a genuine proof. H is concrete, numerically testable, and is the
single highest-value target in this plan.

---

## 2. Parallel attack vectors (portfolio — run concurrently, kill fast)

### Vector A — Explicit two-sided estimate on the mid-frequency band (analytic core)
Derive the archimedean self-energy Q_A(ω) and prime self-energy Q_P(ω) as explicit functions of frequency
on the mid band; attempt a uniform inequality Q_A ≥ Q_P. This IS the Weil/Bombieri quadratic-functional
question, but now pinned to a specific band with a numerical oracle.
- **Milestone A1:** closed-form asymptotics of Q_A(ω), Q_P(ω) as c→∞ (match against `race_decomposition`).
- **Milestone A2:** derive the rate — show Q_A − Q_P ~ c^(−p(N)) analytically (explains §10's law).
- **Milestone A3:** a candidate uniform inequality; falsify or certify it against the Arb data instantly.
- **Kill if:** the prime self-energy provably exceeds the archimedean on some mid mode for large c
  (would be evidence *against* RH — we'd chase that even harder).

### Vector B — Manifest-positivity factorization τ = M*M (operator/SUSY line)
This session's SUSY lead (A12/A13: τ^½ quasi-local, 98% band energy) says a factorization may exist.
If τ = M*M with M *explicit and zeta-structured*, positivity is manifest and the limit argument may
follow from M's structure.
- **Milestone B1:** compute the Arb Cholesky factor M(N,c); test whether M's entries carry recognizable
  structure (prime positions, log-weights, a difference operator) that persists as c→∞.
- **Milestone B2:** if structured, conjecture the limiting M∞ and prove M∞*M∞ = τ∞ ⪰ 0 directly.
- **Kill if:** M is generic/unstructured across c (no persistence) — factorization is then non-informative.

### Vector C — Uniform-finite negative subspace via a structure function (the linchpin H)
Directly attack Sub-goal H. Candidates to test numerically then prove:
- a de Branges structure function E(z) making the space a canonical H(E) where positivity is a phase
  condition;
- a Sonine / Klein-Gordon-type weight that fixes n_neg;
- Connes' Sonine-space + prolate-spheroidal cutoff picture (his "spectral realization" uses exactly the
  low-frequency prolate modes — which is *precisely* our low-harmonic negative subspace!).
- **Milestone C1:** verify numerically that Connes' prolate/Sonine weighting makes n_neg constant in c
  (our §12 low-harmonic subspace strongly suggests this is the right frame).
- **Milestone C2:** if so, port the finite-core+tail reduction into that space and attempt closure.
- This is the deepest *plausible-to-finish* line; it connects our raw numerics to Connes' existing
  positivity framework without needing anyone else.

### Vector D — Function-field transport (Hodge-index / intersection-signature)
Function-field RH is *proven* (Weil/Deligne) via positivity from the Hodge index theorem on a surface.
Our auto-generated Hodge/intersection probes gesture here. The dream: an arithmetic surface / Arakelov
intersection form whose signature yields our τ's positivity (Connes–Consani F₁ geometry, Bombieri,
Deninger). **Highest ceiling, lowest odds, longest horizon** — a background research line, not the main
bet. Milestone D1: determine whether our τ's inertia (n_pos, n_neg growth) matches a Hodge-index signature
prediction for a candidate arithmetic surface. If the inertia numbers line up, escalate; if not, park.

### Vector E — Rigorous computer-assisted closure (plays to our strength)
Some hard theorems fall to a finite computation once a reduction is in hand (cf. Hales/Kepler, formalized).
If Vector C fixes the core size, E *executes the proof*: rigorously certify the finite core in Arb, bound
the tail, and **formalize the whole chain in Lean** (`RiemannLean/` already exists) so the proof is
machine-checked end-to-end — no referee, no paper needed, the Lean kernel IS the verifier.
- **Milestone E1:** Lean-formalize the finite-core PD certificate (import the Arb bound as a checked fact).
- **Milestone E2:** Lean-formalize the tail bound. E3: the combination lemma ⟹ RH in Lean.

---

## 3. Milestone ladder (proof-gated — each rung needs a *proof pointer*, not a claim)

1. **Kernels explicit** (A1, B1): Q_A(ω), Q_P(ω), M(N,c) in closed/computable form, matched to data. ~weeks.
2. **Rate derived** (A2): analytic c^(−p(N)) — we *understand* the collapse. Big milestone.
3. **Right frame found** (C1): a weighting with n_neg uniformly bounded (test Connes prolate/Sonine first).
   **This is the make-or-break rung.**
4. **Core+tail reduction stated** (§1 + C2 + E1/E2): PD-for-all reduced to a finite certified check + an
   explicit tail lemma.
5. **Closure**: certify the finite core, prove the tail, combine — Lean-checked (E3). = RH.

Each rung either advances or *kills a vector* — no rung is "in progress" forever without a pointer.

---

## 4. Infrastructure (already ours; scale it)

- **Certified compute:** Groskin Arb LDL^T (`external/groskin_2607.02828/code/`), our margin/near-null/
  subspace tools (`operators/*margin*.py`, `neg_subspace_AB.py`, etc.). Oracle for every candidate bound.
- **Fleet:** Olympus GPU (large N,c certifications), HP (parallel sweeps), Mac (dev). Free, no quota.
- **Autonomous engine** (`engine/`): now crux-directed (prime-side, mid-frequency). Runs FREE breadth
  24/7 as a *sensor* — generates & smoke-tests structural probes on the real τ; escalates only real leads.
  Not the prover; the scout.
- **Lean** (`RiemannLean/`): the end-game verifier — makes the finished proof referee-free.
- **Knowledge:** full analytic number theory (explicit formula, Bombieri's functional, subconvexity,
  zero-density), Connes' Sonine/prolate spectral picture, de Branges spaces, RMT as sanity oracle,
  function-field/Deligne as the proven template. All usable; none require contacting anyone.

## 5. How we don't fool ourselves (hardcoded)
- Every "done" needs a proof pointer (commit, Arb JSON, Lean build). Numerics = REPORT; a rigorous
  uniform bound or Lean check = HANDLE. A finite PD certificate is NEVER "RH proved."
- Adversarially re-run and cross-check each claimed bound before it advances a rung.
- Kill criteria are written per vector; a vector with no live milestone gets parked, not nursed.

## 6. Immediate next action (this is where we are on the ladder: rung 1)
**Vector A1 + B1 in parallel:** compute the archimedean and prime self-energies as explicit functions of
frequency across the mid band, AND the Arb Cholesky factor M(N,c), and check both for c→∞ structure.
Then **Vector C1**: test whether Connes' prolate/Sonine weighting makes n_neg constant in c — the §12
low-harmonic negative subspace is a strong hint this is the natural frame, and rung 3 is the make-or-break.
