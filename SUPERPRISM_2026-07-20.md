# SuperPRISM — RH analytic coercivity target (2026-07-20)

Prior state read first (per rule): LEDGER §1-36, OBSTRUCTION_MAP, RH_RESULTS.md (two secured unconditional
theorems, Lean-verified tail lemma), NEXT_SESSION_RH.md, `rh-dbn-criticality-crack` /
`connes-decisive-diagnostic` memories. The retracted DBN-toy route (§34c) and the single open target
(★) A(g)+B(g) ≥ P(g) bandwidth-uniformly are the starting point — NOT re-derived, extended.

Trigger: don't re-map τ (routes 1-5 dead, FE sectors don't isolate the obstruction, real DBN
disconnects from the prime side, §35). Find genuinely NEW angles the 36-section ledger has not touched.

## FRAME — invert the problem ≥5 ways

1. **Extremal-g inversion.** Instead of "prove ∀g," ask what the extreme rays of the admissible
   (even, Paley–Wiener/bandlimited) cone are. If Weil positivity is classically known to reduce to a
   smaller generating/dense subfamily (Krein–Milman for this specific cone), checking that subfamily
   suffices. → classical literature (Barner, Yoshida, Poitou) explicitly studies sufficient conditions
   on restricted test-function families.
2. **Domain-swap inversion.** Position-space kernel comparison (Ω̌(x) vs prime comb) is already
   route 2 — DEAD (prime comb exceeds arch kernel at every prime≥5, pointwise). Confirmed covered;
   not re-run.
3. **"Who already knows and doesn't realize it" inversion.** Classical Weil-criterion literature
   (Weil 1952, Barner 1981 quasi-positivity, Yoshida, Graham–Vaaler explicit majorants, Poitou/
   Odlyzko-te Riele computational tests) — NONE of these names appear anywhere in LEDGER.md (checked by
   grep). This whole classical toolbox is untouched territory for this project.
4. **Pre-mortem inversion.** Every past route failed via LOCAL/pointwise domination (2-4) or a fake
   coercivity injection (§34, circular). Never tried: using the ACTUAL modular/theta structure of the
   prime side (Jacobi theta transformation / Poisson summation, i.e. Riemann's own 1859 derivation
   route) to get an exact identity for a *specific, classical* test-function family, where positivity
   might become a modular-form nonnegativity statement — the Viazovska "magic function" analogy
   (informally noted by Bombieri et al. for Weil positivity, never built here).
5. **Strip-to-pure-problem inversion.** The negative core is a SPECIFIC FIXED window: Ω(r)<0 only for
   |r|<6.29 (§34c, exact). Everything past that is positive and grows like log r. This turns "coercivity
   of an infinite operator" into "does the explicit, bounded, real-analytic function Ω restricted to
   [-6.29,6.29] admit a Positivstellensatz/SOS certificate once combined with B, dominating whatever
   admissible g can put there" — a finite real-algebraic-geometry question, in principle attackable with
   interval-arithmetic-certified SOS (Lasserre/SOS with RIGOROUS, not just numerical, bounds). sos_cert.py
   already computes numerical eigendecompositions (trivial SOS); a certified LOW-RANK version tied to
   the OBSTRUCTION_MAP's "23-dim near-kernel" has not been attempted.

## DIVERGE

### Figure panel (blind, each names a mechanism)
- **Analytic number theorist (Weil-criterion specialist):** "Test the classical Gaussian/Poitou family
  in closed form against your non-circular A/B/P — nobody in this repo has done the canonical
  test-function sweep."
- **Special-functions person (digamma/Gamma asymptotics):** "Ω has an exact Mittag-Leffler / Lorentzian
  series (digamma_poisson_decomp.py already half-built this). Convolve it against ANY test function in
  closed form via the Voigt/Faddeeva identity — this is a general TOOL, not a one-off."
- **Convex-optimization / SOS person:** "The obstruction lives in a bounded, explicit window. Get an
  interval-certified low-rank SOS decomposition of the negative core, using the 23-dim near-kernel as
  the ansatz rank, not full-rank eigendecomposition."
- **Physicist (semiclassics/quantum chaos):** "Ω(0)≈-5.37 is basically a fixed 'binding energy'; ask
  whether the rank-2 pole term B by itself is EXACTLY the Lee-Yang / ferromagnetic compensator for the
  core — check B vs core-restricted P as a standalone finite inequality, no tail needed."
- **Harmonic analyst (Beurling–Selberg / Paley–Wiener):** "Graham–Vaaler give explicit EXTREMAL
  minorants/majorants for the SIGNUM function used in zero-counting proofs — adapt their extremal
  construction machinery (not their majorant conclusion, which is route-2-dead) to build an extremal
  g for (★) directly, via duality/linear programming on the test-function cone."
- **Skip-tracer / residue hunter:** "Where is unused positivity hiding? B(g) is a rank-2 POSITIVE term
  (both eigenvalues of the boundary block are positive — has this ever been checked to EXACTLY cancel
  a piece of the negative core algebraically, term by term, rather than just added numerically?"
- **5-year-old:** "If the bad part only lives between -6 and 6, why not just add up how bad it can ever
  get (one number) and how much good B gives you (one number) and see which is bigger?" → literally
  the finite-core-vs-B calculation above, restated as the simplest possible experiment.

### Technique deck
- **Substitute:** swap the finite bandlimited basis (e_hat_n, box functions) for a genuinely smooth,
  classical test family (Gaussian) with closed-form transforms on both sides.
- **Combine two dead ideas:** route 2 (pointwise Beurling-Selberg, dead) + the DBN tail-coercivity
  finding (§34, real part) → don't add an artificial flow; use natural Gaussian WIDTH as the
  "coercivity dial" instead of a fake heat-kernel parameter — no circularity because A(g_w) is directly
  the quadratic form evaluated on a REAL admissible g, not an added diagonal.
- **10×/0.1×:** push w (Gaussian width in frequency) to extremes — w→0 concentrates almost entirely on
  the negative core (never isolated this cleanly before); w→∞ probes deep into the coercive tail.
- **Follow-the-residue:** B(g)=2exp(1/4w²) blows up as w→0 exactly when P(g) also grows (many primes
  contribute) — check the RATE race between B's exponential growth and P's growth as w→0; this is a
  genuinely quantitative, closed-form question.
- **Constraint-break:** classical Weil-criterion papers use test functions with COMPACT support in x
  (needed for the original Weil 1952 setting); our operator's admissible class is bandlimited in r
  (Paley–Wiener dual). Gaussians satisfy neither exactly but approximate both — check whether the
  admissibility gap itself is where danger could hide (a technical audit item, not a research item).

### Banded spectrum
- **[GROUNDED]** Gaussian test-function family, closed-form A via quadrature + Voigt cross-check,
  closed-form B, direct P — sweep w from deep-core to deep-tail. *(Selected — executing now.)*
- **[CLEVER]** Finite-core Positivstellensatz: restrict to |r|<6.29, get an interval-certified SOS
  bound for Ω's core contribution combined with B, compared against the WORST-CASE core-restricted P
  (bounded because g bandlimited ⟹ P's contribution from the core alone is controllable). Sane cousin
  of "prove uniform coercivity" — proving it just on the finite dangerous window is a strictly smaller,
  well-posed target.
- **[CLEVER]** Check whether τ's near-null/negative eigenvectors (23-dim kernel, OBSTRUCTION_MAP) are
  spanned (numerically, then symbolically) by the first ~23 Lorentzian terms P_n(r) from the
  Mittag-Leffler decomposition — if so, the obstruction has an explicit LOW-DIMENSIONAL basis, a genuine
  structural simplification (not yet checked; digamma_poisson_decomp.py stops short of this).
- **[EXOTIC]** Adapt Graham–Vaaler's extremal-function LP/duality machinery (built for Beurling–Selberg
  majorants) to construct the extremal g for (★) directly via linear programming over Fourier
  coefficients — turns "prove ∀g" into "solve one (possibly infinite, but truncatable) LP," with
  infeasibility/optimal-value sign as certificate. Sane cousin: run the LP at finite bandwidth (same
  λ²=14 basis already in sos_cert.py) as a warm-up sanity check before attempting the infinite version.
- **[MOONSHOT]** Viazovska-style "magic function": seek an explicit modular-form-derived g (built from
  theta series / Eisenstein series under Jacobi's transformation) for which A(g)+B(g)-P(g) has an exact
  modular-form-coefficient interpretation whose nonnegativity is a KNOWN modularity fact (e.g. a cusp
  form Fourier coefficient bound, or an Eisenstein series positivity). No concrete construction yet;
  flagged as the most theory-native "new math" direction if the Gaussian family reveals a modular
  signature (e.g. approaches an exact identity at special w tied to θ-function self-duality points).
- **[ABSURD]** Genetic/evolutionary search over a large parametric test-function family (splines,
  wavelets) to hunt for a numerical near-violation of (★) as an ADVERSARIAL check — sane cousin: this
  is exactly what beurling_adversary.py already does at finite N; extend its search space to the
  Gaussian-plus-perturbation family once the closed form exists, as a stress test of the new machinery.

## CONVERGE

**Winner, executing now:** the Gaussian test-function family (`operators/gaussian_family_positivity.py`).
Rationale: fully non-circular (only Ω, the rank-2 pole term, and the prime sum — no zero ever enters,
matching connes-decisive-diagnostic's non-circularity discipline), closed-form on two of three terms,
classical (ties to real literature nobody in this repo has invoked), and — critically — it is the FIRST
tool in this project that can dial continuously into the core-dominated regime (w→0) in a single sweep,
rather than being stuck at whatever the finite bandlimited basis happens to resolve.

**Queued next (not yet run, in priority order):**
1. Near-null-mode / Lorentzian-basis rank check (CLEVER) — cheap, fast, directly answers whether the
   23-dim kernel has an explicit closed-form spanning set.
2. Finite-core Positivstellensatz / interval-certified SOS on |r|<6.29 (CLEVER) — the well-posed,
   strictly-smaller-than-RH target; natural next step if the Gaussian sweep sharpens where exactly the
   margin is tightest.
3. Graham–Vaaler LP-duality extremal-g construction (EXOTIC) — heavier lift, worth scoping once (1)-(2)
   report back.
4. Magic-function/modular search (MOONSHOT) — parked until the Gaussian sweep either shows a modular
   fingerprint (a special w with exact cancellation) or doesn't.

Killed / not re-run: DBN toy (circular, §34c), Loewner/Pick (=route 4, §30), moment/Hankel (dead, route 1),
pointwise Beurling-Selberg in x-space (dead, route 2), pointwise-3T² prime √-phase sum (research-hard,
parked at §31/§33).
