# SuperPRISM 2 — RH analytic coercivity target (2026-07-20, second run)

Second full prism on the same target, run same-day as `SUPERPRISM_2026-07-20.md` (prism 1), deliberately
from a different register. Pre-read (search-before-build, verified by grep): LEDGER §1–36, OBSTRUCTION_MAP,
ATTACK_QUEUE (A5 queued-not-run; A6 DONE), RH_RESULTS.md, VECTOR_D_FUNCTION_FIELD.md,
uniform_positivity_target_2026_05_19.md, nyman_beurling/ (v9 definitive, dead-end), the 07-18 physics
deep-dive ("Lee–Yang is the wall in physics clothing"), structure-factor/Dyson correction. Confirmed ZERO
ledger hits: Borcea–Brändén, Kurasov–Sarnak, Alon–Cohen–Vinzant, total positivity, Pólya frequency
sequences, Knauf, Lagarias half-plane positivity, family-averaged positivity. Confirmed:
`auto_jensen_signature.py` / `auto_selberg_spectrum.py` are corrupted placeholder stubs, not prior results.

## The two constraints every prism-2 idea must pass (new meta-rules)

**(M1) The A6 filter.** ATTACK_QUEUE A6 (DONE 2026-07-18): τ_χ for L(s,χ₄) sits EXACTLY on the PSD
boundary (min eig ~−3e-14) — no pole, different Γ-factor, same knife-edge. Positivity is UNIVERSAL across
the L-family. Therefore any winning mechanism must be family-generic. ζ-specific stories (pole-term B as
the compensator; ζ-only magic functions) are structurally disfavored *by our own ledger*.

**(M2) The quantifier/certificate rule** (pre-mortem: every failed route was a re-coordinatization of ★
attacked with an analytic inequality). A prism-2 winner must change at least one of:
- the QUANTIFIER structure ("for all g / all n" → "on average / for almost all / over a family"), or
- the CERTIFICATE TYPE (analytic inequality → combinatorial bijection, algebraic identity, or
  classification theorem).
A new PROBE (new test family, new basis) is cartography — allowed, but never counted as a proof path.

## FRAME — five fresh inversions

1. **Classify-the-solution-space (Dyson made rigorous).** Don't prove ζ's measure is nice — use the fact
   that "nice" now has a classification. Kurasov–Sarnak (2020) construct positive Fourier quasicrystals
   (FQs) from Lee–Yang polynomials; Alon–Cohen–Vinzant (2024): EVERY 1-D positive FQ is the zero set of a
   Lee–Yang polynomial restricted to a torus line. Under RH the zero-counting measure Σ_ρ δ_γ is a positive
   crystalline measure whose spectrum is {±k·log p} — and the Euler product is literally a function of
   torus variables z_p = p^{−it}, one per prime. Dyson's 2009 "classify the quasicrystals" program is no
   longer a metaphor: the finite-frequency case is a THEOREM, and the explicit formula puts ζ exactly at
   its infinite-variable limit. The ledger's own Lee–Yang finding ("ξ has the reflection symmetry, lacks
   coefficient positivity" — A5) now has a target-shape theorem on the other side of the wall.
2. **Operator-preserver inversion (handle (ii) modernized).** Handle (ii) = Hermite–Biehler for ξ̂ with no
   zero locations. Borcea–Brändén (2009) classified ALL linear operators preserving stability/real-
   rootedness, with a finite SYMBOL TEST (T preserves stability ⟺ the symbol T[(z+w)^n] — or the
   transcendental exp-symbol — is stable). ξ(½+iz) = ∫Φ(u)cos(zu)du is the image of the positive theta
   kernel Φ under an explicit linear transform. Pólya's old program ("cosine transforms of positive kernels")
   stalled for lack of exactly this classification. Reframe: don't prove real-rootedness of ξ̂ — compute the
   BB symbol of the theta cosine-transform (or of the Riemann–Siegel truncation operators) and identify
   PRECISELY which BB stability condition it fails/needs. Worst case: the sharpest restatement of handle
   (ii) that exists. Best case: the transform factors through a known preserver class.
3. **Family-average inversion (change the quantifier, keep the form).** Our two secured theorems live on
   the n-axis (all n ≤ 2πT; almost-all n ≤ 3T²). A6 says the χ-axis carries the same knife-edge. Average
   (★_χ): A_χ(g) ≥ P_χ(g) over all characters χ mod q: orthogonality kills all primes except p^k ≡ ±1
   (mod q) — the prime side THINS by a factor ~φ(q), while the archimedean side barely moves. The averaged
   inequality should be provable UNCONDITIONALLY for bandwidth up to ~log q, and the individual-vs-average
   deviation is controlled by the large sieve (route 3's tool, now pointed along the χ-direction where it
   is strong, instead of the n-direction where it hit the wall). Candidate "Theorem 3": positivity of the
   Weil form for ALMOST ALL χ mod q, unconditionally — a third finite theorem on a NEW axis.
4. **Pre-mortem inversion.** "We already failed — why?" Because probes multiplied but the certificate type
   never changed. Encoded as meta-rule (M2) above; prism 1's winner (Gaussian family) is hereby graded
   honestly: cartography of the core-dominated regime, not a proof path even in principle (a 1-parameter
   slice of ∀g).
5. **Hider's-hand / LP-duality inversion.** Universality + exact knife-edge (margins → 0 at every scale, in
   every family) is the signature of an optimization at its optimum with complementary slackness: (★) as
   the infinite LP "minimize (A+B−P)(g) over normalized admissible g." RH ⟺ optimum = 0; the near-null
   modes we keep measuring ARE the optimal primal; the zero-measure Σδ_γ IS the dual certificate; Weil's
   explicit formula IS strong duality. Proving RH = exhibiting the dual certificate without knowing the
   zeros. Not a method by itself, but it explains WHY every numeric shows collapse-to-0 margins (an optimum,
   not a coincidence) and suggests perturb-primal→track-dual experiments (thin the primes ε, watch the dual
   measure move — quantitative prime–zero repulsion).

## DIVERGE — panel (eight minds, disjoint from prism 1's seven)

- **Algebraic combinatorialist (LGV/total-positivity):** "ξ ∈ Laguerre–Pólya has an Aissen–Schoenberg–
  Whitney/Edrei-side characterization: total positivity of the coefficient Toeplitz matrices. The 2×2
  minors are the Turán inequalities — PROVEN unconditionally for ξ (Csordas–Norfolk–Varga 1986); 3rd-order
  Turán also proven (Dimitrov–Lucas 2011). So there is an EXISTING unconditional ladder k=2 ✓, k=3 ✓,
  general k ⟺ (with reality) RH — and nobody has mapped where the finite minors get TIGHT, the way this
  repo mapped τ's margins. Better: total positivity has combinatorial certificates — Lindström–Gessel–
  Viennot planar networks. A planar-network representation of ξ's (or any completed L's) coefficient
  minors would prove positivity by BIJECTION, a certificate type never tried on ★ — and Euler products
  are born as products over primes, i.e. natural network layers." *(Passes M1: works for any completed
  L-function. Passes M2: new certificate type.)*
- **Stability theorist (Borcea–Brändén):** frame 2 above; deliverable = the computed BB symbol of the
  theta transform, and the multiplier-sequence test (Pólya–Schur) for the Taylor-coefficient sequence of Ξ.
- **Analytic number theorist (families):** frame 3 above; deliverable = the χ-averaged explicit formula
  with the orthogonality-thinned prime side, plus a lit-check (Murty et al. have explicit-formula-on-
  average results; verify what's genuinely open before proving).
- **Spectral geometer (control experiment):** "Selberg zeta for a compact hyperbolic surface SATISFIES its
  RH (the Laplacian exists — positivity is a theorem). Build the SAME τ = A+B−P object there (geodesic
  lengths as 'primes,' Weyl term as archimedean symbol) and diff the two mechanisms term by term. In the
  solved world, WHAT carries the positivity that ζ's world lacks? That difference is a precision pointer
  at the missing structure — a control experiment this repo has never run (the existing Selberg script is
  a corrupted stub)." *(Cartography, but of a new kind: comparative anatomy against a SOLVED sibling.)*
- **Statistical mechanic (Knauf):** "The Farey/number-theoretic spin chain (Knauf 1993+) is PROVEN
  ferromagnetic, with free energy tied to ζ(s−1)/ζ(s). It is the one physical model where the coefficient
  positivity Lee–Yang needs is a theorem. If A5 ever runs, run it THERE first, where the missing hypothesis
  holds." *(Queued behind A5; integrates with, not duplicates, the 07-18 physics dive.)*
- **Arithmetic geometer (Vector D bridge):** "ACV's classification realizes positive FQs as zero sets on
  torus-line closures — a torus with one circle per frequency. For ζ that torus is Π_p S¹, which IS the
  idele-class-space skeleton Connes uses, and a Lee–Yang polynomial in variables {z_p} is a function on
  exactly the object Vector D says is missing (Spec ℤ ×_{F₁} Spec ℤ's 'points' = pairs of prime phases).
  The FQ classification may be the combinatorial shadow of the missing surface." *(Moonshot; grafts prism-2
  frame 1 onto the already-scoped Vector D dictionary rather than opening a new lane.)*
- **Historian (Riemann's Nachlass):** "The Riemann–Siegel formula sat unread for 70 years because everyone
  theorized instead of computing what Riemann actually did. The modern echo: the Z-function truncations —
  non-harmonic cosine sums Σ n^{−1/2}cos(θ(t)−t log n) with POSITIVE coefficients — are empirically
  near-real-rooted. Real-rootedness of positive-coefficient non-harmonic cosine sums is precisely
  Hermite–Biehler territory (frame 2's toolkit). Compute where their stability actually fails as the
  truncation grows; that failure locus is handle (ii) made finite." 
- **Adversarial statistician:** "Positive-coefficient Nyman–Beurling: rerun the dead v9 machinery
  restricted to NONNEGATIVE dilation coefficients. If the Báez-Duarte distance still →0, cone-positivity
  is free and NB links to the Weil cone; if it stalls, you've localized exactly what signed cancellation
  buys — either way a cheap, sharp measurement on infrastructure that already exists." 

### Technique-deck draws (kept)
- **Combine-two-dead:** NB (dead) × Lee–Yang wall (positivity-missing) → the positive-coefficient NB probe
  above.
- **10×:** push the TP-minor ladder to k ~ 40–100 in high precision (arb/mpmath) — margin cartography in a
  basis where k=2,3 are unconditional THEOREMS, so tightness-vs-k directly measures distance from known
  ground to RH.
- **Assumption-reversal:** rather than proving ξ ∈ LP, compute what LP-membership FORCES on the thinned
  χ-averaged prime side (feeds frame 3).

## CONVERGE — graded queue

> **EXECUTION UPDATE 2026-07-20 (Opus, verify-before-ship).** Item 1's CORE was already done — see
> ATTACK_QUEUE **A14** (`operators/toeplitz_tp_test.py`): the Ξ-coeff Toeplitz matrix IS Totally
> Positive to n=15 (Neville all-positive), TP2 ratios 2.15→1.02 decreasing to 1 from above. Crucially
> A14's payoff (Euler product → TP closure) **FAILS** for exactly the reason M1 predicts: ξ=arch×ζ but
> ζ(½)<0, so the per-prime factors are not locally positive — positivity is product-EMERGENT, not
> factorable. This is my own prism's meta-rule confirmed by a result already in the repo. **Self-
> correction: I proposed a largely-completed item without checking ATTACK_QUEUE's A14 line first**
> (I greped "Toeplitz/total positiv" against LEDGER+STATUS but not A14's outcome in ATTACK_QUEUE — the
> SAME miss I flagged prism 1 for on A6). The genuinely-open residue of item 1 (still worth doing, not
> yet done): the **LGV planar-network combinatorial certificate** and the **variation-diminishing /
> sufficient-conditions** TP handles A14 explicitly lists as "unexhausted." But A14+A6 together now make
> item 1's whole family (local/factorization certificates) structurally disfavored — so item 1 is
> DEMOTED below item 2. Meta-lesson banked to PROCESS.md: prism pre-reads MUST include ATTACK_QUEUE
> A-item outcomes, not just LEDGER/OBSTRUCTION_MAP.

1. ~~TP/Turán minor ladder~~ **[CORE DONE = A14]** — residue: LGV planar-network certificate +
   variation-diminishing handle only. Demoted (A14+A6 disfavor the whole local/factorization family).
2. **χ-family almost-all positivity ("Theorem 3" candidate)** [THEOREM-SHAPED, days not hours] — write the
   χ-averaged Weil form, thin the primes by orthogonality, bound deviations by the large sieve along the
   χ-axis. Lit-check first. Would pair with the two n-axis theorems as a third finite unconditional result
   on a NEW axis. Passes M1 (it IS the family) + M2 (quantifier change).
3. **BB symbol of the theta transform** [STRUCTURAL, research-grade but finite] — compute the Borcea–
   Brändén / Pólya–Schur symbol for the cosine-transform-of-Φ operator and for RS truncations; state
   exactly which stability condition is missing. Sharpest-available form of handle (ii). Passes M1+M2
   (classification-theorem certificate).
4. **Selberg-zeta control experiment** [CARTOGRAPHY, new kind] — comparative anatomy vs a solved sibling.
5. **Positive-coefficient NB probe** [CHEAP CURIOSITY] — one evening on existing v9 code.
6. **ACV ⇄ Vector D bridge** [MOONSHOT] — parked until 1–3 report.
Killed for prism 2: anything ζ-specific (M1), any new probe sold as a proof path (M2), Knauf before A5,
re-running structure-factor numerics (done), NB as-is (dead).

---

## Comparison: prism 1 (Sonnet, this morning) vs prism 2 (Fable, now)

**What prism 1 did well.** Correct frame discipline (didn't re-open routes 1–5); found one genuinely
useful new PROBE (the Gaussian/Voigt family — first tool that dials continuously into the |r|<6.29 core;
scan currently running); a clean queue (near-null Lorentzian rank check is still worth doing); flagged
Graham–Vaaler/classical-literature blindspot correctly.

**Where prism 1 fell short — caught by prism 2's re-read of the ledger:**
- **Missed the A6 filter.** Prism 1's "skip-tracer" idea (B exactly compensates the negative core) is
  undercut by our OWN ledger: L(s,χ₄) has NO pole term B and sits on the same knife-edge (A6, 07-18).
  Its "magic function" moonshot as stated was ζ-specific — same filter. Prism 1 was written without
  re-reading ATTACK_QUEUE's A6 outcome; prism 2 promotes it to a standing meta-rule (M1).
- **No quantifier/certificate discipline.** Prism 1's winner is a probe; prism 1 nowhere asks "can this
  category of idea even terminate in a proof?" Prism 2's M2 rule grades every idea on that axis first.
- **Shallower literature reach.** Prism 1 surfaced classical names (Barner, Yoshida, Poitou, Graham–
  Vaaler). Prism 2 surfaced four modern, ledger-untouched toolkits with classification-strength theorems
  behind them: Borcea–Brändén preservers (2009), Kurasov–Sarnak FQs (2020), Alon–Cohen–Vinzant (2024),
  and the ASW/Edrei total-positivity ladder with its two already-proven unconditional rungs.

**Independent convergence (the signal PRISM looks for).** Both prisms, from disjoint panels, converge on
the Hermite–Biehler/stability axis — prism 1 via "handle (ii)" and the magic-function instinct, prism 2
via BB preservers, FQ classification, and the TP ladder. Both independently reach for LP/duality (prism 1:
Graham–Vaaler extremal LP; prism 2: knife-edge = complementary slackness). Read: the stability/HB axis is
where the two runs agree the treasure is buried; the duality view is the shared explanation of WHY every
margin we ever measure collapses to zero.

**Net verdict.** Prism 1 = probe-and-map (better cartographer). Prism 2 = classify-and-certify (better
strategist: meta-rules M1/M2, a theorem-shaped target on a new axis, and certificate types this project
has never fielded). Execution order: finish prism 1's running Gaussian scan (cartography, cheap), then
prism 2's #1 (TP ladder — anchored by two known unconditional theorems), scope #2 (Theorem 3 lit-check),
with prism 1's near-null Lorentzian rank check kept alive as the bridge between the two (its 23-dim
near-kernel is the primal-optimal object in prism 2's duality frame).
