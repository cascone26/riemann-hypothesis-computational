# The RH Corner: a localization map (Riemann8 session, 2026-08-09)

*Honest consolidation of LEDGER §884–§926. Every quantitative claim below is HANDLE (code-run, reproducible)
unless marked REPORT. The purpose of this document is to state precisely (1) what the remaining content of RH
has been reduced to, (2) the complete map of tactics ruled out this session and why, and (3) the honest open
core. It banks a real result — the sharpest localization we have — without overclaiming a proof.*

---

## 1. The setting

Weil's criterion: RH ⟺ the explicit-formula quadratic form **Q(φ) = pole(φ) + arch(φ) − prime(φ) ≥ 0** for all
even band-limited test functions φ, where on a Fejér band-limited span (window `[a,b]`, bandwidth `4W`, `N`
centers) Q is an explicit `N×N` generalized-eigenvalue problem against the Gram metric `G`:

- `arch` = compression of multiplication by `Ω(r)/2π`, `Ω(r) = Re ψ₀(¼+ir/2) − ln π` (archimedean density);
- `prime` = `(1/π) Σ_{log n < 4W} Λ(n)/√n · [compression of cos(r log n)]` (exact for bandwidth-`4W` φ);
- `pole` = rank-≤2 PSD term from ζ's pole at s=1.

**RH-corner** = the regime where the least generalized eigenvalue (the *margin*) is razor-thin.

## 2. The localization chain (§884–§912)

Each step is verified and reduces the object:

1. **Q = M_σ + pole, exactly** (§889). `M_σ = arch − prime` is the compression of multiplication by the *exact*
   symbol `σ(r) = Ω/2π − P(r)` (‖arch−prime − M_σ‖ = 1e-15). The pole is a genuine rank-≤2 term, not an error.
2. **One negative eigenvalue** (§889–890). At the corner `M_σ` has *exactly one* negative eigenvalue `−μ`;
   `#neg(M_σ) ≤ rank(pole)` is forced by RH; the pole lifts it. `σ` is pointwise-negative on 7–40% of every
   window, yet `M_σ` has a single negative direction — the arithmetic near-coincidence of arch and prime.
3. **Scalar form** (§891). `RH-corner ⟺ g = 1 − μ·wᵀA⁺w ≥ 0`, `A = M₊ + p·vvᵀ`, leading form `p(v·w)² ≥ μ`.
4. **The razor is off-diagonal / near-null** (§900–903). The leading `π_w/μ` has ~4% slack, but the M₊ near-null
   cluster eats it: the true razor is an irreducibly **~2–3-dimensional coupling** in the near-null (prolate-
   plunge-related, ~0.53 overlap, *not* clean Landau–Widom) subspace. No single scalar closes it.
5. **Dual / zeros picture** (§908–909). The margin `= 2 Σ_γ w(γ)²`. The worst direction `w` **annihilates its
   in-window zeros exactly**; the entire margin is `w`'s Paley–Wiener tail spilling onto the zeros *just above*
   the window (in-window 0%, below 0%, above 100%).
6. **Edge effect** (§910). The corner exists only at the **bottom edge of the zero spectrum**, where spillover is
   one-sided (upward). The moment a window has zeros below it, the negative eigenvalue vanishes; the interior is
   robustly positive.

## 3. The invariant statement (what everything reduces to)

> **RH ⟺ an unconditional lower bound on the band-limited Weil margin, uniform in bandwidth `4W` and height `t`.**

The strategic survey (§920, five programs, source-verified) confirmed **every** major RH program reduces to this
same statement: Connes–Consani (Weil positivity, the CCM μ_λ≥0 prolate lemma), de Branges (the ω→0 canonical-
system extension), function-field/Hodge (intersection-form positivity), and the positivity criteria (Li /
Bombieri–Lagarias / Nyman–Beurling). They are the same wall from different angles.

## 4. Tactics ruled out this session (the elimination map)

| Tactic | Verdict | Reason |
|---|---|---|
| Scalar / operator bound `floor − pnorm` | dead (§914, §919) | fails at the corner; misalignment correction closes only ~75% |
| Finite bottom region | tempered (§911) | corner region grows with `4W` (not a fixed height range) |
| Concentration / Szegő (band-limited) | robust-region only (§888) | diagonal symbol wrong-sign at corner; off-diagonal effect |
| Low-prime / first-few-primes | dead (§917) | razor is a global cancellation over the *full* prime set |
| Two-term `arch ⪰ p=2` (survey's #1 pick) | **refuted (§925)** | fails at `4W≥6`; primes cannot be peeled — cancellation is essential |
| Density-preserving (PNT main + fluctuation) | bottoms out (§926) | smooth part vacuous at the window; content = the error term = RH |
| de Branges *global* positivity | blocked (§920) | Conrey–Li 2000 counterexample at the 34th zero |
| Function-field / Hodge translation | blocked (§920) | Deninger 2022: no real-coefficient arithmetic Weil cohomology exists |
| Arb-certification of the corner | radius-blocked (§894) | razor margin (~1e-15) sits below achievable interval radius |

## 5. What IS proven / certified (the real, banked positives)

- **Robust region, unconditional** (§875): `floor(t) > pnorm(4W) ⟹ Q ⪰ 0`, with the exact floor law
  `floor(t) ≈ (1/2π)ln(t/2π)` (arb-certified). Certified frontier extended to `4W ≈ 5` (§871–873).
- **Critical-subspace positivity, arb-certified** (§884): `Q ⪰ 0` on the strong-prime subspace at a deep corner,
  min pivot `2.3e-5 ≫` radii — a rigorous zero-free certificate that arch dominates prime on the critical subspace.
- **Closed-form archimedean weight** (§898): `Ŵ(u) = −π e^{u/2}/sinh(u)`, verified; arch/prime pointwise ratio
  `= 2π/log n`, crossing 1 at the scale `u = 2π`.
- **Archimedean density is Cauchy–Stieltjes** (§897): `dΩ/d(r²)` completely monotone, unconditional from ψ₀'s series.
- The Li coefficients are comfortably positive (`λ_n ~ (n/2)ln n`, §921) — the razor is specific to *band-limiting*.

## 6. The honest open core

The remaining mathematical content of RH, localized: **on the bottom-edge, band-limited critical subspace, the
archimedean form dominates the prime form — uniformly in bandwidth and height — by a margin the pole exactly
lifts to zero.** This is a ~2–3-dimensional near-null (prolate-plunge) positivity that is a *global* prime
cancellation, razor-thin, sub-interval-radius, and identical (verified) to the open lemma every leading program
leaves unproven. It is **universal across L-functions** (§928, §930): each L is razor-thin at its own spectrum
bottom edge, so any proof via this route must be universal — it cannot use ζ-specific arithmetic (the pole, the
specific low primes).

**The nature of the barrier (§931, the sharpest statement).** The function-field case, where RH *is* proven, has
the *identical* Weil form — a sum-of-squares over the zeros (a Gram matrix `VV*`), PSD ⟺ zeros on the line, with
the *same* sampling razor (min eigenvalue → 0 as the band-limit reaches the zero count; negative if a zero moves
off the line). The **only** difference is the proof of zeros-on-line: function fields supply it via the
**Hodge-index geometry** (intersection form on `C×C`), which hands them the Gram/PSD structure for free. Number
fields have the same zero-free form (`arch+pole−prime`, ⟺ RH) but **lack the surface/Hodge-index** — Deninger's
own 2022 no-go rules out a real-coefficient arithmetic Weil cohomology. **So the razor is the analytic shadow of a
missing geometric positivity, and no analysis of the arch-vs-prime balance can supply it.** This is why every
analytic lever this session bottomed out: the missing ingredient is structural/geometric (a number-field
Hodge-index / Hilbert–Pólya / Connes-arithmetic-site positivity), not analytic. That is the genuine frontier.

**This map is the contribution: RH's core is now an explicit, finite, reproducible object; a complete account of
why each classical and modern route stalls at exactly it; a universality constraint on any solution; and — via the
function-field comparison — a precise identification of the barrier as the missing geometric structure rather than
an unproven analytic inequality.**

*Provenance: LEDGER §884–§926; PROCESS.md arc; ANCHOR_PAPER.md §8–§10; workflows wf_045cb265, wf_5685d3c6,
wf_e696c870, wf_259382b8. Discipline: every agent-sourced number independently re-verified; multiple over-claims
(mine and agents') caught and corrected in-ledger (§888, §889, §897b, §902, §918, §922→925).*

## 7. Where to go next (2026-08-09 assessment; leverage × tractability)

Ranked, honest. The two things that scale-and-speed genuinely help are the partial wins; the real prize is a theory bridge.

**A. The condensed-math bridge to the missing geometry (RH itself — the real prize, deep theory).**
The barrier is a *real-coefficient* obstruction (Deninger 2022 no-go: no real-coefficient arithmetic Weil cohomology).
Clausen–Scholze condensed / analytic geometry (solid/liquid real vector spaces, analytic stacks) is the framework
built precisely to make real/analytic coefficients behave homologically — it re-proved Serre duality, GAGA, Riemann–
Roch "analysis-free." **No published work connects it to arithmetic Weil positivity / ζ.** The precise open question:
does a Weil cohomology with *condensed/solid* real coefficients for Spec ℤ exist, evading Deninger's no-go (which is
about *ordinary* real coefficients)? This is the genuinely unexplored, highest-leverage direction — but it is a
sustained theory program, not an accelerable computation. **Do:** a deep, source-verified feasibility study of whether
Deninger's argument is coefficient-agnostic or specifically breaks for condensed coefficients.

**B. A true universality law (accelerable partial win — but likely = §875 + conductor).**
Compute the actual first zeros of a few hundred L-functions (LMFDB / mpmath, embarrassingly parallel), place each L's
edge precisely, and extract the edge-margin invariant. HONEST CAVEAT (§933): the mechanism is the universal floor-vs-
pnorm balance with the conductor as an additive floor shift (floor += ln|d|/2π), so the "law" is likely §875
universalized, not a new invariant. Still a clean, publishable universality statement (ζ = the hardest, lowest-floor case).

**C. Massive certified-frontier extension (accelerable, real, unconditional — but does not reach the corner).**
Arb-certify Q ⪰ 0 over a large explicit (4W, height) grid at scale. A genuine unconditional partial theorem. The corner
stays out of reach (sub-interval-radius razor, §894).

**What does NOT work (established, not guessed):** brute-forcing the razor (§884–930), any prime-decomposition (§917,
§925, §926, §929), scalar/operator bounds (§914, §919), SDP/Lasserre (boundary-critical), closed-form SOS (§934),
comfortable criteria (Li/Báez-Duarte are infinite-limit, §921/§927). All bottom out at the same semi-local positivity = RH.

---

## 8. The Π₁ / Krein–Langer reframing — the new frontier (2026-08-09, LEDGER §946–§954)

*Added after the §945 fragility filter selected "encode RH as a robust integer." Executed and cross-verified.
Every quantitative claim is HANDLE (own code, reproducible); literature claims are REPORT (3 independent
source-cited recon agents). This section supersedes §7's ranking as the session's live lane.*

**The reframing (HANDLE, two independent realizations).** The §944/945 finding that the Weil MARGIN is fragile
(knife-edge, generically negative) motivated moving to the INERTIA — an integer, robust under perturbation.
Result: **RH ⟺ #neg(Weil form Q) = 0**, and the entire burden is that the arch−prime part **M_σ = arch − prime
carries exactly ONE negative square** (κ=1) that the ζ-pole term absorbs.
- §946 (function-field, RH_ff PROVEN ground truth): #neg(T_N) is a topologically robust integer — 0 iff RH,
  insensitive to violation magnitude (margin swings 5 orders while inertia stays fixed). The inertia is the
  protected integer the §945 filter demanded; Hodge-index is what computes it =0 there.
- §947–948 (number field, ζ): #neg(Q)=0 across all windows; #neg(M_σ)=1 ALWAYS (never more), despite the symbol
  σ=Ω/2π−P being negative on 28–43% of each window (Szegő predicts ~fN negatives; the defect grows to ~13).
- §953: κ≤1 survives ARBITRARY sampling (irregular/over-dense-4×/clustered/wide/high-band) ⟹ M_σ is a genuine
  Pontryagin **Π₁** reproducing kernel (Krein–Langer sense: #neg≤1 for every finite Gram).
- §954: **exact** (compression-free Suzuki closed-form assembly) — #neg(Q)=0 and κ=#neg(Q−pole)=1 confirmed
  exactly across a,K. Bulletproofs the finding against any quadrature/band-limit doubt.

**It is genuine RH-content, not generic (HANDLE, §949).** Only the REAL von Mangoldt primes pin κ=1 uniformly;
EVERY perturbation — even δ=0.001 on the log-primes — makes κ grow with bandwidth (→9). Matches the classical/
Groskin fact #neg(truncated Weil form)=½·#off-line-zeros: κ IS the RH-violation counter.

**Why every prior route bottomed out here — the PARITY diagnosis (REPORT, Deninger 2204.02714 re-read).**
Deninger's no-go is NOT a positivity/margin obstruction — it is a PARITY/EXISTENCE one: functional-equation root
numbers force ODD-multiplicity symplectic Galois structure, but complexifying a REAL vector space forces EVEN
multiplicity → odd≠even. Coefficient-SPECIFIC to ordinary real vector spaces (Deninger himself proposes
quaternionic/antilinear-involution coefficients as the escape; realized by Clausen–Scholze condensed/solid
modules; NO published bridge to Spec ℤ ζ). **Unifying hypothesis:** the ONE (odd) negative square of M_σ IS the
analytic shadow of this odd-multiplicity symplectic parity obstruction. An odd negative square cannot be resolved
in a real/positive (Π₀, Hilbert) framework — which is exactly why (a) de Branges' Hilbert-space approach failed
(Conrey–Li 2000, counterexample at the 34th zero — the object is Π₁ not Π₀) and (b) every real-analytic lever
this session bottomed out. The correct framework is a Pontryagin Π₁ space (Krein–Langer / Kaltenbäck–Woracek) —
the analytic MIRROR of the quaternionic/condensed escape. Sharpens §7's "analytic shadow of missing geometry"
to a precise PARITY statement.

**The concrete route, with its one named gap (REPORT, deep-extraction agent).** Krein–Langer: a kernel is Π_κ ⟺
K=K₊−K₋ (both PSD, rank K₋=κ). Generalized-Nevanlinna theorem (arXiv:1011.2081): a Nevanlinna function g with a
single simple pole off ℝ has Nevanlinna kernel K_g=(g(z)−ḡ(w))/(z−w̄) with EXACTLY one negative square. ζ's pole
at s=1 is a simple pole. Route: (1) **[the gap]** reformulate M_σ as the Nevanlinna kernel of an explicit ζ Weyl/
m-function g; (2) count g's UHP poles = 1 (off-line zeros would add poles) ⟹ κ=1 by Krein–Langer; (3) the pole
term absorbs the one negative square ⟹ Q⪰0 ⟹ RH. Steps 2–3 are existing machinery invoked by name; step 1 (the
Weyl/Nevanlinna m-function construction) is the crux — and is the SAME "technical m-function duality" left OPEN
in the Krein-string lane (§511–608). Existing assets to attack it: `suzuki_exact_Q.py` (exact screw assembly),
`krein_string_zeta.py`/`weyl_*` (canonical-system machinery), Suzuki 2206.03682/2606.09096 (screw function =
RH-equivalent positivity, canonical-system-ready).

**Honest status.** "M_σ is globally Π₁" ⟺ RH (a repackaging). The value is the TOOL-TRANSFER: a real positivity
inequality (which bottomed out) becomes a Nevanlinna POLE-COUNT (argument principle) — a complex-analytic handle
with its own machinery (Krein–Langer pole-counting, canonical systems). Not a proof; a genuinely new, verified,
literature-novel attack surface with a precisely-named remaining gap. **Next construction:** build the ζ Weyl
m-function from the exact screw assembly and test the pole-count = κ correspondence directly.

### 8b. Cross-Selberg confirmation: κ = pole order (LEDGER §958–958b, HANDLE)

A new falsification test of "the one negative square IS ζ's pole," independent of §928 razor-universality:
build the Weil form arch−prime for real Dirichlet characters (all ENTIRE — no s=1 pole) and compare κ to ζ.
Result (each L at its own bottom edge, bandwidth sweep, correct primitive characters + auto-parity from χ(−1)):
- **ζ (simple pole): κ = 1, capped and robust across bandwidth (Π₁).**
- **Every entire χ (mod 3,4,5,7,8,11): κ = 0** — genuinely PSD, no negative square to lift.
So **κ(arch−prime at the spectrum's bottom edge) tracks the presence of the s=1 pole** across the Selberg class:
ζ's single negative square is the pole; entire L-functions have none. (Two apparent counterexamples were
construction bugs — invalid Jacobi symbol on even moduli, a mis-set gamma parity — caught by verify-before-ship
and fixed; the clean result holds.)

**CORRECTION (LEDGER §969, HANDLE) — it is pole LOCATION-count, not pole ORDER.** The earlier "κ = ord_{s=1}(pole)"
was an overclaim from testing only orders {0,1}. Direct test of ζ^k (pole of order k): since −(ζ^k)'/ζ^k = k·ΣΛ and
ξ^k has Γ-factor^k, the form arch−prime for ζ^k = k·(arch−prime for ζ) — a positive scaling, inertia-invariant —
so **κ(ζ)=κ(ζ²)=κ(ζ³)=1**, not 1,2,3. The clean invariant is **κ(arch−prime) = #pole LOCATIONS at s=1 ∈ {0,1}**
(a pole only ever sits at s=1; higher multiplicity from ζ-powers just scales the whole form). The one negative
square = the single pole *location*, robustly, not its order.

Honest refinement of §8's "κ=1": it is also a BOTTOM-EDGE property (ζ has κ=0 in the deep interior [60,90]+),
so the pole contributes exactly one *edge* negative square. Robustness re-confirmed (LEDGER §964–967): κ stays
pinned at 1 as bandwidth grows (the growth flows into a near-null plunge; #pos saturates), arch alone is PSD so
the negative is prime-manufactured, and κ=1 is universal at saturation across the critical line (resolution
threshold grows ~linearly with height T). Strengthens the "one negative square = the pole" identification.

**Correct-code note (LEDGER §973):** the cross-Selberg test must use PROPER primitive Kronecker characters —
`operators/kappa_equals_pole_order.py` used `jacobi(n,q)`, INVALID for even modulus (q=4,8), which spuriously
gave χ_8→κ=1 (caught by a workflow verifier). The corrected `operators/kappa_pole_location_FIXED.py`
(χ_4=(-1/n), χ_8=(2/n) even, χ_8'=(-2/n) odd, (·/q) for odd q) gives the clean result: ζ→κ=1, all seven entire
χ (mod 3,4,5,7,8×2,11)→κ=0 at 4W=5,6.

**Explicit Krein–Langer factorization (LEDGER §973, machine-precision verified across a 4-lane workflow +
independent adversarial re-verification):** M_σ = arch−prime = K₊ − K₋ with rank(K₋)=1, both K₊,K₋ PSD,
reconstruction ‖M_σ−(K₊−K₋)‖_F = 1.6e-14; Pontryagin signature (79,0,1); and the POLE direction is **98.45%
aligned** (cosine 0.9845) with the single negative direction — corroborated independently by the exact-form
margin probe (pole-lift ≈ |μ|+0.03, LEDGER §972c) and the RMT lane (pole-norm/|λ_min| = 1.052±0.041). The pole
is near-exactly tuned to absorb the one negative square; κ(Q=arch+pole−prime)=0 to ±1e-13.

**Open bridge (LEDGER §974):** the m-function / de Branges pole-count route (argument principle, orthogonal to
Weil positivity) is NOT yet built — the naive E(z)=W(θ,z) from the screw string has NO UHP zero (winding 0),
so it does not realize κ=1 as a pole-shadow zero. The correct canonical system (whose m-function inherits the
explicit-formula pole as its single negative square) is the real unbuilt work. (Round-2 confirmed: building
m(z)=vp/vm and counting κ by UHP SAMPLING gives κ_sampled=10–13 ≠ κ_class=1 — a category error, the finite-sample
Nevanlinna-Gram count misses the indefinite-metric structure; vm's pole was not located. Route stays open.)

### 8c. Geometric containment: RH ⟺ arch−prime is PSD on pole⊥ (LEDGER §975–977, HANDLE)

The single negative square lives in an EXPLICIT 1-D space. As the window → the full critical line, the pole's
principal direction converges *exactly* onto M_σ's single negative eigenvector (cos → 1, lift/|μ| → 1; §975,
resolution-stable). Consequently the compression of M_σ=arch−prime to the orthogonal complement of the pole
direction is **PSD** (κ=0, machine precision, all a,K; §976) — the one negative square is contained entirely in
`span(pole)`. This is non-trivial: a *random* direction's complement leaves the negative at full depth (control).
And it is **faithful to RH**: perturbing the primes off the arithmetic point (δ=0.002) makes κ explode to ~22 AND
the negatives escape pole⊥ (§977). So:

> **RH  ⟺  κ(arch−prime) ≤ 1  ⟺  arch−prime ⪰ 0 on pole⊥  ⟺  neg-space(arch−prime) ⊆ span(pole).**

Via the explicit formula this is Weil positivity restricted to test functions vanishing at the pole — a concrete
finite-codimension PSD statement about an explicit arithmetic operator, distinct from (and sharper than) the §891
alignment margin. Non-circular linear algebra alone (S-lemma/interlacing) bounds #neg(M_σ) ≤ rank(pole) = 2; the
containment sharpens this to ≤ 1 by using *which* 1-D subspace of the rank-2 pole the negative occupies. Proving
the containment for all windows unconditionally is the (reformulated) RH.

### 8d. Universality across degree (LEDGER §978–979, HANDLE)

κ = #pole-locations at s=1 is verified across **23 degree-1 L-functions** (ζ, ζ², ζ³ → κ=1; 19–20 primitive
Dirichlet χ, moduli 3–61 → κ=0; zero exceptions) AND now a **degree-2 GL(2)** L-function: the Ramanujan Δ
(weight-12 level-1 newform, Hecke eigenvalues λ(p)=τ(p)/p^{11/2}, Γ_ℂ(s+11/2)) is entire and gives **κ=0, PSD
(min +0.29)**, robust across bandwidth — once the archimedean density is scaled by the **degree** (a degree-d L
has zero-density ∝ d; scale=1 gives spurious negatives). Two construction pitfalls, both fixed and documented:
the float η-product corrupts τ(n) past n≈113 (intermediate coefficients exceed 2⁵³ — use exact integers), and the
degree-scaling of the arch factor. The reframing is thus not ζ-specific — it is a structural feature of the
Selberg class across degree.
