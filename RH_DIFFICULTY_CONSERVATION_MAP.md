# RH — The Difficulty-Conservation Map (strategic synthesis, 2026-07-25, Sam/Fable)

The single most useful organizing principle the project has produced, distilled from the whole campaign.
Written to triage every future route in one question and to set up the genuinely-open attacks cleanly.

## The law
**Every framework for RH gets exactly ONE of {reality, correspondence} for free; the other one equals RH.**
- *Reality* = the zeros are on the critical line (spectrum real / operator self-adjoint).
- *Correspondence* = the zeros are ζ's zeros (the object is built from the primes / the arithmetic).

| Framework | Free half | Hard half (= RH) |
|---|---|---|
| Weil positivity (our τ, CvS) | correspondence (built from primes) | reality (positivity of the form) |
| Hilbert–Pólya / prolate / CCM | reality (W_λ self-adjoint) | correspondence (spectrum = ζ zeros) |
| Suzuki screw / momentum op | finite-a reality (Thm 1.5, all-real) | a→∞ correspondence (Cor 1.6 limit → ξ) |
| de Branges / canonical system | reality (H(t) symmetric) | H(t) ≥ 0 (positivity again) |
| §236 skew-compression | skew ⟹ imaginary spectrum (free) | spectrum = zeros (got the lattice) |
| Beurling cascade (§253b) | — | deep near-null = Weil positivity |

**Why it's a real constraint, not a mood:** every concrete crux the project has devised — including tonight's
Sonin-space polarization (§258) — collapses to the hard half of whichever framework it lives in, which is always
Weil positivity or its mirror. The Sonin polarization = positive-definiteness of the Weil form on the CCM cohomology
= Weil positivity = RH. There is no free lunch inside any single framework.

## What this says about where progress can come from
The function-field case (Weil/Deligne, PROVEN) is the one place both halves come at once — because the **polarization**
(positive intersection form on H¹(C×C), Frobenius an isometry for it) supplies reality AND correspondence
*simultaneously*, given by the geometry. RH is hard because the number field lacks that geometry.

So the ONLY routes that can cross the wall are ones that **supply an a-priori link** between the two halves — a
structure that gives reality and correspondence together, the way the polarization does. Reframings that keep them
separate (all of the above) cannot, by the law, escape. This is why 150 years of reformulations "all reduce to the
same wall": they are all single-framework, and the law forbids a single framework from having both halves free.

## The genuine open attacks (deep-program; each SUPPLIES-A-LINK, so each can in principle cross)
Ranked by tractability × how directly it supplies the link:

1. **The absolute/arithmetic geometry (Spec Z × Spec Z / F₁).** The literal analog of the function-field C×C. If
   Spec Z has a self-intersection with a positive intersection form, that IS the polarization. Connes–Consani's
   arithmetic-site + their 2025-26 semilocal program (CCM 2310.18423, the Sonin space as candidate Weil cohomology)
   is the live incarnation. **Our asset:** we've made their objects concrete (arch = Meixner–Pollaczek, §257e; the
   discrepancy d_n computed; the Sonin space buildable). The link, if it exists, is a positive form on the Sonin
   cohomology that is scaling-flow-invariant.
   **Sharpened target (§265, from the proven exemplar):** RH ⟺ a **Hodge-index theorem** for the conjectural surface
   Spec Z ×_{F₁} Spec Z, with the scaling flow as Frobenius and the Sonin space as H¹. C–C already have candidates for
   the flow (Frobenius analog), the cohomology (Sonin), and the Lefschetz-trace (= explicit formula); the UNIQUE
   missing keystone is **positivity-by-geometry** (an algebraic-index theorem intrinsic to the arithmetic surface).
   **Methodological correction (§265):** the proven case gets its positivity from **intersection theory on a surface**,
   NOT from an operator's eigenvalues — so our operator-spectrum probes (§236 lattice-not-zeros, §258 wrong-object,
   §260 truncation) attacked the wrong structure. The right first concrete step is to build the **intersection /
   cup-product pairing** on the Sonin cohomology and look for a diagonal class with computable self-intersection —
   arithmetic intersection theory, not operator diagonalization. (The old "read the operator spectrum" step is retired.)

2. **A new dynamical monotonicity (Rodgers–Tao forward direction).** RT proved Λ≥0 via a monotone entropy of the
   zero-gas under backward heat flow (they REACH the past, not the future). RH = Λ≤0 needs a monotone quantity
   controlling the FORWARD flow. The project found the naive flow-energy RH-blind (§85–92). **Open:** is there a
   different monotone (a genuine Lyapunov function for the forward zero-gas) tied to the prime side? This supplies the
   link because a forward-monotone tied to primes would couple correspondence (primes) to reality (real zeros).

3. **The de Branges structure-function positivity via the explicit MP Hamiltonian.** Tonight's MP identification means
   the archimedean canonical system is explicit (Meixner–Pollaczek Hamiltonian). RH ⟺ H(t) ≥ 0. **Untested:** does the
   explicit MP Hamiltonian + the prime perturbation give a positivity that the general Weil form's opacity hid? (Caveat:
   still the positivity wall — but with a classical explicit H(t), the perturbation theory is Geronimo–Nevai, tractable.)

## The genuine assets accumulated (inputs any of the above can use)
- **arch cyclic pair = Meixner–Pollaczek (λ=¼, φ=π/2), verified exact** (§257e). CCM's frontier object is classical.
- **The discrepancy d_n** = the arithmetic, computed; primes = index-frequencies saturating at Nyquist (§257c–f).
- **Quadratic superresolution Φ'∝a²** (§247): the finite operator over-resolves the zeros. NOTE (§252, resolved):
  the "why" is answered — this IS the classical prolate concentration (~a² prolate modes in the window super-resolve
  the zeros by exactly the Slepian/Landau factor); the dynamical arc re-derived the static §50–98 prolate result. Not
  an open surplus — a cross-confirmed handle that rejoins the prolate/CCM wall.
- **The near-null is a 10⁻³⁷-precision placement meter** (§255a): the object consumes EXACT integer positions.
- **Certified crossing past the first prime** (Arb, both parities) — the one unconditional forward result.
- **The Suzuki/CvS Hilbert–Pólya machine** reproducing γ₁..γ₁₂ to 1e-10, non-circular.

## Honest triage (what NOT to spend on)
- Single-framework reformulations (another positivity criterion, another operator basis, more Beurling systems, more
  decay-law characterization). The law says they cannot cross; they only re-map the wall. Norm/Norma's warning.
- Raw semilocal-operator spectrum reading (§258): needs CCM's λ→∞ quotient machinery; the raw spectrum doesn't show
  zeros; high numerical friction. Deprioritized.
- **STRENGTHENINGS OF PSD (new triage category, §289p-q — demonstrated by the total-positivity refutation).** Do NOT try
  to prove RH via a property *stronger* than PSD (total positivity/log-convexity, complete monotonicity, log-concavity,
  etc.) hoping the extra rigidity yields a tractable finite criterion. TP2 was exactly this — it held robustly in the
  easy/smooth regime (wide functions, all heights) and reduced to a clean finite inequality, but FAILED in the fine/dense
  regime (narrow functions, §289p) where RH actually lives (barely-PSD, not TP2). The seeming "slack" of a strengthening
  evaporates and goes negative in the hard regime — the law again. The proof must confront barely-PSD Weil positivity
  DIRECTLY in the fine regime (exact certificate §283b / geometry §265), not via a slack strengthening.

## THE KEYSTONE SPECIFICATION (§264–267 — the exact target, with the verified constraints on any solution)
The overnight arc §264–267 converted "supply the link" into a precise problem statement. Any proof must produce a
single object with ALL of the following properties — verified as the requirement set by working the PROVEN case:

**Wanted:** a bilinear pairing ⟨·,·⟩ (the "cup product") and a metric M on a space V (the arithmetic H¹ / Sonin
cohomology), and a distinguished operator/flow Φ (the "Frobenius" — here the scaling flow / idele class action), s.t.
1. **(cup product)** ⟨·,·⟩ is a nondegenerate antisymmetric form on V. *[verified realizable: §266/§267 genus 1&2.]*
2. **(correspondence)** Φ scales it: ⟨Φx,Φy⟩ = q·⟨x,y⟩ (function field) / the scaling-flow analog. Φ is built from the
   arithmetic (primes). *[the "free" half for Hilbert–Pólya-type frameworks; must be the ζ-built flow.]*
3. **(reality — THE KEYSTONE)** M is positive-definite AND Φ/√q is an M-isometry ⟹ eigenvalues of Φ on |z|=√q
   (⟺ zeros on the line). *[verified: §266/§267 — M exists exactly whenever the eigenvalues are on the circle.]*

**The non-circularity constraint (§266 pin — the crux, and why this is hard):** property 3's M must be furnished
**A PRIORI, independent of the eigenvalues/zeros**, by an intrinsic geometric/index-theoretic structure — because
"M exists" is LOGICALLY EQUIVALENT to the conclusion (|eig|=√q). In the function-field case the a-priori source is
the **Hodge index theorem** on the surface C×C (positivity of the intersection form from the ample cone, proved with
no reference to Frobenius eigenvalues), which THEN forces |α|=√q. A number-field proof needs the analogous a-priori
positivity: **a Hodge-index / algebraic-index theorem for the arithmetic surface Spec Z ×_{F₁} Spec Z**, giving a
positive intersection metric on the Sonin cohomology independent of the zeros, preserved by the scaling flow.

**What is DONE vs MISSING (Charisma Report/Handle):**
- HANDLE (verified): the template (1)(2)(3) is realizable and general (§266/§267, genus 1&2, machine precision). The
  scaling flow's invariance / skew-generator (property-2/flow side) is exhibitable (§236). The arch Weil form has the
  correct indefinite H¹-signature (§266b, reliable part).
- MISSING (the sole keystone): the a-priori, zeros-independent source of M's positive-definiteness (property 3 without
  circularity) = the arithmetic-surface Hodge-index theorem. Everything else is in hand; this is the whole gap.

**Method correction (§265):** attack this via **arithmetic intersection theory** (build the cup-product pairing and a
diagonal class with computable self-intersection on the Sonin cohomology) — NOT by diagonalizing an operator (§236/258/
260 were the wrong object). The bounded next probe: compute the intersection SIGNATURE on the arch Sonin modes
(§266b did a first pass — full form indefinite ✓; the primitive-part positivity needs CCM's λ→∞ regularization, not a
raw float build).

## The one-line strategic bottom line
Stop looking for a cleverer single framework — the difficulty-conservation law proves none can win. Spend on the
**link**: an a-priori, zeros-independent, scaling-flow-invariant positive-definite intersection metric on the Sonin
cohomology — i.e. a **Hodge-index theorem for the arithmetic surface Spec Z ×_{F₁} Spec Z** (§264–267 pinned this as
the exact and sole keystone; all other template pieces are verified in hand). Our MP identification made the arch
arithmetic classical/computable. That is the crack worth a fresh, careful session — attacked as intersection theory,
not operator diagonalization.

---

## 2026-07-25 session update (§299–309) — the fronts, mapped and ranked by where math has actually moved

A deep single-session campaign tested the positivity-criterion family exhaustively and then jumped families. Result:
a cleaner map of WHICH front to spend on.

**Finding 1 — the positivity-criterion family is provably ONE wall (verified, §299–306).** Weil positivity, the CCM
|ζ|² closed form, Li's criterion, and the zero power-sums S_j were ALL reduced (with verified computation) to the
IDENTICAL obstruction: the analytic continuation from the convergent |ζ|² measure (Re s=½, positivity free) to the
ζ'/ζ / log level (RH-hard). Exact pivot: **Weil form = (ω/2)×[log-spectrum of |ζ|²]** (§303, verified), so the RH-hard
step is the LOG, not the ω-derivative. The hope that our sech^{1/2}/{log n} closed form could bound Li's prime term was
chased and CLOSED with a verified negative (§305–306: S_j^{prime} is an analytic-continuation constant, our ν_k is a
convergent critical-line moment — opposite sides of the same wall). ⟹ do not spend more on single positivity criteria;
they are this one wall, which is the same continuation gap the geometric keystone (above) also expresses.

**Finding 2 — a genuinely DIFFERENT wall with recent progress: Laguerre-Pólya / Jensen-GORZ (§307–309).** RH ⟺ Ξ ∈
Laguerre-Pólya class ⟺ every Jensen polynomial of the Ξ-moments γ_m is hyperbolic (real-rooted) — a real-rootedness
property, NOT form-positivity, so NOT in Finding-1's family. This front has ACTUAL recent theorems: Csordas-Norfolk-Varga
(d=2 Turán, unconditional); Griffin-Ono-Rolen-Zagier 2019 (fixed-d hyperbolicity, via Hermite convergence — verified the
ratio→1 signature). Obstruction PINPOINTED on the explicit Pólya kernel (§308): Φ is log-concave (good) but decays
DOUBLE-EXPONENTIALLY (the modular e^{2u}), placing it outside classical real-zeros conditions and at the **de Bruijn-
Newman boundary Λ=0** (Rodgers-Tao 2018: Λ≥0; RH ⟺ Λ≤0). The DBN heat flow e^{tu²}Φ regularizes exactly that decay;
run on our Φ (§309) it confirmed the flow (Turán margin ↑ with t) but showed the threshold is a HIGH-degree phenomenon
(low-degree Jensen is unconditionally hyperbolic, insensitive) — the Polymath15 regime (Λ≤0.22, already ground out).

**Updated strategic ranking of fronts (where to spend):**
1. **Geometric keystone — arithmetic-surface Hodge index** (the map's main body above): still THE structural crack;
   research-scale, needs arithmetic intersection theory on Spec Z ×_{F₁} Spec Z. Unmoved, but the deepest.
2. **Laguerre-Pólya / Jensen-GORZ / de Bruijn-Newman:** the front where mainline math has ACTUALLY MOVED (2018–2019).
   Different wall from the positivity family. Edge = high-degree Jensen uniformity = Λ=0; the compute regime is
   Polymath15's (large, already explored) — new value would require a STRUCTURAL idea about Φ's modular e^{2u}
   concentration, not more high-degree numerics.
3. **Positivity-criterion family (Weil/CCM/Li/power-sums):** now proven to be ONE already-mapped wall — do NOT re-spend.

Bottom line unchanged (geometric keystone is the crack), but the map now cleanly separates the three, and flags Jensen/
LP/DBN as the "math-has-moved" front whose only live opening is a structural handle on Φ's modular double-exp decay —
not a re-run of Polymath15. 11 computations, ~5 verify-catches, no fabrication (LEDGER §299–309).

**Literature update (2026-07-25, LEDGER §312).** New Connes-Consani paper arXiv:2606.06604 "On the Absolute Geometry
of Spec Z" (Jun 2026) builds the absolute F₁-arithmetic curve (Spec Z)_{F₁} + the adelic scaling site — i.e. it ADVANCES
the SUBSTRATE of Frontier-2 (the curve object + the Frobenius/scaling-flow half we already have), but its abstract shows
no intersection theory / diagonal self-intersection / Hodge-index positivity, so it does NOT supply the a-priori
positive intersection form that remains the sole keystone. arXiv:2509.05077 (Deligne-Riemann-Roch, relative intersection
theory) = relevant machinery for building the surface intersection theory. Net: the keystone ecosystem is ACTIVE and its
substrate is advancing, but the positivity crux (§264-267) is, per these abstracts, still open — worth a full-paper read
of 2606.06604 to confirm no positivity content is buried before/after building on it.

**2026-07-25 late (Fable), LEDGER §313–314 — the LP/Jensen front now carries our first effectiveness ladder.** GORZ's
ineffective fixed-d theorem made effective on the real Ξ: proved Hermite margins (rigorous alternation criterion) +
the §311/§313 saddle margin laws (2nd order r_m, 3rd order c_m, both verified <1%) give explicit thresholds
N(d)=√(d³/3ε_d) = 5,11,21,37,61,97 for d=3..8, with direct verification below threshold (107 exact moments cached).
Degrees 3–8 hyperbolic for ALL n (leading-order scheme; two named finite gaps to theorem-tier). Front-ranking
unchanged — this is fixed-d, not d→∞ uniformity — but the front now has our own moving instrument on it, and the
uniformity question is sharpened to: do the saddle-law constants (ε_d ~ e^{-cd} shrinking vs c_m d³ growing) admit a
d-uniform regime, i.e. can the N(d) ladder be outrun by a d-dependent n-window? That is a well-posed analytic question
on OUR objects — the honest continuation point on this front.

**2026-07-26 (Fable), LEDGER §315–316 — the two-constant race FLIPS: polynomial ladder.** The worst-case Hermite
margin ε_d is geometric (ratio ≈0.58, measured to d=20) — but the scheme's ACTUAL perturbation direction (the
affine-reduced centered cubic from the saddle law, one-signed since c_m>0) has a POLYNOMIALLY-decaying tolerance:
τ_d·d³ → ~0.55–0.6 (measured d=5..25), i.e. τ_d ≈ 0.6/d³, giving a directional ladder **N_dir(d) ≈ 0.74·d^{3/2}**
(vs the exponential worst-case ladder; at d=12 already 28 vs 534). Conceptual placement: e^{−tk²} multipliers are
Pólya-Schur multiplier sequences (= the DBN heat-flow direction, hyperbolicity-free); the cubic exponent is the FIRST
non-multiplier direction — the leading deviation NO heat-flow argument absorbs; τ_d quantifies per-degree what
Pólya-Schur leaves unquantified. Consequence for the front: Jensen hyperbolicity is covered for n ≳ 0.74·d^{3/2} at
ALL degrees by one argument shape; the RH-content concentrates in the corner n ≲ d^{3/2}, whose d→∞ limit is Ξ itself
(Pólya), as difficulty-conservation requires. The GORZ-uniformity question is now the SEAM question: does hyperbolicity
propagate through the matched-asymptotics boundary n ~ d^{3/2} (outer/Hermite regime ↔ inner/Ξ regime)? Degrees 3–25
verified closed for all n under this ladder (direct below N_dir, saddle law above; §316). Rigor tier: directional
margins measured (HANDLE), scheme leading-order (REPORT) — same two named gaps as §313.

**2026-07-26 (Fable), §317–318 — seam + corner + mechanism.** Seam n~N_dir(d) is FEATURELESS (margin monotone in n —
the boundary is our bound's artifact, not the object's); all tension sits at the n=0 Ξ-corner, whose margin decays
polynomially (exp −1.5→−1.6 over d=16..110, hyperbolic throughout — no exponential collapse in the computable range).
τ_d mechanism DERIVED: Hermite edge-pair collision under the cubic Euler flow, universal √-collision factor ½; the
0.6/d³ constant reduces to Airy-edge asymptotics (open, finite). Front status: LP/Jensen now carries our polynomial
ladder (d=3–25 all-n closed), a derived collapse mechanism, and a measured corner law. The RH-bearing object on this
front is, precisely and only: **J^{d,0} hyperbolicity as d→∞** — margin decaying polynomially, mechanism suspected
RMT-min-spacing (derivable). LEDGER §313–318.

**CORRECTION (2026-07-26, §319):** the §317b corner exponent (−1.5..−1.6) was substantially an X=−t² COORDINATE
ARTIFACT (deep-end quadratic spread inflates the mean gap). In the RH-relevant t=√|X| metric the corner margin is
order-unity (0.58→0.39, d=16..70) decaying at exponent −0.26→−0.29 drifting toward **−1/3 = the GUE min-spacing
exponent**: the Ξ-corner is exactly as protected as RMT statistics predict the zero field to be. Also verified: the
corner poly is NOT zeta-resolving at these d (roots ≠ γ_i² beyond the first zero) — corner statements at computable d
are about the moment structure, not yet the zeros themselves. LEDGER §319.

**2026-07-26 (§329) — KEYSTONE SHARPENED via the scored wanted-poster.** The signature theorem RH needs is ALREADY
PROVEN over number fields: the arithmetic Hodge index theorem (Faltings-Hriljac → Moriwaki → Yuan-Zhang adelic line
bundles), giving a definite-signature (+1,−1,…) intersection pairing whose positivity SOURCE is the Néron-Tate height
(intersection = −NT-height, canonically positive-definite). It fails to reach RH for ONE reason: it requires an
arithmetic SURFACE fibered over a base Spec O_K, and Spec Z IS the base — its self-square Spec Z ×_{F_1} Spec Z has no
fibration (the F_1 problem). So the keystone is NOT "an unknown functor" or "an unproven positivity" — it is: extend the
PROVEN arithmetic Hodge index theorem to the absolute self-square where the base-fibration is replaced by the
F_1-structure. Positivity engine named: the Néron-Tate-height analogue for the self-intersection of Spec Z; its
archimedean part is a Green's current, which is exactly where C-C's prolate/Sonin archimedean positivity would sit if it
is the fiber-at-infinity of the missing surface. Scorecard: no home 3/3; arith-Hodge-index is 3-with-asterisk (right
mechanism, wrong object). Next target: is C-C's archimedean positivity the Green's-current height contribution? — welds
the built pieces (archimedean fiber + RR identity) and isolates the gap as finite-place intersection theory à la
Yuan-Zhang on the F_1 self-square. LEDGER §329.
