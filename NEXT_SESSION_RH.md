# NEXT SESSION — updated 2026-07-28 (Riemann7 long autonomous stretch: LEDGER through §378)

**Grounding doc for next table: `ROUNDTABLE_GROUNDING_2026-07-27.md`** — current through §375, 6 addenda.

THE ARC SINCE §360 (all HANDLE unless noted):
- **Density law (§366-367):** resonant-core count R ≈ ρ(T)/w² (Weyl zero density / w², clean rel-rms 0.11);
  the PROVEN elementary-certificate reach tracks Nyquist 1/ρ(T) to 5-11% across T=29..440. Certificate geometry
  density-locked on all 3 axes.
- **s↔γ (§368):** critical resonant scale ~T/2π, RS-order but O(1)-drifts (not pinned).
- **Q3 lit-gate (§370):** arch matrix = named Gaussian localization / cross-Toeplitz operator; collapse =
  BGLy-2009 / Seip critical-density frame theorem (rate (1−a²)²); UNIFORM sub-critical bound = OPEN (=Q1).
- **§372→372b (self-correction):** measured Weil≈0.5·overlap in the sparse interior, briefly claimed Q1 reduces
  to the pure frame bound — REFUTED at criticality (frame stays healthy lmin≈3 while Weil cliffs to 0 + bent-prime
  flips NEG). Arithmetic IS the corner, frame blind to it.
- **THE CORNER TRIO (§373-375):** collapse pins to EXACT RvM count (§373, n*/RvM→0.99); NO fixed-direction
  zero-free certificate (§374, tight direction is arithmetic/circular); that direction IS the zero-interpolating
  mode (§375, φ vanishes on zeros, BRS §336 realized). Corner mapped to its core, reaffirms wall (Q5).
- **Universality (§376-378):** built+validated a SECOND instrument (Fejer/sinc², compact-freq, finite prime
  side, l2_fejer.py, rel resid 4e-7). RvM-collapse + zero-interpolation corner + δ² off-line floor-lift all
  REPRODUCE — the phenomena are intrinsic to the explicit formula, not Gaussian artifacts.
- **k=3 (§356-361b):** pairwise-exhausted at log order, M-S-2004-anchored; exclusion program (§364-365) parked
  (mechanism real but no primes-only certificate, §374 obstruction).

STATE: RH-on-this-front = positivity of the Weil form along the zero-interpolation mode at critical density —
arithmetic, conservative (Q5). The concrete deliverable is the elementary theorem grid (§359-362): zero-free
interval-arithmetic PD certificates, 4 bands × 3 widths, density-governed boundary.

OPEN LEVERS (ranked): (1) the M-S pair-rigidity → corner-inequality connection (§369, the one ARITHMETIC angle
on the corner not yet constructed — sharpest); (2) extremal Beurling-Selberg test functions to optimize the
certificate CONSTANT (won't remove the wall per §377, but sharpens reach); (3) the uniform sub-critical frame
bound (§370, genuinely open frame-theory frontier, our per-config certs are a constructive foothold).

Standing: 15 rh-working-rules; ask-don't-assert; Report/Handle; LEDGER append-only. Verify-catch classes this
arc: 3× log-vs-constant-order anchoring, arb ball-equality/non-dyadic/ball-digamma, §372 over-claim (caught same
session by pushing to criticality). Tools: rt5_cT/cT2/sgamma/bgly/uniform/ratio/cliff/nullproj/nullshape/
interp_offline; l2_fejer (2nd instrument); arb_elementary_cert (--w/--h/--band), arb_hybrid_cert (--w/--n0).

---
# NEXT SESSION — updated 2026-07-28 early (Riemann7 autonomous arc complete: LEDGER §343-360c)

**Grounding doc (paste into next table): `ROUNDTABLE_GROUNDING_2026-07-27.md`** — current through §360, 3 addenda.

THE THEOREM MAP (all interval-arithmetic proven, zero-free):
  collective (integral-grade): w=3: s=3.5,3.4 (§340) · w=2: s=3.9,3.8 (§350) · w=1.5: s=4.1,4.0 (§360)
  ELEMENTARY (zero numeric integration): w=3: s=4.5→4.0 · w=2: s=4.5→4.3 · w=1.5: [run landing at close]
  Arithmetic inputs: a few hundred prime powers + in-script-proved ε + RS 1962. Plus the new 2-line
  monotonicity lemma (Ω↑ on (0,∞) via termwise ψ′-series sign).

The autonomous arc also settled: Q4 (Bragg law on the projection object, 4 bands; fixed-point excluded;
involution → the last enemy's dual is the PNT MAIN TERM §358); Q6 k=3 (pairwise-exhausted at log order,
M-S-2004-anchored to 4 decimals §356-357); saints' k=2 arrow (all three pre-registered tests ARROW §352,
intercept matches 1−γ−log2π); Q2 partial (R = C(T)/w², C band-dependent; slot factorization failed §354).

Open fronts, ranked: (1) grounding Q1 — the w→0 statement of the hybrid; (2) k=3 NEXT-order probe vs the
M-S O(H^{3/2}) remainder (pre-register first!); (3) derive C(T); (4) h-refinement ladder if deeper elementary
rungs wanted (available, not strategic); (5) elementary chain at higher bands (T-flat per §348, should be cheap).

Standing: 15 rh-working-rules; ask-don't-assert; Report/Handle; LEDGER append-only; 4 verify-catch classes
added this arc (log-vs-constant order anchoring ×3, arb ball-equality/non-dyadic/ball-digamma).

---
# NEXT SESSION — updated 2026-07-27 (Riemann7 close: the collective arc, LEDGER §343-350)

**Grounding doc for the next table: `ROUNDTABLE_GROUNDING_2026-07-27.md`** (supersedes 07-26; paste that file).

State in one paragraph: RT#5 ran (Q1 spine/Q5 guard), its five registered specs were executed, and the arc
they opened ended in a NEW ARCHITECTURE: the hybrid/collective chain — R ≈ 40/w² resonant primes exact
(quadratic count law, deconfounded §346) + k=5 directional bound on remaining primes < n0 ~ e^{10.5/w} (bulk)
+ PNT-smooth tail integral + fluctuation charge ε‖C‖ (partial summation nearly tight, §347; Chebyshev-grade
ε suffices, §348) = FULL REACH to the Weyl rank limit at every width tested, T-axis FLAT (§348), and
THEOREM-GRADE at w=2 s=3.9 via interval arithmetic with ε ≤ 0.0325 proved exhaustively (§350). Exponential
bulk / polynomial structure. The wall now = w→0 uniformity of the polynomial structure part (grounding Q1).

Next work moves (any session):
1. Grounding Q2 — DERIVE the 40/w² law (envelope × Bragg-slot density × band factor); test its T-prediction.
2. s=3.8 arb width fix (prec-320 retry was in flight at close; check arb_hybrid_38.out; if still wide, split
   the arch integral or raise degree) + arb at w=1.5 s=4.0 (n0=3000 — bigger ψ check, same machinery).
3. Named lit constants: swap the in-script proved ε for RS explicit bounds (headroom 18-1693x, bookkeeping).
4. Grounding Q4 — the s↔γ calibration for the involution branch.
5. §341 minorant-ification of the arch entries INSIDE the hybrid (makes the whole B_lower formula-only).

Standing: rh-working-rules (15) unchanged; ask-don't-assert endings; Report/Handle labels; LEDGER append-only.

---
# NEXT SESSION (updated 2026-07-26 evening, LEDGER through §342) — READ THIS + ROUNDTABLE_GROUNDING_2026-07-26.md FIRST

## CURRENT STATE in one paragraph
The certificate front moved decisively today (§334-342): the dense-regime Weil-form collapse is identified and
measured as a CRITICAL-DENSITY SAMPLING phenomenon (rank saturation at the Riemann-von Mangoldt count; PSD
throughout; sign-blind controls show no threshold; bent primes only detectable at the edge). A rigorous
zero-free certificate chain (2-block gap bound + adaptive exact prime core) reaches the Weyl rank limit on
EVERY band and EVERY width tested, tracking true margins at 95-99.9%; Arb-rigorous PD certified at the plunge
floor (s=3.5, 3.4); a fully a-priori formula-only version certifies to s=4.2. The uniformity wall now appears
as a measured COST LAW: exact core grows ~e^{c/w} toward "all primes" as the test class completes. Lit-gated:
the synthesis appears unpublished; all components (BGL, Landau-Widom, BRS/Kulikov, Groskin) published.

## For the NEXT ROUNDTABLE: paste ROUNDTABLE_GROUNDING_2026-07-26.md (complete banked block + 5 open
## questions Q1-Q5, of which Q1 = beat the e^{c/w} core growth and Q2 = the frame-theory uniformity
## transplant are the sharpest).

## Next work moves (any session)
1. Q1/Q2 from the grounding doc (the substantive attacks on the measured wall).
2. Mechanical: close a-priori s=4.2 -> 3.9 ([0,7] piecewise minorant + smaller h, §341 loss budget).
3. Mechanical: Arb-rigorous chain (not just endpoint LDL) + the [50,80] s=1.8 row.
4. Ingredients: localized-weight Toeplitz floor; Bogoya-Boettcher-Grudsky lambda_4(A) (§340b flags).
Tools: operators/{l2_edge_sweep,directional_cert,resonant_cert,band_sweep_cert,kr_frontier,arb_dense_cert,
apriori_chain,arch_floor,w_sweep_chain}.py — all validated; arsenal/ gates (t01 prefilter before any new route).
Rules: rh-working-rules (all 15). Verify-catches today: 2 (sieve bug via cos(c·log6) residual; arb x**2
non-dyadic-ball NaN). Charisma labels used throughout: §335-342 measurements HANDLE w/ proof pointers,
identifications REPORT.

---
_(older handoffs below — SUPERSEDED, retained for provenance)_

# NEXT SESSION — read this block FIRST (updated 2026-07-26, LEDGER through §333b)
# (Everything below the "=== ARCHIVE ===" line is older context, kept for provenance.)

## WHERE WE ARE (the §325–333 arc — the sharpest the keystone has ever been)

RH has been localized to **ONE homeless clause at ONE place**: the **signature** (the inequality), at the
**archimedean place ∞**. Identity and per-place estimate are BUILT; only the a-priori uniform *sign* is missing.

**The scorecard frame (Roundtable #2, LEDGER §327):** criteria certify, objects PROVE. The explicit formula IS the
Lefschetz trace identity (banked since the 1950s) — over F_q, RH lived entirely in the SECOND theorem (Hodge-index),
not the trace. A home must be (1) PRIOR (structure before the zeros), (2) GENERATIVE (rebuilds the primes),
(3) SIGNATURE-BEARING (the inequality, not the identity). NO home scores 3/3.

**The decisive results:**
- §328 (check-8): Connes-Consani positivity is NOT signature-bearing — it's a spectral ESTIMATE (prolate/Toeplitz,
  all-primes uniformity open) OR an Euler-characteristic IDENTITY (Riemann-Roch/Jacobian, no polarization). Verified
  by reading their 2020-2026 papers.
- §329: the SIGNATURE theorem RH needs is ALREADY PROVEN over number fields — the **arithmetic Hodge index theorem**
  (Faltings-Hriljac → Yuan-Zhang adelic line bundles), signature (+1,−1,…) with positivity = −(Néron-Tate height).
  It fails RH ONLY because it needs a curve FIBERED over a base Spec O_K, and **Spec Z IS the base** — its self-square
  Spec Z ×_F1 Spec Z has no fibration (the F_1 problem). Keystone reframed: "extend a PROVEN theorem past its fibration
  hypothesis," not "build an unknown functor."
- §330-331: the archimedean fiber of that missing surface = C-C's prolate positivity = the arch Weil form's positive
  PRIMITIVE part after removing the low-r indefinite "fiber" band. Ω(r)=Re ψ(¼+ir/2)−log π is INDEFINITE (<0 for
  r<6.29). Computed signature (12,4); C-C's Sonin compression = Hodge-index's "restrict to primitive." Welds to our
  own §298/§299-302 (RH-difficulty localizes to the same low-r subspace). NOTE: this REPRODUCES established C-C-M
  Sonin/prolate structure (labeled honestly).
- Roundtable #3 (§332, +Grothendieck +Witten): the FIRST object to construct is the **normal bundle of the diagonal**
  (Ω¹(Spec Z)=0 in schemes; the arithmetic differentials δ_p in Borger's Λ-geometry supply it). The homeless
  coordinate = the **cross-prime pairing ⟨δ_p,δ_q⟩**. VERIFIED: δ_p (=(n−n^p)/p) don't commute (Grothendieck's hand
  values exact), commutator is a rich antisymmetric table, twisted-Leibniz/Hirota holds. Prediction: the signature
  class CONCENTRATES AT ∞ (initial Λ-structure ⟹ coboundary at finite places) = an APS eta / Weil-metaplectic-index of
  the scaling flow — CONVERGES with our §331 arch numerics. Sharper weld test: OS reflection-positivity of the prolate
  kernel (x→1/x) → self-adjoint generator by theorem.

**THE GENUINELY-NEW OBJECT (the one thing with no literature to rediscover):** the **OFF-DIAGONAL SIGNED cross-prime
pairing ⟨δ_p,δ_q⟩** and its concentration at ∞. Everything around it is built: Burnol (2004) + C-C-M (2024,
arXiv:2310.18423) own the Sonine/prolate/self-dual archimedean structure; Borger-Buium (0908.2512) own the multi-prime
differential structure incl. the DIAGONAL arithmetic Laplacian (Σδ_p²); the arithmetic Hodge index (Yuan-Zhang) owns
the signature-via-height mechanism. The off-diagonal SIGNED pairing between different primes' differentials is the
homeless part (Grothendieck's boundary, confirmed at abstract level; needs a full 0908.2512 read to be certain).

## DO NEXT (in order — READING before COMPUTE; the RH neighborhood is exhaustively worked, verify-before-build is
## non-negotiable — 4 rediscoveries caught this session: §325 LP, §327 control-grid, §328 C-C, §333 OS/Sonine)
1. **Read Borger-Buium 0908.2512 IN FULL** ("Differential forms on arithmetic jet spaces"): does the arithmetic-
   Laplacian / de-Rham structure equip the dp differentials with a SIGNED cross-prime (off-diagonal) pairing, or only
   the diagonal Laplacian Σδ_p²? If off-diagonal signed pairing is ABSENT ⟹ that is the object to construct.
2. **Read Burnol JTNB 2004 (Sonine spaces) + arXiv:2310.18423 (C-C-M 2024 semilocal prolate)** BEFORE any OS/reflection
   test — they own the x→1/x self-dual + prolate spectrum split. Extract what's established; then the OS-positivity→
   generator FRAMING (if new at all) sits on their space.
3. **Construct the cross-prime signature object correctly** (only after 1-2): a symmetric SIGNED pairing on {δ_p}, test
   whether it's a coboundary at finite places (Grothendieck predicts yes) and whether its class concentrates at ∞ as a
   Weil-index/APS-eta of the scaling flow. Do NOT guess it from the commutator values — that's the "beautiful dead
   computation" trap.
4. Owed lit-greps if pursuing the ∞ side: Weil 1964 metaplectic cocycle (is the ∞-class this?), APS 1975 eta / Callan-
   Harvey inflow, Segal-Wilson (cocycle-as-positivity), Borcea-Brändén / Lee-Yang (QFT home exemplar).

## KEY FILES: LEDGER.md §325-333b (full detail); RH_DIFFICULTY_CONSERVATION_MAP.md (keystone map, current); PROCESS.md
## (Roundtable #2/#3 arcs); memory [[rh-roundtable-literature-gate]] [[connes-decisive-diagnostic]] [[riemann-project]].
## TOOLS this arc: operators/{link_*, lp_*, dstar_experiment, joint_limit, borel_*, mixture_lp*, arch_signature}.py;
## moment caches operators/_gam{96,106,130}.pkl.

## DISCIPLINE (hard-won this session): the multi-Fable Roundtable (localhost:7796) is a proven tool — run it for
## adversarial pre-novelty review + strategy. But EVERY agent claim is Report-tier: verify the load-bearing computation
## yourself and literature-gate before booking novelty. The LP/Jensen family is BANKED-CLOSED (kernel-local, RH-free,
## d~T² published to ~10^20). Do not reopen it.

=== ARCHIVE (pre-§325 context below) ===

# NEXT SESSION (updated 2026-07-25 OVERNIGHT + MORNING STORMS, LEDGER through §283b, STORM #1/#2 complete — read RH_DIFFICULTY_CONSERVATION_MAP.md FIRST)

## STORM #1/#2 RESULT (2026-07-25, HONEST — corrected from a fork's overclaim; see LEDGER §296-298)
The forks did genuine work but a fork over-framed it as a "breakthrough / RH follows." Corrected reality:

**STORM #1 (certificate) — genuine findings (HANDLE), honest significance:**
- Extracted the saturating/extremal Weil function (min-eig ~1e-18 on an oversampled band). Coefficients [0.4409,
  -0.7882, 0.4103, -0.1246, 0.0221, ...], geometric decay q~0.154, super-exp spatial decay. Zero-vanishing verified (1e-15).
- Classical theta (θ₃,θ₄,E₄) RULED OUT (honest negatives).
- **RECOGNITION (§296):** the extremal's ZEROS are EXACTLY the zeta zeros (γ₁..γ₆ to 5 digits; spacings = zeta gaps). So
  the extremal = xi(1/2+ix)·envelope = the **de Branges Hermite-Biehler E-function** of the space H(E) that IS the Weil-
  distribution space (§294, arXiv:2301.00421). Fork's "BRS kernel" guess = this (BRS build interpolation from exactly E).
- **HONEST significance (NOT a breakthrough):** recognizing the extremal as the de Branges E-function is the KNOWN
  analytic wall (§294: Weil≈deBranges≈Nyman-Beurling). It does NOT give the certificate — proving E's Hermite-Biehler
  positivity (de Branges' H(t)>=0) IS RH. So this validates/pins the object (real, useful) but "RH follows" is FALSE.
  The certificate = proving de Branges positivity = the open problem, not a formula-lookup.

**STORM #2 (Castelnuovo/geometry) — closed to research-scale, honestly (§297-298):**
- A fork tested Castelnuovo-Severi using the PSD Weil form as the intersection pairing -> artifact (PSD can't be
  hyperbolic; §297 verify-catch). Corrected test with the arch form (§298): arch is POSITIVE in the zero-band (Omega>0
  for r>6.3), also not indefinite. => the Weil space does NOT realize the indefinite intersection form Castelnuovo needs.
- VERDICT (correct, right reasoning §298): the geometric/Castelnuovo route genuinely needs arithmetic SURFACE Riemann-Roch
  (Spec Z x Spec Z, Deninger/Connes-Consani, OPEN) -- NOT bypassable via finite Weil numerics. Research-scale keystone.

**NET of the storm (honest):** both routes converged on the KNOWN wall. #1 = de Branges positivity (analytic wall);
#2 = arithmetic intersection theory (geometric wall). Difficulty-conservation confirmed again: both halves are the wall,
the escape needs them supplied together (geometry). No breakthrough; a validated, sharply-pinned pair of walls + the
extremal object (de Branges E) concretely in hand. The forward work is genuinely research-scale (de Branges HB positivity
/ arithmetic Riemann-Roch), not a finite computation.

---

## TP2 / TOTAL-POSITIVITY LEAD (§289-289p) — REFUTED as a complete route; real only for WIDE functions (READ THIS)
STATUS: DOWNGRADED / REFUTED (§289p). The total-positivity/log-convexity route is NOT a complete path to RH. The Weil
kernel is log-convex (TP2) for WIDE test functions (w>=2) but log-convexity FAILS for NARROW ones (w~1.1-1.7, D_full<0,
§289p). Since TP2 is STRICTLY STRONGER than PSD and PSD (=RH) holds throughout, the kernel is PSD-but-not-TP2 for narrow
functions => TP2 cannot certify Weil positivity. The robust TP (§289-289n) was WIDE-functions-only; it does NOT extend to
the fine/narrow (dense-analog = RH-hard) regime. So the REAL leads remain the exact certificate (below) and the geometry;
the TP2 finding is a genuine phenomenon for smooth functions but NOT a route. (Kept for the record: it was rigorously
stress-tested and honestly refuted -- a clean example of verify-before-ship on a promising-looking lead.)
--- original (now-refuted) writeup below, kept for provenance ---
EMPIRICAL DISCOVERY (robust for WIDE functions ONLY -- see §289p refutation above): the full Weil form appears TOTALLY
POSITIVE (all minors >= 0) and TP reached into the dense regime (spacing 4-5) where elementary certificates failed --
BUT this was width-dependent (w>=2) and FAILS at narrow width. Confirmed general across bands/heights FOR WIDE functions.
WHY (partial): arch=Meixner-Pollaczek (§257e) => arch kernel is Cauchy-Stieltjes/Hankel => TOTALLY POSITIVE (classical);
empirically the prime perturbation PRESERVES TP -- likely because arch dominates prime (§276) enough to stay in the (open)
TP cone. TP => PSD, and TP is rigid + richly theorized (Karlin "Total Positivity", sign-regular / variation-diminishing kernels).
THE ROUTE (if it holds): prove the Weil kernel is TOTALLY POSITIVE (arch-MP-TP + prime perturbation stays sign-regular) =>
PSD => Weil positivity => RH. Reaches the DENSE core the certificate/geometry routes target, via a DIFFERENT (TP) mechanism.
REDUCED to a concrete inequality (§289f-h): the binding constraint is 2nd-order (log-convexity/TP2); prove the Weil-kernel
2x2 minor >= 0, i.e. [net prime minor-effect]/[arch log-convexity margin] < 1, uniformly, HARDEST at the first zeros (x~14).
Both sides explicit closed forms (arch margin from MP/Cauchy-Stieltjes log-convexity; prime effect from Sum Lambda(n)/sqrt n
[...] bounded via §280/§282 tail control). 3-step proof path: lower-bound arch margin, upper-bound prime effect, show ratio<1.
FINITE-INEQUALITY FORM (§289k-l, sharpest): the differential defect D=K*Kxy-Kx*Ky gives, at the first zeros,
D_arch=0.0921 > D_prime=0.0824 (ratio 0.894, matches the discrete h->0 limit -- cross-check PASS). D_prime is FINITE:
2,3,5,7-dominated, super-exp tail (§280 FT-decay), converged by prime 7. So the ENTIRE keystone = prove the finite
inequality D_arch > D_prime: (1) D_arch>0 lower bound from classical MP/Cauchy-Stieltjes log-convexity; (2) D_prime upper
bound = explicit ~4-term small-prime sum + Gaussian-FT tail bound; (3) show (1)>(2), tightest at first zeros (~0.1 slack;
higher zeros have larger D_arch => more slack). Every piece explicit/finite -- the proof of this finite inequality is the
research-scale crux, but the target is now a concrete finite 2nd-order inequality, not infinite-dim positivity.
HONEST MARGIN (§289i/j, CORRECTED -- verify-catch on an earlier coarse-scale over-estimate): the CONTINUOUS (h->0) ratio at
the first zeros is ~0.89-0.91, NOT the ~0.62 that coarse h=2 suggested. So the slack is only ~0.1 (TIGHT) -- but it does NOT
-> 1, so TP2 is a genuine sufficient condition WITH slack, not RH-restated-with-zero-margin. Higher zeros have MORE slack
(arch margin grows); the first zeros are the binding case (~0.90).
MAKE-OR-BREAK (research-scale but NAMED + CONCRETE): prove the minor-ratio < 1 uniformly (tight ~0.1 at first zeros). THE
THEORY: why does the prime kernel preserve log-convexity? (variation-diminishing / sign-regular composition; the arch
log-convexity margin dominating the oscillatory prime by ~0.1). HONEST: empirical only, no proof; TP2 likely STRICTLY
STRONGER than RH (PSD !=> TP2), tight margin => could fail at an untested config even if RH holds -- a sufficient-condition
CONJECTURE, robust where checked (§289-289g). Tools: operators/weil_tp*.out. Proof pointers: LEDGER §289-289j; origin §288c
(Groskin tail-TP) + §257e (arch=MP). Named tools: Karlin total positivity, log-convexity, sign-regular kernels.

## STATE OF THE ATTACK (synthesis of the §264–283b night — read this first)
**One target, sharply stated:** RH <=> Weil positivity <=> **lambda_max(prime rel. arch+pole) < 1, a-priori** (§277).
Measured comfortably < 1 everywhere (~0.29 worst genuine, 71% margin, §278); validated as GENUINE RH content (off-line
zeros break positivity, §279); highest at LOW height (§281), dominated by SMALL primes 2,3 (§280).
**Decomposed (proof skeleton, §281-282):** finite HARD CORE (low zeros x primes 2,3, ~0.29) + TAIL#2 (height, ~C/log→0,
arch growth, verified §282) + TAIL#3 (higher primes, super-exp suppressed, §280). Both tails controllable; the core is
the irreducible part.
**Structure of the core's difficulty (§283):** a DUAL-DOMAIN problem — F^2>=0 in r-space vs its FT SAMPLED at primes
{log n} in u-space (multiplier-vs-sampling). NOT pointwise-checkable (that shortcut is dead).
**Three ranked attack directions on the core (all rejoin it, per difficulty-conservation):**
  1. **Fourier-LP/SDP + certificate function (Cohn-Elkies-Viazovska), §283b-286** — newest, modern, high-ceiling; PARTIALS
     ALREADY IN HAND: explicit closed-form certificate proves Weil positivity for well-separated Gaussians (§284-285,
     arch-growth beats small-primes + Gershgorin). PRECISE TARGET (§286): the crude certificate reaches probe-spacing
     ~2-widths; an EXACT "magic function" certificate must close the spacing-2-to-0 gap (the dense regime = RH-core, where
     positivity currently comes only from the zeros' reality). This is the Viazovska analog: crude LP loose, exact modular
     certificate closes it. Automorphic (arch=MP/Gamma, §257e) or SDP-on-instrument-matrices. Most actionable + proof-like.
     **LITERATURE (§288, VALIDATED frontier-aligned) — engage these specific papers next:**
       - arXiv:2607.02828 (2026) "finite Guinand-Weil dictionary + archimedean tail order for the truncated Weil quadratic
         form" — DIRECTLY parallels our instrument+tail work (§271-282); their partial positivity + the arch<->p-adic
         "dictionary" barrier = our dense-regime barrier. Read for their restricted-class positivity technique + tail bounds.
       - arXiv:2005.02996 "Fourier interpolation with zeros of zeta and L-functions" (Radchenko-Viazovska-style) — the
         certificate/magic-function machinery for zeta zeros. The likely source of the exact-certificate construction.
       - Carneiro-Milinovich-Soundararajan (Fourier optimization + SDP for Guinand-Weil bounds) — the established SDP setup.
       - Connes "Weil positivity and Trace formula, the archimedean place" (Selecta) — our §276 arch-dominance mechanism.
     Our independent night reproduced this frontier's framework, partial results, AND barrier (§288) — the approach is the
     right bet; next real work = adopt these papers' techniques rather than re-deriving.
     **TWO ADOPTABLE TECHNIQUES (§288c, from Groskin 2026):**
       (A) INTERVAL-ARITHMETIC LDL^T factorization (mpmath.iv) -> upgrade our float PSD checks (§273-287) to CERTIFIED
           (rigorous) positivity. Our partials §284-285 become certified, not just measured. Concrete, near-term, bounded.
       (B) The arch tail is a "TOTALLY POSITIVE CAUCHY-STIELTJES increment" -- and THAT IS BECAUSE arch=Meixner-Pollaczek
           (§257e; MP is a Cauchy-Stieltjes/Hamburger-moment object). This is the RIGOROUS mechanism behind arch-dominance.
     **NEW LEAD (§288c) — TOTAL-POSITIVITY / sign-regularity angle:** since arch=MP => arch is totally-positive
     Cauchy-Stieltjes, ask whether the WHOLE Weil form (arch + prime perturbation) is a totally-positive / sign-regular
     structure. Total positivity is a strong PSD-implying property with rich theory (Karlin, sign-regular kernels). If the
     prime perturbation preserves enough sign-regularity, dense-regime positivity could follow WITHOUT the zeros. A genuinely
     new angle on the dense core, distinct from the certificate route, suggested by Groskin's tail + our §257e. Worth real
     exploration. Proof pointers: LEDGER §288c, §257e.
  2. **Geometric keystone (§265)** — the a-priori positivity = a Hodge-index theorem for Spec Z x_F1 Spec Z (Sonin
     cohomology intersection positivity); the deepest, Connes-Consani's program.
  3. **Beurling-Selberg extremal (§276b)** — classical 1D precursor of (1); gives partials, gap-to-1 = RH.
**Asset:** a machine-precision-validated (7e-16) self-contained Weil instrument (operators/weil_*.py, §271-273c) that
computes all of the above exactly, float-free. **Compass:** RH_DIFFICULTY_CONSERVATION_MAP.md (Keystone Specification).

## §283b — NEW MODERN ANGLE: Weil positivity is a FOURIER-LP/SDP problem (Cohn-Elkies-Viazovska)
§283 showed Weil positivity is a DUAL-DOMAIN problem: a positive function F^2 whose Fourier transform is constrained at
the primes {log n} -- NOT a pointwise inequality (the prime term SAMPLES FT[F^2] at the primes; that shortcut is dead).
This is EXACTLY the shape of Fourier linear-programming / SDP extremal problems -- the Cohn-Elkies/Viazovska method that
solved sphere-packing in dim 8 & 24 by constructing exact "magic function" certificates (modular forms). RH's a-priori
keystone = an exact certificate function certifying arch >= prime for all test functions.
- Sub-direction STRUCTURAL (research-scale, high-ceiling): is there an automorphic/modular "magic function" for Weil
  positivity? The arch=Meixner-Pollaczek/Gamma structure (§257e) + primes-via-{log n} suggests an automorphic certificate.
- Sub-direction COMPUTATIONAL: feed the validated instrument's exact matrices (Barch, Bprime, Bpole) to an SDP solver
  and study the dual certificate / extremal structure (note: finite-matrix PSD-check just re-confirms §273; the power is
  the CONTINUOUS certificate with exact sign constraints -- needs a real Fourier-LP setup, cvxpy/mosek).
This is a genuinely new, powerful toolkit for the project with a track record of cracking "hopeless" problems. Honest
caveat: Fourier-LP gives only bounds unless the certificate is EXACT (an exact RH certificate = the proof); none known.
Proof pointers: LEDGER §283-283b. Ties §276b (Beurling-Selberg = 1D precursor), §257e (automorphic arch).

## §281 — PROOF-STRATEGY SKELETON: the keystone's hard core is FINITE (concrete, actionable)
The operator-norm keystone (lambda_max(prime rel arch+pole) < 1 a-priori) DECOMPOSES (§278-281) into:
  (1) a FINITE HARD CORE: low-height test functions (first zeros, r~14-20, where arch Omega~0.8 is smallest) x the
      SMALL PRIMES (2,3 dominate, §280). Measured lambda_max ~ 0.29 there (71% margin). This is a nearly finite-dim,
      concrete positivity -- a candidate for a FINITE positivity/SOS certificate.
  (2) HEIGHT TAIL: lambda_max DECREASES with height (§281: 0.29 -> 0.11 as r0: 20->300), because arch Omega~log(r/2pi)
      grows. Control: show arch log-growth dominates the bounded prime interaction for r > some r0.
  (3) PRIME TAIL: higher primes are super-exponentially suppressed for smooth test functions (§280, g(log n)~e^{-w^2(log n)^2/2}).
      Control: bound the tail Sum_{p>P} by the FT decay.
=> PROOF STRATEGY: prove lambda_max<1 on the finite core (concrete, ~0.29, big margin) + bound the two tails analytically.
This is qualitatively more tractable than the monolithic "uniform over all test functions." NOT a proof, but a real
skeleton with the hard part isolated and measured. Caveat: the "core" is effectively-finite (super-exp decay), not
strictly finite -- a rigorous proof must still handle the infinite tail, but the tails are the CONTROLLABLE part.
Proof pointers: LEDGER §278-281. Tools: operators/weil_*.py (validated instrument).

## §277 — THE KEYSTONE AS AN OPERATOR-NORM BOUND (the precise target the skeleton attacks)
Weil positivity = RH <=> the prime perturbation's operator norm RELATIVE TO the archimedean+pole backbone is < 1,
uniformly over all test functions (on the full-rank/primitive subspace):
   lambda_max( prime  rel.  arch+pole )  <  1,   where
   prime[phi]=2 Sum_n Lambda(n)/sqrt(n) g_phi(log n),  arch[phi]=(1/2pi)int|phi_hat|^2 Omega,  Omega=Re psi(1/4+ir/2)-log pi>0 at height,
   pole[phi]=|phi_hat(i/2)|^2+|phi_hat(-i/2)|^2.
MEASURED with the validated instrument (§277/§277b): lambda_max ~ 0.20 (low band) to 0.38 (wide probes, tightest genuine)
to 0.28 (height) -- comfortably < 1, robust arch-dominance. (lambda_max=1 ONLY in rank-deficient configs = #probes>#zeros
= the free-side resolution limit, NOT positivity failure -- §277b.) NO RH-critical regime exists (RH true; §276b).
=> THE KEYSTONE (research-scale): prove lambda_max < 1 UNIFORMLY, A-PRIORI (without the zeros). This is the operator-norm
form of "arch backbone beats the prime fluctuation" = the arithmetic-surface Hodge-index positivity. Concrete attack
directions: (a) is it a Beurling-Selberg extremal bound (the sharp constant for prime-vs-arch)? (b) does the MP/de Branges
structure (arch=Meixner-Pollaczek, §257e) give the norm bound via a classical OP inequality? (c) the geometric route
(Sonin cohomology intersection positivity). All rejoin this one bound. Tools: operators/weil_lowner.out,
weil_lmax_sweep.out. Proof pointers: LEDGER §277-277b.

## §271–275 — A VALIDATED SELF-CONTAINED WEIL INSTRUMENT + first results (a real tool now exists)
Built and machine-precision-validated (7e-16 vs the direct zero-Gram, §273c) a from-scratch high-precision Riemann-Weil
quadratic form: Q(phi)=|F(i/2)|^2+|F(-i/2)|^2 + (1/2pi)int F^2 Omega dr - 2 Sum Lambda(n)/sqrt(n) g(log n), Omega=Re
psi(1/4+ir/2)-log pi (§271-verified). Tools: operators/weil_selfcontained.py, weil_zerobasis.py, weil_margin.py,
weil_resolution.py. GOTCHAS baked in: (a) test functions must OVERLAP the zeros (gamma>=14) or the form is ~0
(zero-blind, §272); (b) mp.quad needs BREAKPOINTS at each Gaussian bump or it under-resolves far-from-origin bumps (§273c);
(c) the +/- signature is generically ~half-negative and basis-invariant (Sylvester) -- NOT where the zero info is.
FIRST RESULTS: §274 Weil-positivity MARGIN (min eig) converges to ~0.878 as the band widens (robust positivity, no
near-critical direction for well-separated probes). §275 the zero-count lives in the RANK: rank(band-restricted Weil form)
= N(band) (# zeros in band), saturating regardless of #probes (extra probes -> exact-null directions). Parallels
dim H^1 = #Frobenius eigenvalues (the §266-267 template's cohomology-dimension slot).
LOADED FOLLOW-UPS (de-risked, trivial with the instrument): (1) sweep band height T, confirm rank grows as
N(T)~(T/2pi)log(T/2pi) -- the prime-built form reproduces Riemann-von Mangoldt counting; (2) push probes DENSER than the
local zero spacing and study the margin-collapse threshold as the zeros' rigidity limit; (3) the REAL RH use: the arch+prime
form is a FIXED prime-built object whose PSD-ness = RH (Weil positivity) and whose rank = N(band) always -- study how the
arch+prime split realizes the positive rank structure (that split is where RH lives). Proof pointers: LEDGER §271-275.

## §264–266 SHARPENING — the link is precisely specified
The difficulty-conservation law (below) stands. Tonight's overnight arc pinned the ONE permitted attack to an exact target:
- **§264:** a forward-direction monotone MUST be prime-coupled — any functional of only the zero positions is RH-blind
  (off-line response is Pythagoras, not arithmetic). Kills the "cleverer position-functional" sub-lane.
- **§265:** anchor-to-exemplar (elliptic-curve Frobenius = zeta_E zeros, |α|=√p verified) ⟹ the proven case gets its
  positivity from INTERSECTION THEORY on a surface (Hodge index), NOT operator eigenvalues. So §236/§258/§260's
  operator-spectrum probes were the WRONG OBJECT. Target: RH ⟺ a Hodge-index theorem for Spec Z ×_{F₁} Spec Z
  (scaling flow = Frobenius, Sonin space = H¹).
- **§266 (the pin):** verified the exact 2×2 template — cup product J, Frobenius scales it by p (correspondence),
  invariant pos-def metric M with F/√p an M-isometry (reality), all machine-precision. KEY: "M exists" ⟺ |eig|=√p is
  CIRCULAR; the function-field proof is non-circular ONLY because the Hodge-index theorem furnishes M **a priori from
  geometry, independent of the eigenvalues**, then forces them onto the circle.
- **⟹ the link = an a-priori, zeros-INDEPENDENT, scaling-flow-invariant positive-definite intersection metric on the
  Sonin space.** Flow-invariance (b) we can already exhibit (§236 skew generator); the a-priori positivity source (a)
  is the sole missing keystone and must come from geometry (an index theorem), not from diagonalizing an operator.
- **LOADED NEXT PROBE (bounded) — the PRIMITIVE-subspace signature vs λ (corrected form, see §270/§270b):** §266b found
  the RAW Weil-form signature robustly INDEFINITE (surface-H¹-like ✓); §270 found it generically ~half-negative
  (nneg/nmodes~0.5) with NO zero-count structure. §270b (Sylvester's law): the signature is BASIS-INVARIANT, so a
  Hermite/Gaussian basis CANNOT change it — do NOT bother rebuilding in another basis expecting a fix (that was a
  misconception, now retired). THE GENUINE OPEN QUESTION: compute the signature specifically on the PRIMITIVE/constrained
  subspace {f(0)=fhat(0)=0} and track its negative-count vs λ — does nneg there STABILIZE to a fixed count as λ grows?
  A fixed primitive-negative count independent of λ WOULD be a zero/prime-tied invariant (the RH-relevant "degree"). The
  obstacle is real: at finite λ + float the primitive modes sit at the near-radical floor (§266b fragile), so this needs
  a regularization that survives — higher-precision mpmath at small λ, or CCM's λ→∞ quotient-by-range(E) machinery (§258),
  NOT a basis swap. Use operators/dbn_dynamics/weil_gram_constrained.py's kernel. Proof pointers: LEDGER §264–270b.
- **§269 vein RESULT (closed, feeds the above):** the de Branges coefficient level sees primes at the p⁻¹ scale
  (§269e: 1−f(p)~4/p, the arch measure moment-averages the oscillatory p⁻½ Weil term away) — structurally BELOW the
  Weil level, so the free half / coefficients are provably RH-blind. Confirms §268 from a second angle. ⟹ the signature
  experiment above must be at the WEIL (log-derivative) level, which weil_gram_constrained.py already is — good.

---
## THE STRATEGIC HEADLINE (§256–260 + §264–266) — full context below

## THE STRATEGIC HEADLINE (§256–260, autonomous night) — read this before choosing any direction
**The DIFFICULTY-CONSERVATION LAW** (RH_DIFFICULTY_CONSERVATION_MAP.md): every framework gets ONE of {reality,
correspondence} free; the other = RH. Verified across all our frameworks; every concrete crux (incl. the Sonin-space
polarization) collapses to Weil positivity or its mirror. ⟹ **single-framework reformulations provably cannot cross
the wall** — STOP spending on them (more positivity criteria, operator bases, Beurling systems, decay characterization).
The crack REQUIRES supplying the reality↔correspondence LINK a-priori (function-field polarization analog). Live
incarnations, ranked: (1) Spec Z×Spec Z / F₁ = CCM arithmetic-site / Sonin cohomology + a scaling-invariant
polarization [most concrete; needs CCM's regularization machinery — naive builds are truncation-dominated, §260];
(2) a NEW forward-direction dynamical monotonicity (RT reach only Λ≥0); (3) explicit MP canonical-system H(t)≥0.
**Key new asset (§257e):** CCM's arch cyclic pair = MEIXNER–POLLACZEK (λ=¼) EXACTLY — their frontier arithmetic is
classical/computable (discrepancy d_n; primes = index-frequencies → Nyquist). The Geronimo–Nevai closed form of d_n
is the concrete unlock. Tool: operators/ccm_discrepancy.py.
Honest: Lane B thoroughly explored + consolidated; the genuine cracks are the deep link-programs (fresh careful
session, not marathon-forcing — tonight hit 7 verify-catches). Lane A (arc-Toeplitz ‖D‖ finite theorem) still bankable,
room's lane.

## (earlier) LEDGER through §254 — the ARC COMPLETE + §253b/§254 blocks

## ARC COMPLETE — the phase-function / superresolution lens (§241-252)
The Suzuki spectral-flow arc reached its root. Deliverables: exact-Q pipeline (1e-13, a>1 wall gone, 700×);
phase function Φ=angle(−u/w) explains θ-universality; **measured a²-superresolution Φ'∝a²** DERIVED (§251) as
φ'=1/max(ε,δ) with Suzuki resolution ε~1/a², criticality-locked (§248/251); **converges to the ALREADY-
ESTABLISHED prolate near-null** (§98: cluster=prolate angle spectrum, ~a² modes) — §252. So it's a new MEASURABLE
characterization of the known prolate wall, NOT a new crack. 2 self-caught overclaims (2π² gauge §248b, off-line
discriminator §250 — Suzuki self-adjoint ⟹ zeros always real). The φ'=1/max(ε,δ) discriminator is real but self-
adjoint-blind. **§253b DONE — DECISIVE POSITIVE RESULT on the orthogonal survivor.** The near-null telescoping cascade
(§176b: h(n) ≡ −residual, ratio≡1.000, real primes cascade to λ_min~10⁻³⁷) is SPECIFIC to the integers.
Both Beurling systems ({2^k}, shifted-position) HALT: ratios scatter (1.35→12.76), λ_min stays O(1)/
indefinite (RH-false positivity failure). Ran LOCALLY — connes_cvs is vendored (cvs_venv/), NOT Olympus-only —
c=9 dps=100, ~2 min. Tool: operators/dbn_dynamics/beurling_cascade_cvs.py. This is the arithmetic-coherence
fingerprint separating RH-true from RH-false = "what a proof must consume" (§182). A GENUINE live crack-candidate
(not a rejoin), likely = prolate near-null (§252) sampled at integer frequencies ω_n=1−log n/L.
**§254 NEXT: derive the telescoping identity h(n)=−R_{n−1}** — equivalent to the super-geometric decay law of
h(n)=Λ(n)/√n·⟨ψ|sin(2πxω_n)|ψ⟩ for real primes. Test the prolate-sampling explanation (measure h(n) vs ω_n vs
prolate eigenvalue decay). Whatever the decay needs beyond PNT is the RH-content Beurling lacks. Bottoms out in
the prolate near-null core — but now with a CLEAN falsifier-discriminated measurable.

## (prior) staircase/dilation arc (§241-246, one night, Fable + exact pipeline)
The live frontier moved. Suzuki operator D_{a,θ} spectral flow in a: branches PARK at zeta ordinates,
θ-robustly (§241). Exact-Q pipeline (suzuki_exact_Q.py — closed-form pieces, validated 1e-13, a>1 wall
GONE, λ_a = machine-zero knife-edge through a=3.0, ~700× faster) (§242). Universal law v_min = c(a)·γ/a,
zero- AND θ-independent; c: 0.083→0.014 over a∈[1.1,3.9] ~ a^{−1.44} (§243). Per-prime force ledger:
kinematic = free-scaling −z/a exactly; restoring = collective prime interference (§244). Two-law split:
best-parked offset ~ e^{−1.26a} (missing-prime tail) vs sliding ~ γa^{−2.4} (window drag) (§245).
MECHANISM: parking = asymptotic DILATION-COVARIANCE of the prime comb; c(a) = the covariance defect;
RH-discriminator = an off-line zero's width breaks the c(a)-collapse at its ordinate (§246).
§247-248b DONE (with a self-correction). Phase function Φ (e^{iΦ}=−u/w, |u|/|w|=1.000000 exact):
θ-universality EXPLAINED (all extensions = level sets of one Φ). **Robust law: v ∝ γ/a², Φ' ∝ a², E-zero
depth ∝ 1/a² — QUADRATIC SUPERRESOLUTION** (vs linear explicit formula). §247's "2π²/1/π²" constants
RETRACTED §248b — they were λ-shift(=0.2)-artifacts; the a^{−2} EXPONENT is the gauge-invariant law,
prefactor C(shift) is not. Cross-zero uniformity (equal weight, ratios 0.87-1.10) + θ-universality are
gauge-free. β-deformation (§248): 3-channel discriminator {bulk=2a}∧{peaks a-stable at γ}∧{Φ'∝a²} all
fail together off criticality (λ_a max EXACTLY at β=1). Knife-edge: entire low spectrum of A_a machine-zero
(fully degenerate near-null) — zeros live in the E-zero/deficiency structure, NOT resolved eigenvalues.
**§249 NEXT — derive the a² superresolution exponent** from the deficiency element v̂₊(z) asymptotics,
shift-explicit (C(shift) ties to Suzuki Cor 1.6's own RH-tuning parameter). The exponent is theorem-shaped;
don't chase a "nice" constant (that was the §247 error). RH-shape: off-line pair = Φ'-peak that floors
(width → ζ-width) instead of shrinking 1/a². Tools: suzuki_phase_density.py, suzuki_exact_Q.py (1e-13).

# (prior) NEXT SESSION (2026-07-24 LATE, LEDGER through §238 — read §234-238 first, then §209-220)

## CURRENT STATE (2026-07-24 late — TWO LANES running in parallel)

**LANE A — BANKING (finite unconditional theorem, ~1 lemma from done).** A fleet "room" (Ein/Isa/Rie/Hawk)
made the tail lemma rigorous overnight. Both sectors' low blocks NOW Arb-certified: EVEN K0=23 (mode-1 pivot
[4e-8±4.35e-9], re-run on Mac too), ODD K0=24 (mode-2 [1.12e-5±3.96e-8]) — the odd sector's FIRST real
computation (the old "piece 3 done" was a mislabeled even-block minor — vacuity caught + repaired §239-class).
Tail architecture ("4′" split at cut M*) closed-form: u_c=−1 exact, β=4w_p/π, diagonal law log(kb/2π)=2πN'(kb)
(Isa), coupling HS-bound ½√log2 (Isa). ONE thing gates the block size M* (∈[~50, ~10⁵]): the β-channel operator
norm ‖D‖=0.743432. I reframed it (§Sam room turn): the Toeplitz kernel sin(θm)/2m has symbol EXACTLY (π/2)·𝟙_arc
(measure 2θ) — so ‖D‖=(π/2)·‖Herm[modulated arc-Toeplitz]‖, a classical Kac–Widom object; answer key 0.47331
normalized. Rie owns route 1 (sawtooth-symbol), Ein owns M*+interval-Schur composition. §222/§225 rescoped to
even-only then odd repaired: both-parity certified reach [0.549,0.64] now GENUINELY two-sector (odd certifies to
0.70). Watch: gate tools exit-0-on-fail bug (Hawk) + un-audited DONE rows (piece 5 pending).

**LANE B — RH FRONTIER (mine, §234-238).** The §234 pivot: proven mechanism is Hilbert-Pólya (skew Θ), not Weil
positivity. Arc:
- §235-AKM: read AKM 2410.20758 first-hand. Skew-adjointness needs ONE axiom — flow ISOMETRIC on leaves; then
  Θ=leafwise Lie-derivative on Ker(Δ_ℱ), spectrum⊂iℝ, reality AUTOMATIC (no positivity). Metric we need = the
  scaling-isometric leafwise metric (named condition, not blind SDP).
- §236: RAN it (hp_scaling_compress.py). Compress dilation gen Θ₀=diag(2πin/L) to L² near-null → exactly skew
  (1e-14, AKM mechanism realized) BUT spectrum = bent LATTICE, not zeros. Obstruction pinned to the METRIC.
- §237: obstruction NAMED = log-growing Riemann–von Mangoldt density incompatible with any fixed lattice; unifies
  §232 (no scale-invariance) + Isa's diagonal law + AKM Weyl density as ONE fact. (Caveat: secular root-finder
  dropped sin(zL/2), roots muddied — only the density/sign facts are robust.)
- §238 HONEST consolidation: density-mismatch was ALREADY known (RH_MAP.md:27, berry_keating_deficiency.py "must
  be NONLOCAL periodic-orbit BC"). AKM's genuine delta = a PROVEN 2024 template for that nonlocal operator + §236's
  finite witness. LIVE EDGE: build the candidate OPERATOR (dilation + prime periodic-orbit BC, spectrum→zeros with
  log-growing density, skew by construction) — the "no candidate operator" the project marked missing. connes_
  decisive gives zeros as a VECTOR; the missing object is the OPERATOR. §239 next.

---

# (prior) NEXT SESSION (2026-07-24 pm, LEDGER through §220)

## CURRENT STATE (2026-07-24 'reimann6' session — second-bridge hunt, cluster closure, certified extension, 2×2 reduction)
The definitive orientation. Read LEDGER §209–§220 for detail; this is the map.

**What got CLOSED this session (all with concrete diagnostics, not vibes):**
- **Operator/prolate bridge cluster** — B1 frame (CM 2310.18423 has NO lower-bound thm, workflow confabulated it, §212),
  B2 determinant (exponential/prolate spectrum not Fisher-Hartwig, §211c), B3 SUSY (superpotential is prolate
  ground-state log-deriv not ζ'/ζ, §213). All reduce to the prolate near-null (§213b). Extends the A1-A11 exhaustion.
- **GM Lens B** (§215) — the last classical goal-relevant sub-problem: GM's large-value method is a σ≈3/4 phenomenon,
  the mollifier barrier is σ=1/2 Kloosterman — no transfer. Structural negative.
- **Norm-budget RELAY** (§218) — quantified (arch=−prime to 6 digits on the near-null) but delivers only the leading
  cancellation; the RESIDUAL (~1e-8) is the RH content = the prolate/Fuchs defect. The relay reduces to the residual wall.
- **The "race" / Onset Lemma** (§219, §219b) — REAL new structural fact: each prime enters as a rank-one NEGATIVE ε³
  square in the universal direction ℓ_k=k (κ≈21.45). BUT the drain is ⊥ the near-null (overlap 0.005), so it does NOT
  control the margin — the race frame doesn't crack it.

**The ONE forward RESULT (verified §223):** first CERTIFIED crossing past the first prime — multiprime cert built +
gate-verified (`operators/track_crossblock/arb_certify_multiprime.py`, reproduces the p=2 cert exactly). Even-4×4 +
odd sector both rigorously PD across the p={2,3} window (even: 11 pts a∈[0.55,0.62], stable ~3-4e-8 margin at
NLERCH=50000; odd: full window [0.55,0.80]). Extends the theorem from the p=2 endpoint (a=½log3=0.5493) into the
p={2,3} window. NOTE: §217b's "even fails past 0.57" was a NLERCH-precision artifact, RETRACTED (§223). Scope: still
FINITE-GALERKIN (mode-tail lemma for p={2,3} open/delicate), weaker than the minimal theorem; entries have no closed
form (§223) so it's Arb-certification not symbolic proof.

**The CRUX, razor-sharp (§220):** the residual reduces EXACTLY to a 2×2 determinant — margin sign = sign(det S),
S = Schur{1,3}|{5,7}, det S = S₁₁S₃₃ − S₁₃². Product-emergent (det-arch and det-prime both negative; positive only
from the combination) = the §182 wall at 2×2 resolution. The finite margin is ~1e-9 (high-Arb-precision only). The
RH-limit is the FULL spreading near-null (23-dim §107), NOT the 2×2 — the 2×2 is a finite-window tool.

**IN FLIGHT (2026-07-24 pm):** (a) workflow wadbh9t2j — 6 parallel high-prec Arb certs across [0.565,0.615] to find
the rigorous finite even-block positivity BOUNDARY + symbolic recognition of the entries; (b) agent — grounding B4
(Lax-Phillips scattering), the LAST unexplored sign-carrying mechanism candidate (§182 spec vi). Synthesize both, log,
follow the path.

**HARD-WON META:** everything reduces to signing the product-emergent residual via a mechanism ≠ the (sign-blind)
explicit formula. Operator-side is multiply-exhausted. Forward = the genuinely-new mechanism (B4, or a fresh
non-operator idea) OR banking finite certified reach (Arb-grind the 2×2). Tools: operators/track_crossblock/
{arb_certify_multiprime.py, relay_decomp.py, onset_probe.py}. Verify-before-ship caught 3 agent overclaims + 2 own +
2 confabulations this session — the §207 gate is load-bearing; keep it.

---
_(older §181 handoff below — SUPERSEDED by the above; retained for provenance)_

# NEXT SESSION (updated 2026-07-23, LEDGER through §181)

## CURRENT STATE (§181) — every basin is mapped + walled; STOP basin-hopping (verified, not assumed)
Read the LEDGER, not this doc's older sections, for basin status — they are STALE. Verified closures:
operator-positivity (§122, six faces) · zero-density estimates (§180, route-3, not-an-RH-path) · mollifier proportion
(route-7, θ→1 = RH-hard) · **DBN Λ≤0 (§35-47e — a Λ-DETECTOR, no prime-handle; NOT "dormant/unmapped")** · **LeClair
spectral-flow (§125 — = Lagarias' known RH-equivalent)** · Jensen/moment/de Branges/Beurling-Nyman (OBSTRUCTION_MAP 1-5).
Do NOT re-enter these as if fresh (my §180 nearly re-treaded DBN; search-before-build/read §35-47e,§125 caught it).

**The honest forward set (§181), goal-relevant only:**
1. **§126 Lens B (the one concrete live sub-problem):** can Guth-Maynard's 2024 large-value/energy method (singular
   values of the prime matrix M_W) transfer to push the Levinson-Conrey mollifier past θ=4/7 (proportion past ~42%)?
   §180 verified GM's RESULT is pure density (doesn't touch the mollifier); the METHOD-transfer is open/possibly-
   unclaimed. Hard analysis, not a probe. Real progress if it works (proportion, not full RH unless θ→1).
2. **Non-classical cross-basin ideation** (§125 rule-1; superPRISM on Fable) — the only route to the genuinely-new
   idea the verified wall demands.
3. **Consolidate** the §176 cascade (novel analytical near-null anatomy, §175-flagged unpublished) + the maps.
Per riemann-goal (prove RH, not partials): 1 & 2 are the goal-relevant moves; 3 is banking.

---
_(older handoff below — LEDGER-through-§123, RETAINED for provenance but basin-status is SUPERSEDED by §181 above)_

# NEXT SESSION (updated 2026-07-21, LEDGER through §123)

## STRATEGIC PIVOT (§123): operator-positivity basin EXHAUSTED → enter the ZERO-DENSITY basin
Forks A+B (parallel) settled the strategy. FORK A (§122): attacking the limit-positivity directly is the RH-hard core,
no non-circular handle (arch O(−e^a) ↔ prime +O(e^a) cancellation is vivid but positive⟺reality⟺RH). Our operator-
positivity basin (§1-122, six equivalent faces §121) is exhausted at ONE wall. FORK B (§123, survey): the genuinely-
DIFFERENT basins (not our wall) are, ranked: ① ZERO-DENSITY/Levinson (obstruction=mollifier Cauchy-Schwarz losses,
GENUINELY different; ≥40.88% of zeros proven on-line — real unconditional progress; ACTIVE: Guth-Maynard 2405.20552
2024); ② Spectral-flow/LeClair 2406.01828 (S-matrix unitarity→real zeros, novel, unproven); ③ Hilbert-Pólya (density
mismatch); ④ dBN Λ≤0 (dormant). SAME-WALL (avoid): function-field/Connes arithmetic-site, moments/RMT, Beurling-Nyman.

## Immediate next block — enter the zero-density basin
1. **READ Guth-Maynard 2405.20552** (2024, new Dirichlet-polynomial large-value estimates → improved zero-density
   N(σ,T)) + the Levinson-Conrey mollifier method (how 41% on-line is proven; where the Cauchy-Schwarz loss is).
2. Find the concrete tractable sub-problem: improved Dirichlet-poly bounds, or a mollifier-free gap argument toward a
   higher critical-line proportion. This is a genuinely DIFFERENT obstruction (analytic, not positivity-at-limit).
3. High-risk alternative to keep in view: LeClair spectral-flow (rigorously justify S-matrix unitarity as Re(s)→½).

## What the operator-positivity campaign leaves (durable)
The comprehensive six-face map (LEDGER §121), tools (compact G_a §110, mpmath hi-prec §109, prime-bend §112, semilocal
Jacobi §119), and the honest closure: no operator-positivity face is a shortcut past the a→∞/ω→0/all-primes limit.
Don't re-enter this basin. Self-corrections banked (§92, §117): check definitional properties (|m|≡1) before overlays.

## Orientation
LEDGER §92-123, PROCESS/STATUS 07-21, memory rh-suzuki-screw-function-frontier. Sources in membrane_pulls/. Tools:
operators/dbn_dynamics/ (many). Model: Opus for the zero-density reads + work; Fable = brainstorm only.
Parked: REBOUND Λ-hunt; Rick offline (SpaceEngine specimen).

---
_(superseded detail below from the §106 redirect — kept for provenance)_

# NEXT SESSION (updated 2026-07-21 step-back, LEDGER through §117)

## State after the §117 retraction (read this first)
The §106 redirect (off the k_λ wall onto Suzuki's W/Cor-1.6 door) STANDS, but §113-116's "reality transfers for free
via Herglotz" was RETRACTED in §117: m_a = v̂_+/v̂_- is UNIMODULAR/INNER on R (|v̂_+|=|v̂_-| exactly, forced by Thm 1.5),
not Herglotz; inner limits do NOT forbid C+ zeros. The correct target in this framework is the classical one:
**E_ξ = Ξ + iΞ′ is Hermite–Biehler ⟺ RH** (equivalently Θ_ξ=E_ξ^♯/E_ξ inner in C+).

## The honest map: ONE wall, four faces (all equivalent, all = RH)
1. k_λ ≈ θ_x rigorously (CCM 2511.22755 §7 — the worldwide stuck point)
2. θ(a) selects the right self-adjoint extension in the a→∞ limit (Suzuki Cor 1.6; deficiency (1,1) = limit-circle)
3. E_ξ Hermite–Biehler / Θ_ξ inner (classical)
4. momentum→zeta spectral transition (the prime-bend §112 made qualitative progress on this face)
Everything this session touched lands on one of these faces. Do NOT re-derive the equivalences — they're mapped.

## What genuinely stands from §104-117
- Suzuki Thm 1.5: finite-a reality + simplicity UNCONDITIONAL (his theorem, source-read). The door is real.
- Compact G_a reformulation (§110, g_a_compact.py): breaks the large-a conditioning wall, clean λ_a>0 to a=3 float64.
- High-precision run (§109): §42's negativity = float64 artifact; even-simple robust; two floors separated
  (precision fixed; basis-resolution structural).
- Prime-bend (§112): primes inject zeta-like gap-irregularity (CV 0.393 vs zeta 0.389 at a=2), density-independent.
- S_t(z)/P_t(z) explicit (§110, from [14] eq 1.5-1.6) — the Parseval identity ∫S_x S̄_y dz = g̃(x,y) (RH-equivalent,
  §107b) is fully computable and NEVER attacked analytically. Genuinely untouched face.

## Next block (recommended order)
1. **READ CM 2310.18423** ("Zeta zeros and prolate wave operators" — semilocal prolate = arch+primes; queued twice,
   never read; the paper working the SAME arch↔prime structure as our §104).
2. **READ Suzuki 1204.1827** ("A canonical system of differential equations arising from ζ") — the explicit canonical-
   system Hamiltonian H(t) for Ξ; plus the de Branges subspace CHAIN (total ordering, 2606.09096 §5.7) — the chain IS
   the canonical system in inverse spectral theory, and we have NEVER explored it. HB(E_ξ) ⟺ H(t)≥0 is the classical
   reformulation — find where the primes sit in H(t).
3. Then: attack the **Parseval/completeness identity** (face nobody's worked): ∫S_x(z)S̄_y(z)dz = g̃(x,y) with S_t
   explicit — completeness/frame-theory tools, not approximation.
Method rule sharpened (§117 lesson): before claiming object X has property P, run P's DEFINITIONAL one-line check
first (|m|≡1 would have killed the Herglotz claim instantly). Magnitude checks before interpretive overlays.

## Orientation
LEDGER §92-117 (the full arc incl. retraction §117), PROCESS 07-20/21, memory rh-suzuki-screw-function-frontier.
Sources in repo: suzuki_2606.09096.txt (Thm 1.5 l.385, §7.2 l.2206, §8 l.2420), suzuki_2301.00421 (S_t eq 1.5-1.6),
rt_1801/, connes_2602/. Tools: operators/dbn_dynamics/ (validated: g_a_compact, g_a_W_prime_bend, mpmath_weil_hiprec;
confounded-documented: lp_lc_test, weyl_mfunction). Unpulled: 2310.18423, 1204.1827, 2112.05500.
Parked: REBOUND Λ-hunt (Olympus); Rick offline (SpaceEngine specimen unchecked).
Model: Opus for the reads + analytic work; Fable = brainstorm/step-back only.

═══════════════════════════════════════════════════════════════════════════
# STATE AS OF 2026-07-29 (Riemann7 overnight autonomous, §411-438)
═══════════════════════════════════════════════════════════════════════════

## THE PICTURE (rock-solid, grounded 4 ways)
RH decomposes cleanly (§434): **functional-equation SYMMETRY** [= 1D Arakelov RR on Spec ℤ = the θ/ζ functional
equation = HAVE IT, classical] + **POSITIVITY** [zeros exactly ON Re=½ = the 2D self-square Spec ℤ×_{F₁}Spec ℤ
diagonal Hodge-index positivity = THE OPEN STEP]. Every route (Weil form, Circle-Lemma, Connes-Consani, δ_p,
Nyman-Beurling) reduces to this same open positivity-half.

Proven this session:
- The analytic Weil instrument has NO 1-D shortcut (7 verify-catches killed every apparent one). It faithfully
  re-derives Connes-Consani 2020 (archimedean place = THEOREM); the prime-place combination is the open frontier.
- Size/effective bounds provably NEVER suffice: ρ_edge=||prime||op/archmargin ∈[1.2,2.0] at ALL heights (§430,
  extended to 150 zeros). RH cancellation is irreducible at every γ.
- The wall is GLOBAL UNIFORMITY across bands, not weak local signals (§431: off-line detectable height-indep).
- Independent cross-check: Nyman-Beurling d_N²logN→2+γ−log(4π)=Σ1/|ρ|² matches 150-zero data (§433).

## THE OPEN STEP (Connes-Consani's own decade-open problem)
Construct an intersection pairing + manifestly-positive Hodge form on Spec ℤ×_{F₁}Spec ℤ (self-square). The
F₁/tropical foundation is COMPUTABLE (min-plus, scaling site); the positive Hodge form is the missing algebraic-
geometric piece. Latest tools: C-C 2024 "Riemann-Roch for ℤ" (ℤ=𝕊[X], 1+1=X+X²), C-C 2026 "Absolute Geometry of
Spec ℤ" (arXiv:2606.06604), Morishita 2026 Deninger↔C-C bridge, "Zeta Spectral Triples" (C-C-Moscovici 2026,
targets diagonal self-intersection ↔ zeros).

## LOADED NEXT ATTACKS (not started — need care/steer)
1. A CORRECT toy of the C-C tropical self-square (needs the precise min-plus intersection conventions pinned from
   the papers first — do NOT guess; §423 showed naive re-indexings give regularization artifacts).
2. Euler-product-specificity of positivity: does a functional-equation-only function (Davenport-Heilbronn) have
   the Weil form BREAK, isolating multiplicativity ⟹ positivity? (§411 partial: bent primes break it.)
3. The "Zeta Spectral Triples" diagonal-self-intersection ↔ zeros construction — read the 2026 paper, find a
   computable piece.

## ASSETS: zeros150.npy (150 zeros, γ<319); extended-instrument mode (l2_fejer with ZEROS override); Báez-Duarte
exact-Gram tooling (rt6_baez_duarte.py); primes-from-zeros recovery (rt6_primes_v2.py). LEDGER §411-438 full detail.

═══════════════════════════════════════════════════════════════════════════
# UPDATE 2026-07-29 (morning autonomous session §442-448) — GEOMETRIC ROUTE MAPPED
═══════════════════════════════════════════════════════════════════════════
Drove the C-C tropical/geometric self-square (the loaded #1 attack). Result: the POSITIVITY side is realized;
the gap is precisely the CORRESPONDENCE.
- §442: verified Weil's F_q Hodge-index MECHANISM (intersection form on {f₁,f₂,Δ}, signature (1,2), Castelnuovo-
  Severi D²≤2ab ⟹ |α|=√q). The exact template.
- §443 RESEARCH: scaling site conventions pinned; the tropical/matroid Hodge-Riemann IS PROVEN (Adiprasito-Huh-Katz
  arXiv:1511.02888, signature (1,rank-1)); the C-C gap is narrowly "does the F₁ self-square inherit it." Subtlety:
  ζ's infinitely-many zeros come from the ℕ*-structure (log p ℚ-independent ⟹ infinite genus), not a naive genus-1 H¹.
- §444-448: built the number-field Rosati form A[r,s]=Σ_ρ[e^{ρr}e^{(1-ρ)s}+swap] through a verify-driven arc (missteps:
  too-small-M §444, wrong q-diagnosis §445, conditioning-worry §448 — all resolved). LANDING: the FE-closed Rosati
  VALIDLY detects off-line (well-conditioned, M=20 sees β=0.1), PSD ⟺ RH in the M→∞ limit. The geometric positivity
  criterion WORKS for ℚ.

## THE FRONTIER IS NOW A SINGLE OBJECT: the correspondence Θ (Deninger's phantom / C-C scaling correspondence).
Positivity given the eigenvalues {e^ρ} is fine (§448 Rosati + §443 proven tropical Hodge index). RH ⟺ Θ EXISTS as a
geometric correspondence on the F₁ self-square with Θ·Δ = explicit-formula trace. This is Deninger's open operator /
Hilbert-Pólya. The MODEL where it exists: Selberg trace formula on hyperbolic surfaces (Laplacian = Θ, self-adjoint
⟹ "RH"). ζ's analog = the missing arithmetic dynamical system.

## LOADED NEXT (Θ-attack): (1) Selberg-comparison — what would ζ's "Laplacian/dynamical system" need vs a real
hyperbolic surface. (2) Morishita 2026 Deninger↔C-C bridge (arXiv:2508.15971) — computable piece? (3) intersection-
number consistency: are Θ's explicit-formula-required numbers consistent with a geometric correspondence?
Verify-catch count now ~11 (added §444/§445 missteps, all self-corrected).

═══════════════════════════════════════════════════════════════════════════
# UPDATE 2026-07-29 (session §449-460) — GRAND UNIFICATION + obstruction characterized to the bone
═══════════════════════════════════════════════════════════════════════════
See RH_MAP_2026-07-29.md for the full navigable synthesis. Headlines:
- GRAND UNIFICATION (§449-450): Θ = Berry-Keating xp (1D; density log r matches zero counts exactly). Analytic
  (Weil form) = Geometric (Hodge/Rosati) = Spectral (xp/Connes) — ONE open core: Weil-Connes positivity =
  xp self-adjointness with the arithmetic = my instrument's full≥0.
- 5 RH-equivalent faces all confirmed consistent & tied to the zeros: Weil, Nyman-Beurling (§433), Li (§441),
  Mertens/Möbius (§455), ψ(x)−x (§458).
- SEMI-LOCAL COMBINATION (the C-C open step) characterized (§451-457): non-incremental (partial prime sums
  indefinite), arch dominates prime mode-by-mode (§456), the crux is the off-diagonal coupling which is full-rank/
  spread/razor-tuned with NO structural shortcut (§457). Only a global mechanism (Hodge signature / operator
  self-adjointness) can work — incremental place-by-place proofs provably can't (matches §422).
- Size bounds NEVER suffice: ρ_edge>1 to γ~670 (§459, 500 zeros). Only phase/cancellation reaches it.

## SHARPEST LOADED ATTACKS (for Θ's existence — the sole remaining gate)
1. The F₁ self-square ↔ matroid/tropical-surface BRIDGE: AHK Hodge-Riemann is PROVEN for matroids; does Spec ℤ's
   F₁ self-square (infinite genus = ℚ-independence of {log p}) inherit it? The precise open construction.
2. "Zeta Spectral Triples" (Connes-Consani-Moscovici 2026): diagonal self-intersection ↔ zeros — read for a
   computable piece.
3. xp self-adjointness with the adelic domain (Connes): the operator whose spectrum IS the zeros.
All three are the genuine Hilbert-Pólya/Connes/Deninger frontier — new mathematics, not more computational probing.

## HONEST STATE: the computational investigation is comprehensive. The core is the actual 165-yr open problem
(needs the structural construction of Θ). The map is complete, cross-checked, and the positivity side is
realized/available; the sole gate is the arithmetic object's existence. Assets: zeros150.npy, zeros500.npy,
full instrument suite. This is a precisely-located frontier, not a wall — the next move is the structural build.

---
## UPDATE 2026-07-30 — roundtable S(T) handoff discharged + extension no-go generalized (§466–468)
The CCM archimedean prolate operator was built (§464, ~3% / O(1) residual) and the roundtable then
Arb-certified: deficiency (2,2)/U(2) extensions (§424), Weyl density = RvM exactly & θ-invariant (§425),
δ_n=O(1) (§426), bare fluctuation right-size/wrong-shape vs S(T) (§427). THIS SESSION discharged the 3 loaded
tasks and generalized:
- **§466 (Task 2, NULL):** scanned the whole U(2) extension family grading by pearsonr(F_θ,S) — max|r|=0.146
  over 148 extensions (bare 0.118). NO extension carries arithmetic S(T). operators/u2_ST_scan.py.
- **§467 (Task 3) + §467b/c/d:** the STRUCTURAL reason — deficiency (2,2) ⇒ ‖F_θ−F_θ'‖_∞≤2 (rank-2 bound);
  S(T) unbounded (Selberg/Montgomery); F_0=O(1) (1-D Sturm–Liouville Prüfer bound) ⇒ no extension reproduces
  S(T). Height-robust (r=0.075/0.021/0.029 at γ<95/120/150, vanishing). THREE instruments agree
  (u2_ST_scan, task1_bare_Xrobust, prov/arb_ST.py w/ mpmath arg-ζ = verdict B(SIZE-ONLY)). RH_TASK3_rank2_vs_ST.md.
- **§468 (generalization):** finite-deficiency self-adjoint-extension freedom is a bounded rank-≤n counting
  perturbation and CANNOT create the unbounded arithmetic S(T) — arithmetic content is a BULK invariant, not a
  boundary DoF. Rules out boundary-condition Hilbert–Pólya on arithmetic-free archimedean operators.
  RH_NOGO_extension_arithmetic.md.

**CONCLUSION / REDIRECT:** the pure-archimedean prolate route is honestly closed as a Hilbert–Pólya
realization. S(T)'s prime content must enter through the BULK = the CCM SEMILOCAL measure
dm_S=|∏_v L_v(½−is)|²ds (each prime p adds |L_p|²). This is the live continuation.

**NEXT BUILD (loaded, constructive):** first-order |L_p|² measure-perturbation of the base prolate
eigenfunctions — does S={p} inject a **log-p-periodic** level modulation = S_p(T)=−(1/π)Σ_k p^{−k/2}sin(kT log p)/k
that the archimedean operator entirely lacks? If yes → affirmative evidence the semilocal/bulk route supplies
the arithmetic (the constructive counterpart to the §466–468 no-go).
PREREQ (search-before-build, do NOT skip): re-read CCM 2023 (arXiv:2310.18423) for the EXACT L²(ds) semilocal
operator definition, the base (S=∅) measure, and what the finite-S spectrum represents — the archimedean
operator took 2 false starts (§462/463) before reading the full paper; do not repeat that. Then perturb.

---
## UPDATE 2026-07-30 (cont.) — §469/469b/469c: primes located in the bulk; the SHARP next attack specified
CONSTRUCTIVE progress past the §466-468 no-go:
- **§469/469b (Handle, verified):** built the CCM semilocal measure's Jacobi matrix (operators/semilocal_jacobi.py).
  Base dm_0=|Γ(1/4+is/2)|² reproduces b_n=√((n+½)(n+1)) to 1e-9 (arithmetic-free). Adding prime p imprints
  PURELY off-diagonal (Δdiag~1e-14), size ~ p^{−1/2} (the arithmetic weight; RMS(Δb)/p^{-1/2}≈const, Δb_1≈−0.5p^{-1/2}),
  oscillatory in n. Grid-converged to 5e-7. The primes ARE in the bulk, exactly where the no-go said.
- **§469c (Report, refines §465):** the truncated matrix W=−J²+4π²(4N+1)−¼ (λ=√2, validated J) CONVERGES in K
  (17.13/27.98/36.78 at K=80/160/240 — §465's "M-dependent" was a setup artifact) BUT to the WRONG (Friedrichs)
  extension — those γ ≠ zeros (14.13/21.02/25.01). Prime shift K-stable but O(1)-large (not the physical correction).

### THE SHARP NEXT ATTACK (specified, promising, tractable) — index-space recessive-solution shooting
Diagnosis: the pentadiagonal W DECOUPLES into even/odd index sublattices (±s symmetry) → each is a clean
**3-term recurrence** in steps of 2:  −[b_{n-1}b_n v_{n-2} + (b_n²+b_{n+1}²)v_n + b_{n+1}b_{n+2}v_{n+2}] + d_n v_n = α v_n,
d_n=4π²(4n+1)−¼.  Hard truncation (v_K=0, Friedrichs) gives the WRONG extension (γ=17.1). The zero-realizing
extension = imposing the **RECESSIVE (L²/minimal) solution as n→∞** (the index-space analog of the ODE's
pure-sin-at-∞ BC that gave §464's ~3%). METHOD: for trial α, shoot the 3-term recurrence outward, quantize α by
requiring the minimal solution (continued-fraction / Pincherle, or match to the asymptotic recessive form of the
recurrence). VALIDATE on the BASE first — it MUST reproduce 14.05/20.78/26.11 (§464) from the matrix recurrence,
not 17.1. If it does → apply the SAME recessive-BC shooting to the SEMILOCAL 3-term recurrence (which has NO ODE,
so this is the only route) → get the semilocal spectrum with the correct extension → test its fluctuation vs S_p(T).
This is the concrete path to the spectral realization the delicate extension has been blocking. Intricate (get the
recessive/Pincherle quantization right — do NOT rush into a §462-style false claim; validate against base 14.1 as
the hard gate before trusting any semilocal number).
