# RH Barrier Map — proof-styles × properties-used × witness objects

**Compiled 2026-08-07 (barrier session, LEDGER §802).** Canonical form for every kill:
**[proof-style] + [exact properties of ζ it uses] + [witness object that has those properties but breaks RH].**
A proof-style is barred when a witness exists: no argument using ONLY those properties can prove RH,
because the same argument would "prove" it for the witness, which violates it. Model example:
Davenport–Heilbronn (functional equation, no Euler product, off-line zeros).

Companion docs: `OBSTRUCTION_MAP.md` (mechanism-level, routes 1–7), `RH_DIFFICULTY_CONSERVATION_MAP.md`
(framework-level reality/correspondence law). This map is the **witness-side dual**: OBSTRUCTION_MAP says
*how* each route fails on ζ itself; this map says *why no repair can work* within the same property set.

**Charisma labels:** Handle = verified on our own instrument (proof pointer in repo). Report = literature,
agent-sourced 2026-08-07; citations marked ⚠ have unverified dates/details — re-check before external use.

---

## Kill-type taxonomy

- **K1 witness** — an object with the used properties violates RH. The strongest kill: bars the whole style.
- **K2 premise-false** — the style's premise fails for ζ by direct computation.
- **K3 circularity** — the style's target is provably equivalent to RH; certificates reduce to the thing itself.
- **K4 logic-ceiling** — success would prove a strictly stronger open statement (so the tools aren't strong enough).
- **K5 identity** — the "test" is an identity satisfied regardless; carries no positivity information.

K2/K3/K5 need no witness (the kill is internal). K1 and K4 rows without a pinned witness are the **gaps** —
each gap is a theorem to attempt (list at bottom).

---

## Witness inventory

| W | Object | Has | Lacks | RH-breaking behavior | Tier |
|---|--------|-----|-------|---------------------|------|
| W1 | **Davenport–Heilbronn** F(s)=(1−iκ)/2·L(s,χ)+(1+iκ)/2·L(s,χ̄), χ mod 5, κ=0.28408 | Dirichlet series, continuation, **Riemann-type FE**, growth bounds | **Euler product** | Off-line zeros; infinitely many (D-H 1936). **In-repo: off-line zero verified** \|F(0.808517+85.699348i)\|=6.5×10⁻⁷ (LEDGER §725) | **Handle** (built here) + Report (1936) |
| W2 | **Epstein zetas**, class number h>1 (and linear combos of Euler products) | FE, continuation, growth | Euler product | Off-line zeros (Davenport–Heilbronn 1936; Voronin 1976; Bombieri–Hejhal, Duke Math. J. **⚠ vol 80 = 1995**, agent said 1976 — verify) | Report |
| W3 | **Diamond–Montgomery–Vorhauer** Beurling system (Math. Ann. 334 (2006) 1–36) | **Euler product** (g-primes), pole at s=1, N(x)=κx+O(x^θ) with θ<1, PNT-quality π | **FE** | Infinitely many zeros on σ=1−a/log t, i.e. RH-analogue fails badly despite good counting | Report |
| W4 | **Broucke–Debruyne–Vindas** (Adv. Math. 370 (2020), arXiv:2004.11501) | Euler product, pole, π(x)=Li(x)+O(√x), **RH-analogue HOLDS** | FE | Converse-side control: RH-analogue ⇏ smooth integers — N(x) has Ω±(x·e^{−c√(log x log log x)}) oscillation | Report |
| W5 | **Diamond–Zhang** pair (2007) | Euler product, N(x)=kx+O(x^{1/2}e^{c(log x)^{2/3}}) | FE | **Tunable**: system A all zeros on the line; system B zeros on σ=1−a/log t. Same counting quality, opposite RH verdicts | Report |
| W6 | **Broucke–Debruyne–Révész [α,β]-factory** (arXiv:2309.01567, Trans. AMS 2024) | Euler product, pole, **dialable** ψ error x^α / N error x^β; **prescribed zero set ρ** via ζ_P(s)=E(s)e^{Z(s)} (Thm 3.2) | FE | Witness FACTORY: unconditionally realizes any α∈[0,1), β∈[1/2,1) with zeros where you put them. Constraint (Hilberdink-type): max{α,β}≥1/2 always | Report |
| W7 | **Broucke–Debruyne** (Acta Arith. 207 (2023), arXiv:2211.08716) | Euler product, pole, **Axiom A** (N(x)=κx+O(x^θ), θ<1) ⇒ Carlson-type zero-density bounds APPLY (Révész 2209.01689 Thm 1 needs Axiom A alone) | FE | **Ω(T^{c(1−σ)}) zeros** in [σ,1]×[−iT,iT], σ arbitrarily close to 1 — saturates the density bound; Révész §1.2: Carlson-type estimates "best possible" | Report (closes B5's gap; found by lit-gate 2026-08-07) |

**The anti-witness fact (load-bearing):** there is **NO known object with Euler product + FE + off-line
zeros.** Every known RH-violator drops exactly one of the pair. This is the content of the Selberg-class
conjecture, and we hit it independently on the instrument: §725 found off-line-zeros and Euler-product-loss
inseparable in every available example; §431+§439 demonstrated the two components directly (off-line
placement breaks PSD; permuting log p breaks PSD) — duals through the explicit formula.

---

## The barrier table

| # | Proof-style | Exact properties of ζ used | Witness / kill | Status |
|---|-------------|---------------------------|----------------|--------|
| B1 | **FE + growth/order alone** (symmetry + complex analysis; naive de Branges axioms; "zeros must be symmetric so…") | FE, continuation, order 1, trivial zeros | **K1: W1 (D-H), W2 (Epstein)** — same properties, RH false | **DEAD-proven** |
| B2 | **Counting-precision alone** ("PNT-strength or better error terms force RH") | Euler product, pole, ψ/N error exponents | **K1: W3 (DMV), W5-B, W6** — EP + good counting at every unconditionally-achievable quality, zeros off the line | **DEAD-proven** |
| B3 | **Moment/Hankel** (Li λ_n a Stieltjes moment sequence) | numeric λ_n values | **K2**: Hankel fails PD at N=2 (det=−0.0037, in-repo); (λ_n) is log-concave, the opposite | **DEAD** (premise-false) |
| B4 | **Pointwise/local majorization** (Beurling-Selberg majorant of the prime comb) | local decay data of arch kernel vs prime comb | **K2**: prime comb exceeds arch kernel at every prime ≥5 (in-repo, OBSTRUCTION_MAP #2). Witness form open → **G3** | **DEAD as stated; witness gap** |
| B5 | **Large sieve / zero-density-strength L² information** | mean-value bounds on Dirichlet polynomials | **K1: W7** — Beurling system satisfying the exact hypotheses of the Carlson-type density theorems (Axiom A; Révész 2209.01689 Thm 1, no FE needed) with Ω(T^{c(1−σ)}) zeros up to σ→1. Zero-density-strength information provably cannot force RH (G1 resolved: already published, see BARRIER_G_RESULTS.md) | **DEAD-proven** (witness in literature) |
| B6 | **de Branges / canonical system via FE-structure only** | FE, Hermite–Biehler / phase structure | **K1: W1/W2** (D-H-type satisfies the axioms used); also K2 in-repo: A+B has 4 negative eigenvalues (not a dB inner product) | **DEAD without EP input** |
| B7 | **Positive-proportion / mollifier (Levinson→Conrey ≈41.7%)** | Dirichlet-polynomial mean values (coefficient-level, works without EP) | **K1: W1** — D-H has a **positive proportion on the line** (Balanzario–Sánchez-Ortiz, Acta Arith. **⚠ 2001/2007 — verify**) AND infinitely many off. Positive proportion can coexist with RH-failure ⇒ the machinery cannot close to 100% from these properties | **DEAD as full-RH route** (still the best partial axis) |
| B8 | **Local / per-prime / factorization certificates** | Euler factors individually | **K2**: ζ(½)<0 ⇒ per-prime factors not locally positive (A14, in-repo); metatheorem 2026-07-20 | **DEAD-proven** |
| B9 | **Single positive test function** ("find the magic g with A+B−P≥0 visibly") | one positive h in the explicit formula | **K5**: A+B−P=Σ_γ h(γ) is the explicit-formula **identity** — positive term-by-term iff zeros on-line; an off-line zero enters sign-indefinitely (§37, 30-digit in-repo) | **DEAD** (identity, no information) |
| B10 | **Family averaging** (χ-average; any almost-all argument) | orthogonality across a family | **K4 + K1 (continuous-level)**: a large sieve forcing EVERY χ would prove GRH (§62); explicit witness family BUILT — mean positivity + 0.61σ outlier violating, no variance bound can exclude it (`operators/barrier_g2_averaging_witness.py`, BARRIER_G_RESULTS.md G2) | **Barred; witness built (continuous level)** |
| B11 | **Deformation / heat flow (de Bruijn–Newman)** | Λ≥0 proven (Rodgers–Tao 2018); target Λ≤0 | **K3-adjacent**: heat flow **destroys the Euler product** for λ≠0 (§47d-e in-repo) — the deformation discards the one property no witness can break (see anti-witness fact). §34c: earlier DBN "route" retracted as circular | **Barred as EP-blind deformation** |
| B12 | **Operator/spectral certificates for τ⪰0** | the Weil form itself | **K3**: τ⪰0 ⟺ RH; every finite certificate (Gershgorin, SOS, …) reduces to Weil (OBSTRUCTION_MAP #5) | **Circular; finite-N results real but finite** |

---

## Surviving rows — what any proof MUST use (the point of the map)

Reading the table by complement, a successful proof must jointly use properties that **no witness object
shares**:

1. **Euler product AND functional equation together.** Every known RH-violator drops exactly one (W1/W2
   drop EP; W3–W6 drop FE). Their conjunction is exactly the Selberg-class conjecture's hypothesis, and
   it is the *only* property-pair with an empty witness column. A proof using either alone is dead (B1, B2).
2. **The specific arithmetic, not counting quality.** W5/W6 realize every achievable counting precision on
   both RH-true and RH-false systems — error exponents carry no RH information by themselves (B2, B5).
3. **Global positivity structure, not local/pointwise** (B4, B8, B9 — metatheorem 2026-07-20: positivity is
   irreducibly global/product-emergent).
4. **Individual-isolating structure, not averaging** (B10; §62/§59: the average→individual gap = the
   analytic→geometric gap; function fields isolate the individual via polarization).
5. **EP-preserving methods** — any deformation/relaxation that destroys the Euler product (B11) walks into
   witness territory where RH is simply false.

This is consistent with, and sharpens, `RH_DIFFICULTY_CONSERVATION_MAP.md`: "correspondence" = the EP/
arithmetic input (the half no witness can fake), "reality" = the global positivity (the half every
framework leaves hard).

---

## Gaps = theorems to attempt (one per missing witness)

- **G1 (large-sieve witness, feeds B5).** Build, via the W6 factory (prescribe the zero set ρ in
  ζ_P=E·e^Z), a Beurling system with Carlson/zero-density-strength estimates AND a zero at prescribed
  σ₀>1/2. Would prove: zero-density-strength information can never force RH. Beurling zero-density
  machinery already exists (arXiv:2209.01689, 2407.12746) — this looks buildable.
- **G2 (averaging-wall witness, feeds B10).** Construct a finite family of Beurling systems with mean
  Weil-positivity + large-sieve-strength variance bound, where one member violates its RH-analogue.
  Would upgrade §62's logic-ceiling to a hard witness barrier.
- **G3 (local-data witness, feeds B4/B8).** Two Beurling systems with identical per-prime local data
  (same g-prime multiset up to height X, X→∞ scaling) with opposite RH-analogue verdicts — W5 (Diamond–
  Zhang's tunable pair) is 80% of the way there; pin whether their local data can be matched.
- **G4 (the surviving-row converse).** Prove for a natural class of Beurling-type systems: FE + EP ⟹ no
  off-line zeros (a Selberg-class-conjecture fragment in the Beurling setting, where FE is the rare
  property). Even a narrow class would convert surviving-row 1 from empty-witness-column to theorem.

## Lit-gate verdict (rule 13, run 2026-08-07)

Individual barriers B1/B2/B6/B7's witnesses are **well-known** (D-H 1936 → Bombieri–Hejhal;
Kaczorowski–Perelli 1999 on Selberg-class structure; Conrey 2003 states "we do not know which properties
of ζ are essential to RH"). The Beurling witnesses W3–W6 are published 2005–2024. **Not found in the
literature:** a formalized axiom-lattice / barrier framework for RH in the complexity-theory sense
(relativization / natural-proofs analogue) — the *systematic* map, the K1–K5 taxonomy, and the G1–G4
witness-gap program are candidate-novel as organization. Claim no individual barrier as new.

---
*Living document, same discipline as OBSTRUCTION_MAP: add a row when a style is barred; never delete a
verified barrier; a gap row (G) graduates to a witness row when its object is built.*
