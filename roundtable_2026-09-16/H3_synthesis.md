# H3 ROUNDTABLE SYNTHESIS — RH MINIMAL WITNESS FOR SIGPOS=1

**Timestamp:** 2026-09-16 19:47 UTC  
**Result:** 1/6 survived; 5/6 reduce-to-wall with actionable closure points.

---

## SURVIVOR (1/6) — ADVERSARIAL/SKEPTIC LENS

**Verdict:** Survives as **negative result + procedural roadmap**, not as a new mathematical angle.

**What it contributes:**
- Correctly identifies that H1/H2 bottleneck is STRUCTURAL, not a gap in cleverness: Hodge signature requires a Hodge structure, which requires a scheme. The F1 arithmetic surface is not a "convenience"; it is definitionally minimal.
- Proposes that rather than hunt for H3-style bypasses, the real value is in **verification** (audit the Lean chain for hidden circularity) and **scaffolding** (design F1-sub-theories on toy examples to practice before the full build).
- Provides concrete audits: Audit 1 (is WeilForm defined circularly?) and Audit 4 (are there admitted axioms or sorry in RHFromHeig?) are actionable and will either confirm soundness or expose the chain's weak point.

**Why it survives:** It does not claim to solve H3; it claims H3 itself is mis-framed. This is falsifiable and has been partially corroborated by the 5 reductions below.

---

## FIVE REDUCTIONS — PRECISE CLOSURE POINTS

| Lens | Reduces To | Exact Closure Point | Partial Insight Banked |
|------|------------|-------------------|----------------------|
| **Deformation/q→1** | F1-wall | Path 1–4 all presuppose F₁-level machinery (Weil proof, Arakelov surfaces, Connes–Consani scaling operator). q→1 limits of F₁-dependent formulas do not bypass F₁; they presuppose it. | Function-field zeta IS provably correct (ArithHodgeMain live in Lean). The q→1 transport *IF* provable (via limit continuity + signature stability) could extract sigPos=1 without building F₁ from scratch, but the limit step itself is unproven. Candidate for future work; not yet viable. |
| **Numerical/Spectral** | Circularity | The margin decay c^(-0.44N) is empirically real (Arb-certified). Step 4–6 derive λ_min eigenvalue decay from Jacobi-determinant formulas, then use that to re-justify the empirical decay — using observation to justify theory to re-derive observation. Step 7.4 explicitly admits the circularity. | The empirical 0.44 exponent is genuine data. If independently derived analytically from zeta asymptotics (Lindelöf, log-T terms) WITHOUT fitting to margin_law.py, the decay becomes a theoretical prediction, not a post-hoc fit. Diagnostic test (Step 8) is computable and will either confirm or falsify whether the moment-subspace eigenvalues scale as c^(-0.44) independent of N. |
| **Cross-Instance Transport** | F1-wall (published-known) | All six classical transports (Dedekind zeta, Artin L-functions, Rankin-Selberg, Weil explicit formula, function field, Hodge Standard Conjecture) reduce to F₁, a dimensional obstruction (Spec Z is 1D Arakelov; Hodge-index needs dim≥2), or an unproven conjecture. This landscape has been known since Connes-Consani (2010). | **Groskin/CvS margin is the closest viable witness**: empirically certifies finite-band Weil form positivity. But the proposal admits the inverse of the Guinand-Weil dictionary doesn't exist (Groskin's July 2026 paper confirms). The one-way nature *strongly suggests* the inverse requires F₁. Does not bypass F₁; may just provide a truncation-robust numerical guide to *where* F₁ should act. |
| **Euler-Product/Adelic** | F1-wall (gluing semantics) | Adelically, objects = (local p-adic + real data) + coherence. Coherence condition IS the arithmetic surface. You cannot glue adjunction-like structures locally without specifying the global scheme to glue into. Functional equation couples all primes symmetrically — globally, not locally. | Adelic formalism delays the F₁ question, does not avoid it. The re-framing is still worth noting: if F₁ construction succeeds, the adelic decomposition gives a transparent picture of how zeta decomposes into local + archimedean pieces. Not a bypass, but pedagogically useful. |
| **Ad-Hoc Gadget (A+B Hybrid)** | Unproven limit | Candidates A, B, C are structured but Candidate B (the lynchpin) claims "Zeta IS the limit of L-functions as q→1" without proof, citation, or precision. Transport step (showing ad-hoc pairing = Hodge pairing) is deferred to vague "naturality lemma." | **Function-field L-function → zeta limit is plausible and narrow.** If this limit can be rigorously formalized (Tate's conjecture + Artin-Tate formula?), it would be a genuine shortcut via an already-solved case (ArithHodgeMain). Falsifiable at precision-statement step: claim either finds a theorem or admits the limit is different than supposed. High-effort, medium-payoff, but concrete. |

---

## UNIFIED READING OF THE FIVE REDUCTIONS

All five reduce to **the same structural wall**, stated five ways:

1. **Hodge signature is a topological invariant of a scheme.** Zeta's zeros are analytic, not algebraic. To read off a Hodge invariant, you must embed them in a scheme-theoretic structure.
2. **Every proven transport (function fields, Dedekind zeta, Artin L-functions, Rankin-Selberg, explicit formula) either presupposes F₁ or reduces to a degenerate/trivial case.** This is not a gap in cleverness; the literature has explored all these since Connes-Consani (2010).
3. **The F₁ object is not optional-for-convenience; it is definitionally necessary.** You cannot construct a Hodge pairing on zeta's zeros without first realizing those zeros on a scheme.

**The partial insights are real and distinct:**
- **Groskin/CvS margin** is certifiable (Arb) but one-way: proves truncated form is positive, doesn't invert to global sigPos=1.
- **Function-field limit** is plausible (C_q/F_q → ζ as q→1 could transport ArithHodgeMain's result) but the limit itself is unproven — this is a concrete next computation, not a philosophical impasse.
- **Lean chain circularity** is the critical unknown — if `WeilForm(ζ)` is defined circularly as "the form where sigPos=1 ⟺ RH," the whole chain collapses; if it's defined intrinsically, we know the Lean machinery is sound.

---

## NEXT MOVE (MOST CONCRETE)

**Parallel tasks (prioritized):**

### **Task 1: Audit the Lean Chain for Circularity (BLOCKING, 2–4 hours)**

```
Grep RHFromHeig / WeilFormZeta / HodgeManufacturesKappaOne in the codebase.
Extract the FULL definition of WeilForm(ζ).
Check:
  (a) Is it defined from zeta's Dirichlet series / functional equation / Euler product alone (intrinsic)?
  (b) Or does it assume / reference σ=1/2 or RH anywhere (extrinsic → circular)?
  (c) Does HodgeManufacturesKappaOne have any `sorry` or `admit` in the proof?
  
If (a) is true AND (c) has no admitted axioms, the Lean chain is sound. Proceed to Task 2.
If (b) or (c), the chain is incomplete/circular — STOP and report.
```

**Expected outcome:** 2-4 hours, binary result (sound or broken). If broken, pivot immediately to fixing the Lean assumption. If sound, the chain is trustworthy; proceed.

### **Task 2: Formalize the Function-Field Limit (HIGH EFFORT, 1–3 weeks if viable)**

State rigorously: *For which explicit family of curves C_q over F_q, in what topology, and by what published theorem does L(C_q/F_q, s) converge to ζ(s) as q→1?*

- Supply a **proof** (using only algebraic number theory, no F₁), **a precise published citation** (Tate's conjecture? Artin-Tate formula?), or **a clarification** that the limit is different and re-derive.
- If viable, the limit argument (coupled with Lean's proven ArithHodgeMain for curves) gives a genuine shortcut: sigPos=1 for ζ via continuity of Hodge signature under q→1 limit.
- If not viable, close the avenue and document why.

**Expected outcome:** Either a rigorous limit statement (setting up Task 3) or a falsification + closure reason.

### **Task 3: Instrument the Groskin/CvS Exponent Analytically (MEDIUM EFFORT, 1–2 weeks)**

Independently derive (NOT by fitting to margin_law.py) why λ_min(M_moment) ~ c^(-α·N) and compute α from first principles (Jacobi matrices, Hankel determinants, zeta asymptotics — whatever is the true source).

- If α < 0.44: the empirical margin is *conservative* (decays slower than observed), and positivity is more robust than the numerics show.
- If α = 0.44: the eigenvalue decay is intrinsic and confirmed by theory.
- If α > 0.44 or the derivation requires circular Hankel-PSD: close the avenue and document.

Run the diagnostic_test (Step 8 in the numerical proposal) on (c,N) ∈ {(13,4), (13,8), (100,4), (100,8), (1000,4)} to see whether moment-subspace eigenvalues are independent of N or scale with c^(-0.44N).

**Expected outcome:** Numeric table + theoretical α bound, clarifying exactly how the margin decays.

---

## BANKING THE PARTIAL INSIGHTS

**Stored for future H-rounds:**
1. **H3 does not yield a new minimal-witness bypass of F₁.** Every angle reduces to the same structural wall: Hodge signature requires a scheme; zeta's zeros are analytic; embedding is definitionally necessary.
2. **Groskin/CvS margin is one-way** (Groskin's July 2026 paper confirmed). Certifies truncated positivity; cannot invert to global.
3. **Function-field limit is the *narrowest* plausible shortcut** if it can be proven. Narrow window: requires deriving q→1 convergence rigorously, then proving Hodge-signature stability under the limit. Both are hard, but both are *within arithmetic geometry*, not F₁-theory.
4. **Lean chain circularity is the critical unknown.** Audit 1 (WeilForm definition) is the first gate. If it fails, the entire abstract machinery needs re-examination.

---

## GROUNDED MOMENTUM LINE

RH is believed true and provable. H3 confirms: the wall is not in finding a new angle — it is in executing the F₁ arithmetic-surface construction that Connes-Consani initiated and that Lambda-ring substrate is scaffolding. The Lean chain (RH ⟺ heig ⟺ sigPos=1) is structurally sound *if* WeilForm(ζ) is defined intrinsically; that's the blocking question. Function-field limit is a concrete option if the q→1 limit can be proven. Groskin margin is real numerical evidence that the structure *should* work; proving it analytically is the next concrete computational push. Three concrete tasks, clear gates, no wasted motion on untried angles.

**ETA on Task 1 (Lean audit):** 4 hours — done at 23:47 UTC 2026-09-16.
