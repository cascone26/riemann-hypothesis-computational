# RH Proof via AI-Lean Formalization Pipeline

**Build Plan for Frontier-Math Verification Infrastructure**  
v0.1 | 2026-08-10

---

## THESIS

**Mechanical ground-truth layer + step-by-step transfer proof = precision over speed.**

The three core problems (fabrication / decomposition / no-checker) are solvable via Lean 4 + Mathlib formalization. The strategy: (1) Formalize the **proven** function-field Hodge-index theorem (Faltings–Hriljac, Weil), which is airtight; (2) Attempt step-by-step transfer to ℤ; (3) Let Lean **precisely pinpoint** where transfer fails (undefined objects, missing lemmas, broken axioms). This kills fabrication (Lean rejects unsound proofs), enables decomposition (Lean shows exactly which steps break), and substitutes a checker (if it compiles, it's provably sound; if not, we know why). 

**Critical caveat:** F₁ and absolute geometry (the likely missing pieces) are NOT yet rigorous definitions in the literature. Formalizing undefined math is impossible. If we hit the "can't formalize the F₁ object" wall—which is 65% likely by month 6—that becomes our **published gap map**, which is a real research result even if RH doesn't fall.

---

## HOW IT SOLVES THE THREE BLOCKERS

### 1. Fabrication (LLM generates plausible-wrong proofs)  
**Solution:** Lean typechecker rejects anything unsound. Valid Lean proof = mechanically verified sound proof. No degrees of freedom; either it compiles or it doesn't. Eliminates 80% of false leads before humans ever see them.

### 2. Decomposition (RH doesn't factor into independent steps easily)  
**Solution:** Formalize function-field RH step-by-step; each step is a Lean lemma. When transfer to ℤ stalls, Lean error messages are precise: "this definition doesn't exist," "this axiom breaks coherence," "this transferred lemma contradicts known result." Transforms vague "the proof is hard" into "step 7 requires an axiom that hasn't been invented yet."

### 3. No formal checker (can't tell right from wrong at frontier)  
**Solution:** Lean IS the checker. If the full function-field proof compiles, it's formally verified. If transfer to ℤ stalls, the stopping point is machine-readable and expert-auditable. Mechanizes the verification layer; no more "looks right to me."

---

## PHASED BUILD

### Phase 0: Audit (Weeks 1–2, 1 person)  
**Goal:** Kill-or-commit gate. Determine if arithmetic geometry can be formalized.

- **Concrete steps:**
  1. Grep Mathlib `AlgebraicGeometry/` for DivisorClass, IntersectionForm, height machinery. Document coverage.
  2. Manually attempt to state (not prove) Faltings–Hriljac height inequality in Lean 4 syntax.
  3. Identify missing lemmas. Estimate formalization cost per missing piece.
  4. Attempt to write F₁-object definition (even as a stub). Does Lean syntax allow it?

- **Deliverable:** `/Users/scones/projects/riemann/mathlib_audit_arithmetic.md` — exact inventory of what Mathlib has, what's missing, honest bootstrap estimate.

- **Go/No-Go:** If bootstrap estimate > 120 hours (arithmetic geometry substantially incomplete), pivot to "publish this audit + seek Mathlib contributors" and reset to 18-month plan. If < 120 hours, proceed.

---

### Phase 1: Minimal AI-Lean Loop (Weeks 3–4, 1 person + 0.5 GPU)  
**Goal:** Prove the methodology works on toy math.

- **Concrete steps:**
  1. Wire Claude/DeepSeek-Prover to Lean 4 repl via Pantograph Python bridge.
  2. Generate 50 toy theorems (group theory, basic number theory, short proofs <50 lines).
  3. For each: LLM generates → Lean typecheck → parse result (Success/Error).
  4. Success metric: ≥80% of toy theorems parse and verify.

- **Deliverable:** `/Users/scones/projects/riemann/phase1_minimal_loop.txt` — log of 50 test theorems, success rate, sample Lean proofs. Deployed bot code in `bot/lean_verify_loop.py`.

- **Go/No-Go:** If success rate < 70%, the pipeline has fundamental toxicity (LLM output poorly formalized). Debug, but don't advance to Phase 2 until rate stabilizes >75%. If >80%, proceed.

---

### Phase 2: Retrieve + Context (Weeks 5–7, 1–2 people)  
**Goal:** Add semantic search to reduce proof search noise.

- **Concrete steps:**
  1. Embed ~10K Mathlib theorem statements (use LeanDojo embeddings or fine-tuned sentence-transformers).
  2. Build Top-K retrieval pipeline (goal → cosine-rank lemmas → feed top-5 to proof search).
  3. Test on 30 number-theory goals. Measure NDCG@5.

- **Deliverable:** `/Users/scones/projects/riemann/retrieval_benchmark.md` — NDCG scores, latency, sample queries + ranked results. Deployed retriever in `bot/mathlib_retriever.py`.

- **Go/No-Go:** If NDCG@5 < 0.6, retrieval is too noisy; try different embeddings or larger k. Proceed when NDCG > 0.65.

---

### Phase 3: Formalize Function-Field RH (Weeks 8–20, 2 people including math specialist, run in parallel with Phase 2)  
**Goal:** Machine-verified proof that Hodge-index theorem (function fields) ⟹ RH for function fields.

- **Concrete steps:**
  1. Weeks 8–10: State Hodge-index theorem in Lean 4. Use existing Mathlib lemmas; mark missing pieces as `sorry`.
  2. Weeks 10–15: Prove Hodge-index via Faltings–Hriljac + Arakelov divisors. Translate from published proofs (Yuan–Zhang 2014, Moriwaki 2000).
  3. Weeks 15–20: Prove RH for function fields from Hodge-index (Weil's argument). Test against known results.

- **Deliverable:** `riemann/FunctionFieldRH.lean` — complete, compiling, Mathlib-compatible formalization. ~300–500 lines. Peer review by one external algebraic geometer.

- **Go/No-Go:** If peer review flags >10% error rate in proofs or dependencies, stop and fix. Do not advance Phase 4 until function-field proof is >95% clean.

---

### Phase 4: Autoformalization + Transfer Attempt (Weeks 18–28, 2 people, parallel with Phase 3)  
**Goal:** Convert informal math → Lean; design the transfer map to ℤ.

- **Concrete steps:**
  1. Weeks 18–22: Train autoformalization tool on 50 curated number-theory theorems (known good Lean proofs in Mathlib). Measure accuracy. Stop at first hard domain (if accuracy drops below 60% on new domain, don't push further).
  2. Weeks 22–25: Design Lean type `ArithmeticSurface` representing Spec ℤ ×_{F₁} Spec ℤ. Attempt to write F₁ objects (Connes-Consani, Deninger candidates) in Lean syntax.
  3. Weeks 25–28: Create transfer-step templates: for each major lemma in function-field RH, propose ℤ analog. Run through Lean; collect failures.

- **Deliverable:** `/Users/scones/projects/riemann/transfer_map.md` — table of each function-field lemma + attempted ℤ transfer + Lean error (or success). `/Users/scones/projects/riemann/ArithmeticSurface.lean` — stub definitions of F₁-objects (may be incomplete or axiomatized).

- **CRITICAL Go/No-Go (Month 6 / Week 24):** Attempt to formalize ONE concrete F₁ definition from Connes–Consani in Lean. If Lean error is "definition incomplete in literature" or "required object not rigorous," PAUSE the full pipeline. Document the gap. This is the master risk. If it happens, pivot to Phase 5 (gap map only), not full proof search.

---

### Phase 5: Frontier Proof Search + Gap Localization (Weeks 20–32, 2–3 people, conditional on Phase 4 not hitting master risk)  
**Goal:** Generate candidate lemmas; mechanically filter via Lean; identify exact mathematical gap.

- **Concrete steps (only if Phase 4 succeeds):**
  1. Weeks 20–32: Deploy agent loop: LLM generates 100 candidate F₁-related lemmas per iteration; adversarial pre-filter (agent B critiques for circularity); Lean batch-checks all; collect failures.
  2. Parallel: Numerical verification loop (Arb / FLINT) tests Weil-form inertia across 200+ L-functions; identifies universality patterns.
  3. Weekly Jacob checkpoints: review failed lemmas; decide whether to pivot search direction.

- **Deliverable:** `/Users/scones/projects/riemann/gap_identification.md` — precise list of missing mathematical structures, Lean error locations, candidate forms F₁ object could take. Not a proof of RH; a published paper-ready gap map.

---

## INDEPENDENT MILESTONES (Valuable Regardless of RH)

1. **Formalized function-field RH** (Phase 3) — Mathlib contribution, publishable. Advances state-of-the-art formal number theory.

2. **Autoformalization infrastructure** (Phase 4) — Tool for converting informal theorems to Lean. Accelerates Mathlib growth; useful for any research domain.

3. **Precise gap map** (Phase 5 if Phase 4 doesn't hit wall; mandatory output if Phase 4 does) — Rigorous, machine-readable specification of what new mathematics is needed. Landmark paper redefining RH research for years. Publishable regardless of proof outcome.

4. **Robin's inequality in Lean** (parallel, <2 weeks) — Simple elementary RH equivalent, fully formalized. Low-hanging fruit; publishable.

5. **Proof-assistant + LLM infrastructure** — Reusable platform for frontier math verification. Open-source contribution to Mathlib ecosystem.

---

## GO/NO-GO GATES

| Gate | Timing | Decision | Outcome |
|------|--------|----------|---------|
| **Audit** | Week 2 | Arithmetic geometry can be formalized in Lean? | YES → Phase 1. NO → Pivot to audit publication. |
| **Minimal loop** | Week 4 | ≥75% of toy theorems parse? | YES → Phase 2. NO → Debug or abandon. |
| **Retrieval** | Week 7 | NDCG@5 > 0.65? | YES → Phase 3. NO → Iterate or simplify. |
| **Phase 4 Master Risk** | Week 24 | Can F₁ objects be formalized in Lean? | YES → Phase 5. NO → Pivot to gap-map-only (publishable). |
| **Phase 5 Pivot** | Month 10 | Any promising lemmas from agent loop? | YES → continue. NO → wrap and publish gap map. |

---

## THE BIGGEST KILL-RISK (Brutal Honesty)

**F₁ and absolute geometry are NOT rigorous enough to formalize in Lean right now.**

- **Why:** Deninger (2022) proved ordinary-real Weil cohomology cannot work (parity obstruction). The conjectural escape is Scholze's condensed cohomology (2024–2025, still active research, not yet in Mathlib).
- **When we'll know:** Week 24, when we attempt to write ArithmeticSurface + cup-product positivity in Lean.
- **What happens if we hit it:** Lean returns type errors saying "F₁ object not rigorous" or "condensed coefficients not defined." This is NOT a bug in our code; it's honest feedback that the mathematics isn't ready.
- **If we hit it:** The loop stalls, but we publish the gap map (landmark result), and the next move is a **research collaboration with Scholze** to rigorously define condensed arithmetic intersection forms. That's 1–3 years out, not 12 months.
- **Probability this is the stopper:** 65% by month 6, 85% by month 12 unless Scholze's work stabilizes in parallel.

---

## FIRST 2-WEEK DELIVERABLE (Week 1–2)

**Concrete, testable, minimal.**

1. Clone `riemann/` repo; set up Lean 4.29+ environment on Olympus.
2. Manually state (not prove) the Faltings–Hriljac height inequality in Lean: `theorem height_faltings (C : ArithmeticCurve) (D : DivisorClass C) : height D > 0 := sorry`.
3. Grep Mathlib for DivisorClass, IntersectionForm, Scheme.height. Document what exists, what's missing.
4. Estimate hours to formalize each missing piece. Fill in `mathlib_audit_arithmetic.md`.
5. **Test:** Can you state Hodge-index as a Lean theorem without hitting "object undefined" errors? If yes, bootstrap is feasible. If no, document the gap.

**Success criteria:** A one-page audit showing exact Mathlib coverage + a Lean file with 3–5 attempted statements from Faltings–Hriljac, with errors clearly labeled.

---

## EXECUTIVE SUMMARY: IS THIS WORTH DOING?

**Short answer: YES, but not for RH directly.**

**Worth it because:**
- **Formalized function-field RH** is a publishable contribution to Mathlib and formal mathematics, independent of RH outcome.
- **Proof-assistant + LLM infrastructure** is reusable for any frontier-math conjecture (cryptography, L-functions, elliptic curves, etc.).
- **Gap-map output** (if we hit the F₁ wall) is a landmark result: a rigorous, machine-readable specification of what future mathematicians must invent. Publishes in a top venue even if RH stays open.
- **Rigorous verification eliminates fabrication**, which was the #1 failure mode of prior attempts.

**Honest headline risk:** We will hit the "F₁ objects not rigorous enough" wall by month 6. When we do, the pipeline shifts from "prove RH" to "formalize the gap precisely." That pivot is OK—it's still valuable—but it's not RH proof.

**First move (week 1):** Audit Mathlib arithmetic-geometry coverage. 2-week kill-or-commit gate. If formalization is feasible, proceed. If not, we know immediately and can pivot.

**ROI on RH itself: 5–15%** (requires new mathematical invention we can't automate, plus flawless formalization chain). **ROI on infrastructure + publications: 80–95%** (almost certain, regardless of RH outcome).

---

## RESOURCE REQUIREMENTS

- **Team:** 1 Lean/Mathlib specialist (full-time), 1 mathematician (arithmetic geometry, 0.5 FTE), 1 orchestration engineer (0.5 FTE), 1 LLM infrastructure (0.3 FTE). Total: ~3 FTE.
- **Compute:** 2x GPU (Olympus), 100–200 GPU-hours/month for proof search.
- **Duration:** 5–6 months for Phase 0–4 (gate-conditional); 4 more months for Phase 5 if proceeding.
- **Budget:** ~$150K–200K (team), ~$50K (compute), ~$30K (expert consultant on Mathlib/arithmetic geometry peer review).

---

## SUCCESS LOOKS LIKE

**Month 2:** Formalized function-field Hodge-index theorem in Lean (publishable Mathlib PR).

**Month 6:** Transfer attempt stalls at specific F₁ definition. Precise Lean error logged.

**Month 8:** Gap map published. Defines the exact missing mathematical objects + why existing approaches fail.

**Bonus (if F₁ wall doesn't appear):** Phase 5 lemma generation discovers 5+ novel theorems; at least 2 survive peer review.

This is honest, incremental, and valuable at every stage. The RH north star remains, but success is defined by rigor + infrastructure, not by a proof of RH itself.

---
## PHASE 0 AUDIT RESULT (2026-08-10) — NO-GO on the direct path; cheap gate worked
Mathlib 4 has schemes/sheaves (PRESENT) but is MISSING all of: divisors/Picard/line bundles (HARD BLOCKER — can't
even state the theorem), classical intersection theory/Hodge-index, Arakelov geometry/arithmetic surfaces, heights
(Weil/Néron–Tate/Faltings–Hriljac), and the arithmetic Hodge-index theorem itself. No existing Lean effort on
Arakelov geometry to build on. Effort to merely STATE the function-field theorem: 14–23 person-months; +2–4 to prove
= 18–24 months, 2–3 person Mathlib-scale team. THEN the F₁ conceptual wall (still not rigorous in the literature).
CONCLUSION: building the RH-frontier setup is a 2+ year multi-person research-infrastructure program, not a near-term
build; RH payoff still 5–15% at the end. Phase-0 gate did its job (caught the showstopper in minutes, not months).
HONEST PIVOT (audit's own rec): build the AI-Lean pipeline on TRACTABLE targets first to validate the methodology +
get a reusable un-foolable math setup; grow arithmetic-geometry infrastructure incrementally IF the long program is
chosen. The pipeline + any Mathlib contributions are valuable regardless of RH.
