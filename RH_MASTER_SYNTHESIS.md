# RH — Master Synthesis (2026-07-30/31 session arc)

Honest, durable account of what this session established, with what confidence, and what the frontier is.
Report/Handle labels per Charisma. Every "certified" statement points at a proof artifact.

---

## 1. The one deliverable that is a genuine theorem: a de Bruijn–Newman bound past the record

**Λ ≤ 0.177482 — CERTIFIED** (Handle, conditional as stated below).
- Proof artifact: `~/cobo/roundtable/room/prov/tao_assembly_17748_FIRED.out`; assembly
  `tao_assembly_17748_fire.py` parses+verifies all four gates from the flint-Arb artifacts.
- Gates: **L** exact rational (0.165 + 0.158²/2 = 177482/10⁶); **C** large-N forall lemma (A+B at N*=10⁸
  = 1.766 < 2); **B** 17 flint-Arb F_lo tiles covering [690988, 10⁸] (binding tile 8× threshold, deep
  tiles 0.69–0.85); **A** region-(a) covering cert (worst edge margin +0.04077 > 0).
- Hypothesis: zeros of H_{0.165} confined to |Im| ≤ 0.158 at X = 2T, T = 3×10¹² (Platt–Trudgian, published).
- **Ladder this session:** published record Λ ≤ 0.22 (Polymath15, arXiv:1904.12438) → in-house 0.19862
  (prior roundtable) → **0.18280 (§483)** → **0.177482 (§489)**.
- **This is the ladder's ENDPOINT with the inherited barrier** (§488): the binding region-b F_lo goes
  negative just below Λ ≈ 0.176 at the Polymath15 mollifier exponent (0.2).
- **BUT the floor is a CERTIFICATION artifact, not a hard wall** (§490/§490b/§491): (i) the true mollified
  point-margin |M(s*)| at N0 stays strongly positive (0.70–0.96) far below the floor, down to Λ≈0.05 — so
  positivity is NOT the problem at the barrier point; the difficulty is localized to the CONTINUUM/rectangle
  worst-case. (ii) The region-b mollifier exponent has large numerical slack: re-optimizing it (EXP 0.2→0.03)
  recovers F_lo>0 at Λ=0.1746 (where the inherited exponent gives −0.076), with a legitimate interior
  optimum. So a tighter Λ (~0.174 or below) is reachable by continuum-mollifier optimization.
- **VALIDITY GATE on that improvement:** the region-b bound is "EXACT from arXiv:1904.12438 sec 8.5" at
  (exponent, fudge)=(0.2, 1.03); the fudge co-occurs with the exponent and likely bounds an exponent-tied
  Taylor remainder. Certifying a re-optimized exponent needs re-deriving fudge(EXP) from the paper's b-bound
  lemma — not cleanly fetchable this session (WebFetch/ar5iv weak; crawl4ai needs browser install), and a
  self-contained direct continuum bound is infeasible (the mollifier exists precisely to avoid tiling the
  10¹²-scale oscillation). So: characterized completely, honest lead, NOT a claimed theorem past 0.177482.
- **Honest provenance:** my session contribution = completing gates A/B at the tighter targets + writing
  the verifying assemblies. The barrier construction (regions a/b/c, the large-N lemma, the H_t machinery)
  is INHERITED from a prior roundtable and NOT independently re-derived. Conditional on (i) Platt–Trudgian
  and (ii) the correctness of that inherited Polymath15-style construction. Not an independent re-proof.

---

## 2. The unifying result: what the wall IS (the session's deepest honest statement)

RH ⟺ the Weil quadratic form is positive. This session mapped that positivity **five ways** — KLMN form
bound, Weil form min-eig, Keiper–Li λ_n, the Dirichlet-jump decomposition, off-line-zero detection — and
they ALL collapse to the **same a = 1 zero-margin razor** (§470–474). That is not five results; it is one
**methodology theorem**: *any positivity re-encoding of RH terminates at zero slack, because the positivity
IS the theorem.*

The **function-field Rosetta** (§482–487) explains WHY, by instrumenting the SOLVED case (RH-for-curves,
Weil 1948):
- RH-for-curves is proven by the **Hodge index on the surface C×C**: the intersection form has signature
  (1, ρ−1), the primitive part is negative-definite, giving |Σα_i^r| ≤ 2g·q^{r/2} for all r, which FORCES
  |α_i| = √q (demonstrated by counterfactual: an off-circle zero breaks the bound, §484).
- Cohomologically this is **Rosati positivity**: the norm form N(m+nα) = m² + a_p·mn + p·n² on End⁰(E) is
  positive-definite ⟺ disc a_p²−4p < 0 ⟺ Hasse ⟺ RH-for-E. **RH-for-curves is literally "a discriminant
  stays negative,"** guaranteed by geometry (§487).
- **THE GAP, precise:** the number field has no surface `Spec ℤ ×_{F₁} Spec ℤ` and no finite endomorphism
  algebra with a norm form — so the identical Weil positivity has **no geometric source**. Our a = 1 razor
  IS the number-field Castelnuovo–Severi bound, a theorem in the function field, unproven in the number
  field precisely because the arithmetic surface that would supply it does not (yet) exist.
- Independent lit-honest placement: this razor edge is the Weil-form incarnation of **Rodgers–Tao Λ ≥ 0**
  ("RH, if true, is barely true"). Our contribution is the explicit five-lane map + the Rosetta dictionary,
  not new mathematics on the razor itself.

**So every lane bottomed out at the same place because it is the same missing surface.** This answers the
session's own §451–457 question ("only a global mechanism can work — what does one look like?"): the global
mechanism is intersection-form positivity of a surface.

---

## 3. The wall's full mechanical anatomy (number-field side, §476/§485)

At bounded test-function support, positivity climbs a **support ladder** (Connes–Consani proved the
archimedean rung and the p=2 rung [log2, log3); [log3, …) is the first OPEN rung). This session located the
obstruction mechanically:
- **§485 (symbol):** the Weil form is a multiplier form ∫F(r)²W(r)dr; its symbol
  W(r) = Ω(r)/2π − Σ_{log n<4W} Λ(n)/√n·cos(r log n)/π goes NEGATIVE exactly when p=3 enters (min −0.170
  at r ≈ 17.35), while it is razor-positive at the p=2 rung (min +0.012). The first open rung is the symbol
  sign-change — where archimedean density stops pointwise-dominating the finite primes.
- **§485b/c (compensation):** band-limiting smooths the sub-resolution negative dips (form stays PSD, §476c
  certified), but the accumulating negative area + the **prolate eigenvalue plunge** collapse the margin
  exponentially (~exp(−10·4W)). Honest correction (§485c): the collapse rate is prolate-driven, NOT the
  symbol-area growth (rates differ 5.6×).
- Net: number-field positivity holds at each rung but with **exponentially-thinning margin**; local
  arch-domination (CC's method) reaches exactly log2 and no further; past it the **global mechanism** (the
  missing surface) is required.

---

## 4. Instruments built (reusable), and lanes explored & closed

- **Nyman–Beurling lane (§477–481):** validated closed-form Gram `G[u,v]=A(u,v)/u` (Vasyunin) + closed-form
  `F_n=(γ−1−log n)/n` ⇒ d_N with NO ζ integrals (likely first such numerics). d_N²·log N → C_BD = Σ 1/|ρ|²
  (oscillatory, RH-consistent). Structure: G = D^{−1/2} H D^{−1/2}, H multiplicative-Toeplitz with symbol
  |ζ(½+it)|²/(¼+t²) — the known NB–Burnol framework, re-derived + instrumented. V↔Bettin-Conrey bridge
  `c₀(n/m) = −V(n⁻¹ mod m, m)` derived+verified. Files: `nb/nb_gram.py, nb_dist.py, nb_recip.py`.
- **Function-field Rosetta (§482–487):** `ff/ff_rosetta.py, hodge_mechanism.py, rosati_positivity.py`.
- **Weil-form symbol/rung tools:** `operators/rung_log3_margin.py, arb_rung_certA2.py, arb_zerofree_cert3.py`
  (identity-route Arb certification: full form PD ∧ no-3 form indefinite at the open rung, §476c/g).
- **CCM prolate operator (earlier, §461–469):** built the concrete Θ; showed arithmetic S(T) cannot come
  from a self-adjoint extension (rank-≤2 bounded vs unbounded S(T)) but IS in the semilocal bulk measure
  with exact Euler weights. `operators/u2_ST_scan.py, semilocal_jacobi.py, …`.

---

## 5. The frontier, and the honest forward options

The number-field proof needs a **positivity with a geometric source** — the arithmetic surface
`Spec ℤ ×_{F₁} Spec ℤ` and its Hodge index (Connes–Consani F₁ / Deninger dynamical cohomology). That is a
research program, not a session. Concrete alternatives:
1. **Higher RH-verification height** — the only lever to push Λ below ~0.176. Real, well-defined, large
   compute, incremental.
2. **Better barrier** — Polymath15's own open problem (stalled since 2019).
3. **The F₁ construction** — the real long game.

No Handle for RH anywhere. What stands: one certified conditional theorem past the published record, a
complete honest map of the positivity program, and the frontier named as sharply as it can be without new
mathematics.
