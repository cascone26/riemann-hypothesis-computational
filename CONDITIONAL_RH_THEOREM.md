# A Machine-Verified CRITERION for RH (as of 2026-08-30, cy362 — corrected after Roundtable skeptic review)

**Status: axiom-clean and real — but it is a formally-verified Hilbert–Pólya / Weil-positivity CRITERION,
NOT a reduction of RH to a strictly-weaker problem.** The "input" below is *equivalent to* RH (the heig /
Hilbert–Pólya form) or *plausibly stronger than* RH (the isometry / F₁-realization form). So this is a
verified REFORMULATION of RH, not a conditional proof that leaves only an easier gap.

> **⚠️ Correction (cy362):** an earlier framing called this a "complete conditional proof, reducible half done."
> That over-stated it. The Roundtable skeptic seat correctly flagged that the hypothesis is not weaker than RH.
> The formalization itself is genuine (the theorems exist on the live Olympus forge, conclude Mathlib's real
> `RiemannHypothesis`, and are axiom-clean — the skeptic's "theorems don't exist" claim was from grepping a
> stale mirror). What is real is the *criterion*; what was wrong was calling it a partial *proof*.

---

## The theorem we have (machine-verified, GREEN BOARD, axiom-clean)

> **RH ⟸ [ONE input].** If the ζ Weil explicit-formula quadratic form is *negative-definite on its
> primitive part* — equivalently, is isometric to the Faltings–Hriljac template `arithB NT` with a
> positive-definite height `NT` — then the Riemann Hypothesis holds.
>
> **But note (the honest part):** the Hilbert–Pólya form of this input, `heig` = {∃ symmetric operator T with
> `T·v = s(1−s)·v` at every nontrivial zero}, is **equivalent to RH** (under RH the diagonal operator works;
> if RH fails, an off-line zero forces a complex eigenvalue a symmetric operator cannot have). So the theorem
> `heig ⟹ RH` is a true, machine-verified Hilbert–Pólya CRITERION — a reformulation, not a reduction.

Every step from that input to Mathlib's `RiemannHypothesis` is formalized in Lean 4 + Mathlib on Olympus
(`~/rhforge/forge/`) and passes the axiom gate using only `{propext, Classical.choice, Quot.sound}`
(no `sorry`, no extra axioms). Verified GREEN BOARD 10/10 on 2026-08-30 (cy358): `RHForgeMain`,
`RHReductionMain`, `ArithHodgeTransfer`, `RHFromBoxLimitClean` (+ cy345: `RHFromHeig`, `AbstractRH`,
`HodgeManufacturesKappaOne`, `DiscriminantTraceBound`).

## The axiom-clean chain (input → RH)

1. **Input** (the one open hypothesis): the ζ Weil pairing `P` is isometric to `arithB NT`, `NT` positive-definite.
   - The hyperbolic degree-plane part (`F²=0, S²=0, F·S=1`) is FREE — baked into the `arithB` template.
   - The height-positivity (`NT` pos-def) is BUILT for the arithmetic-curve case (`ArakelovVecRHodge.pairingFun_self_pos`).
2. `ArithHodgeTransfer.arith_hodge_signature_transfer`: isometry to `arithB NT` (NT pos-def) ⟹ `sigPos P = 1`.
3. `HodgeManufacturesKappaOne` / `KappaLeOneSignature`: `sigPos = 1` ⟹ `#neg ≤ 1` (κ=1).
4. `AbstractRH.re_eq_half_of_strip_bound` / `RHFromHeig` / `RHForgeMain`: `#neg ≤ 1` (strip bound) ⟹ RH.

Proven *analog* instances (that the template is realizable): `EllipticArithHodgeMaster` §1083 — weak
Mordell–Weil ⟹ `sigPos=1` for elliptic curves (axiom-clean); `FrobeniusIsometryEigenvalues` §1728 —
`MᴴM=q·I ⟹ |eig|=√q` (the Weil function-field heart). The function-field |λ|=√q chain (§1720–1726) is complete.

## The one open input, in its equivalent forms (all = Weil positivity ⟺ RH)

- **Geometric:** the Hodge index theorem (intersection pairing negative-definite on the primitive part)
  for the **square of the Connes–Consani scaling site** — the F₁ "arithmetic surface" Spec ℤ ×_F₁ Spec ℤ,
  which *exists* in the C-C framework (ξ = its Hasse–Weil zeta; divisors = Schwartz kernels of the Frobenius
  correspondences).
- **Operator:** the CvS truncated-Weil-form even-sector ground eigenvalue stays positive in the limit
  `λ_∞^even > 0` (⟺ the finite-cutoff ground state's zeros converge to the ζ zeros as c→∞).
- **Classical:** Weil's positivity criterion for ζ (`A(g)+B(g) ≥ P(g)`), equivalent to RH.

## Honest status of the input

It is **Connes–Consani's own actively-researched open conjecture** (Connes 2026 §6.4 *heuristic*;
`weilpos`, `truncweil` 2605.20224). Backed by a heuristic predicting a **positive** limit
(`log₁₀ λ_∞^even ≈ −530`) and **113–334 orders of magnitude** of numerical convergence — RH-consistent,
but **unproven** (truncweil: "we do not claim to validate Connes' convergence conjecture").

## What this is / isn't (corrected — read with the top-of-file ⚠️ and §B below)

- **IS:** a genuine, novel, machine-verified **CRITERION** for RH — the entire chain from the input to Mathlib's
  `RiemannHypothesis` (positivity → signature → κ=1 → RH) is axiom-clean. A formally-verified Hilbert–Pólya /
  Weil-positivity reformulation, which no one else has built.
- **IS NOT:** a proof of RH, and **NOT a reduction to a strictly-easier problem.** The input is *equivalent to*
  RH (heig ⟺ RH) or *plausibly stronger* (the F₁-realization). So "our side is done" means the *downstream
  reformulation* is done — it does **not** mean an easier gap remains. The gap is RH itself, as Weil positivity.

## Next moves (require external input — all autonomous/reformulation routes exhausted; see §C)

1. Watch C-C `weilpos`/`truncweil` for a proof of the input (Weil positivity / scaling-site Hodge index) — the
   instant it lands, the forge closes RH, machine-verified.
2. A focused research push on the input itself (Connes–Moscovici prolate exact spectral realization; the
   deficiency-(4,4) self-adjoint extension) — sizable functional analysis, the real Hilbert–Pólya frontier.
3. (Infrastructure, not RH progress) a formalized interval-arithmetic layer in Mathlib → a fully-Lean finite
   Weil-positivity certificate.

_Source: THE_BRIDGE.md cy329–358; memory `rh-solver-os`; forge `~/rhforge/forge/` on Olympus._

---

# Full Verified Inventory & Attack-Surface Map (appended cy368, 2026-08-30)

## A. What is machine-verified (Lean 4 + Mathlib, axiom-clean {propext, Classical.choice, Quot.sound}, 0 sorries)

**On Olympus `~/rhforge/forge/` (the arithmetic-Hodge / spectral chain):**
- `riemannHypothesis_of_boxEigen_realization` (RHFromHeig): heig ⟹ Mathlib's genuine `RiemannHypothesis`. [gated GREEN]
- `arith_hodge_signature_transfer` (ArithHodgeTransfer §1025): isometry to `arithB NT` (NT pos-def) ⟹ sigPos=1.
- `HodgeManufacturesKappaOne`: sigPos=1 ⟹ #neg≤1. `AbstractRH.re_eq_half_of_strip_bound`: #neg≤1 ⟹ RH.
- `hodge_signature_of_adjunction` (HodgeAdjunction): F²=S²=0, F·S=1, neg-def on primitive ⟹ Hodge signature.
- `disc_nonpos_of_quadratic_nonneg` (DiscriminantTraceBound §1726): nonneg quadratic ⟹ b²≤4ac (the connector).
- Function-field |λ|=√q chain §1720–1726 complete; `elliptic_arith_hodge_master` §1083 (weak Mordell–Weil ⟹ sigPos=1);
  `FrobeniusIsometryEigenvalues` §1728 (MᴴM=q·I ⟹ |eig|=√q); `pairingFun_self_pos` (Arakelov height pos-def).
- **917 modules; the RH-spine gated GREEN twice (cy345, cy358: 7/7 then 10/10 clean).**

**On Mac `~/projects/riemann/RiemannLean/` (the criterion + sign-theorem lane):**
- Sign theorem: `sign_theorem_onset` etc. — Bessel-positivity chain + μ-sign algebra, 22 theorems, 0 sorries.
  (Consequence-side: the open step hD assumes the zero-set onset factorization; NOT a proof path.)
- Weil criterion: `weilForm_nonneg_iff` (states RH positivity criterion), `weil_posSemidef_of_dominates`
  (Q⪰0 ⟸ [M⪰0 machine-checked] + [(Q−M)⪰0 = hQM, one labeled transcendental interval step]). 31 theorems.

**Arb (organ-2, rigorous interval → exact-rational LDLᵀ):** finite Weil form certified PD on 6-function basis (exact pivots).

## B. The honest status (Charisma tiers)

- **HANDLE:** the axiom-clean Lean criteria above; the Arb finite-positivity certificates; the sign-theorem skeleton.
- The criterion `heig ⟺ RH` (Hilbert–Pólya) and `Weil positivity ⟺ RH` — machine-verified **reformulations**, not reductions to an easier problem. **This is a genuine, novel, formally-verified result; it is NOT a proof of RH.**
- **The single open input = Weil positivity ⟺ RH = Connes–Consani's own open conjecture** (scaling-site-square Hodge index / CvS λ_∞^even>0 convergence / prolate exact spectral realization — all the same wall). Heuristic + 113–334 OOM numerics, unproven.

## C. Attack-surface map — every route tried, and why it walls (9-seat Opus+Sonnet Roundtable + autonomous cy329-367)

| Route | Verdict |
|---|---|
| Box-moment / Favard→heig | un-circular from ξ, N(T)-law, determinate — but RH lives in a shrinking S(T) fluctuation (RH-blind smooth part) |
| Weil positivity (Groskin margin) | finite-band certified, margin plateaus in N / c^(−0.44N); cutoff-free needs T~c^(0.44N) (computational wall) |
| Alpöge–Furman proportion | ≥2/3 zeros on-line PROVEN unconditionally; past 2/3 = the mollifier θ-road (Bettin–Gonek θ=∞⟹RH, known-hard) |
| de Branges κ=1 / dBN / CCM convergence | ⊂ heig; convergence ⟺ RH |
| F₁ geometric (scaling-site square) | surface exists in C-C framework; its Hodge positivity = C-C's open conjecture; classical route Deligne-parity-blocked |
| Kato perturbation of CvS | prime perturbation Σ_p log p/√p ~2√c DIVERGES — cancellation not smallness → residual wall |
| Condensed math (Clausen–Scholze) | parity is a DUALITY fact, not categorical — relabels, doesn't evade. Non-starter |
| RMT / Katz–Sarnak | family-vs-individual gap; one-directional (proven over function fields where RH already known). Consequence-side |
| Model theory / o-minimality | category error — ζ's zero set violates o-minimality's finiteness. Dead lens |
| Quantum chaos (Berry–Keating, BBM) | every candidate H bottoms at self-adjointness = Weil positivity, relabeled |
| Sign theorem (Selberg S(γ)) | consequence-side; needs the joint law of S(γ) and γ (Selberg marginal-proven, joint-open) |
| Fully-Lean finite Weil certificate | blocked on missing Mathlib verified-interval-arithmetic infrastructure |

**Conclusion:** all routes bottom at ONE uniform wall — Weil positivity ⟺ RH — confirmed deep (not a wrong-framework artifact) by 9 expert lenses. No new proof handle in the reachable space.

## D. The remaining step (genuine external mathematics)

Prove the criterion's hypothesis: Weil positivity for ζ. Equivalent named forms: the scaling-site-square Hodge
index theorem (C-C `weilpos`), the CvS/prolate exact spectral realization (Connes–Moscovici, deficiency-(4,4)
self-adjoint extension), or the classical A(g)+B(g)≥P(g). This is C-C's active research frontier. Not closable by
autonomous search; needs new mathematics or C-C's proof.
