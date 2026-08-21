# Mathlib Name/API Cheatsheet — verified working in our toolchain (Lean 4.33.0)
_Every entry below cost a compile round-trip to discover this session. Consult BEFORE guessing a name.
Append new ones as they're found. This directly kills the #1 time sink (name-churn)._

## Renamed / moved (the guess is WRONG — use the RHS)
| You'll reach for | It's actually |
|---|---|
| `le_or_lt` | GONE — use `le_total a b` (gives `a≤b ∨ b≤a`) or `by_cases h : a ≤ b` + `not_le.mp` |
| `tsum_le_tsum h hf hg` | `Summable.tsum_le_tsum` → `hf.tsum_le_tsum h hg` (summable-of-f is the receiver) |
| `Function.update_same` | `Function.update_self` |
| `Ioo_mem_nhdsWithin_Ioi` | `Ioo_mem_nhdsGT` |
| `div_le_div_of_nonneg_right` | doesn't fit — use `gcongr`, or `div_le_div_of_nonneg_left (ha:0≤a)(hc:0<c)(h:c≤b) : a/b ≤ a/c` (note arg order: numerator-nonneg, then SMALL denom pos, then denom-order) |
| `Nat.Ico_succ_right` / `Finset.range_succ_eq_Icc_zero` | `Nat.range_succ_eq_Icc_zero` |
| `countable_empty` | `Set.countable_empty` |

## Hidden gotchas (name is right, but there's a trap)
- `LSeries.norm_term_eq` : RHS has a **hidden `if n = 0 then 0 else ‖f n‖/↑n^s.re`** — after any positivity
  case-split, `rw [if_neg hn.ne']` (hn : 0<n) to clear it before the div lemma.
- `ContinuousAt.comp` / `ContinuousAt.div` on a `fun t => g (f t)` : the composed goal is a **beta-redex**;
  `.comp`/`.div` won't unify directly. Bind an explicit-typed `have hcomp : ContinuousAt (fun t => …) t := …`.
- Complex `const_cpow` continuity: `fun_prop` FAILS. Use `Continuous.const_cpow (hf) (Or.inl hbase_ne)`
  (base ≠ 0) or `DifferentiableAt.const_cpow`.
- `differentiableAt_id` under a `ContinuousAt`/no-ℂ-typed goal: field 𝕜 is a stuck metavar → annotate
  `(differentiableAt_id (𝕜 := ℂ))`.
- `integral_boundary_rect_eq_zero_of_differentiableOn` : `apply` can't HO-unify `f` → `refine … (fun s => …) z w ?_`.
- `Set.EqOn` goal after `intro` is a beta-redex → `simp only [update_of_ne …]` (or `simp only []`), NOT `rw`.
- `field_simp` recurses into an atom's internal division (e.g. F's `/(s(s+1))`) and demands its ≠0 side goals —
  `generalize <the raw body> = A` first so field_simp only clears the outer denominator.
- `rw [intervalIntegral.integral_const]` then `simp` may CLOSE the goal → a following `ring`/`ring_nf` errors
  "no goals"; use `rw […, smul_eq_mul]; ring` instead of `simp; ring`.
- `simpa` reducing `((σ:ℂ)+t*I - 1).re = 0` : it leaves `σ - 1 = 0`, needs a following `linarith` (not auto).

## Confirmed working idioms (reuse verbatim)
- `‖L(Λ)(s)‖ ≤ ‖L(Λ)(Re s)‖` : §1591 `norm_LSeries_vonMangoldt_le_line`. Monotone in σ: §1628.
- `logDeriv riemannZeta s = -LSeries (↗Λ) s` (Re s>1) : `logDeriv_riemannZeta_eq_neg_vonMangoldt`.
- `riemannZeta_ne_zero_of_one_lt_re` : needs `import Mathlib.NumberTheory.LSeries.Dirichlet`.
- ζ analytic off 1 : `analyticOn_riemannZeta : AnalyticOnNhd ℂ riemannZeta {1}ᶜ`; `.deriv` for ζ'.
- `‖(x:ℂ)^y‖ = x^y.re` (x>0) : `Complex.norm_cpow_eq_rpow_re_of_pos`.
- residue/Cauchy : circle only — `circleIntegral.integral_sub_inv_of_mem_ball` (∮(z-w)⁻¹=2πi),
  `circleIntegral_div_sub_of_differentiable_on_off_countable` (Cauchy formula, countable exceptions).
  Rectangle: `integral_boundary_rect_eq_zero_of_differentiableOn` (VANISHING only — **no rectangle-winding lemma exists**).
- removable singularity : `analyticAt_of_differentiable_on_punctured_nhds_of_continuousAt` +
  `continuousAt_update_same` (⟺ Tendsto punctured) + `differentiableAt_dslope_of_ne`.
- `d/dt Complex.log(f t)` : `HasDerivAt.clog_real` (needs `f t ∈ slitPlane`, i.e. `mem_slitPlane_iff.mpr (Or.inl re_pos)`).
- `d/dt (↑t:ℂ) = 1` : `(hasDerivAt_id t).ofReal_comp`.
- FTC : `intervalIntegral.integral_eq_sub_of_hasDerivAt (hderiv ∀t∈uIcc) (hint : IntervalIntegrable)`.
- interval-integrable-from-continuous : `ContinuousOn.intervalIntegrable (hu : ContinuousOn u (uIcc a b))`.

## FASTER THAN GUESSING: use the search tactics (installed, we've been ignoring them)
- `exact?` / `apply?` — find the lemma that closes the goal by type. Run it BEFORE hand-writing a term.
- `rw?` — find a rewrite. `LeanSearchClient` is present in our Mathlib (`#leansearch "…"` / `#loogle`).
