import SmulFactorEqSplitLFunction
import EqOfCodiscreteContinuous
import SplitFunAnalytic

/-!
# §1481 — `LFunction χ z = P_half(z)·(P_rest(z)·g z)` on the sphere (§1471's `hfac`)

L-analogue of ζ §1391, cleaner (entire — no pole/`hball`).  Rewriting §1476's codiscrete factorization into split
form (§1480), then upgrading to pointwise on the sphere via the generic §1385
`eq_of_codiscreteWithin_of_continuousAt` — using that both `LFunction χ` (`differentiable_LFunction`, continuous)
and the split function (§1390, analytic ⟹ continuous) are continuous at every `z ∈ sphere c R ⊆ ball c R'`
(`R < R'`).  This is exactly the `hfac` hypothesis the generic Borel–Carathéodory factor bound §1471 needs.
EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric Classical Filter

theorem LFunction_eq_split_on_sphere {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {c : ℂ} {R' R : ℝ} (hRR' : R < R') {g : ℂ → ℂ} (hg : AnalyticOnNhd ℂ g (ball c R'))
    (hcod : LFunction χ =ᶠ[codiscreteWithin (ball c R')]
      (∏ᶠ u, (· - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u)) • g)
    (p : ℂ → Prop) {z : ℂ} (hz : z ∈ sphere c R) :
    LFunction χ z
      = (∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter p,
          (z - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat)
        * ((∏ u ∈ (divisor_LFunction_ball_support_finite hχ c R').toFinset.filter (fun u => ¬ p u),
            (z - u) ^ (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) * g z) := by
  rw [smul_factor_eq_split_LFunction hχ p] at hcod
  have hzball : z ∈ ball c R' := by
    rw [mem_sphere] at hz; rw [mem_ball, hz]; exact hRR'
  refine eq_of_codiscreteWithin_of_continuousAt isOpen_ball hzball hcod ?_ ?_
  · exact (differentiable_LFunction hχ z).continuousAt
  · exact (split_fun_analyticOnNhd
      (fun u => (MeromorphicOn.divisor (LFunction χ) (ball c R') u).toNat) hg z hzball).continuousAt

#print axioms LFunction_eq_split_on_sphere
