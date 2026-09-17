import Mathlib.NumberTheory.LSeries.DirichletContinuation
import Mathlib.NumberTheory.LSeries.Dirichlet

/-!
# §1360 — uniform (absolute-constant) bound `‖logDeriv g σ‖ ≤ A` on `[1,2]` (§1359's `A`)

§1359's pole bound carries `A = ‖logDeriv g σ‖` (point-dependent).  Making it an *absolute constant* — needed so
the region `K = 3A+4B+C` has `A` independent of `γ` — is a compactness argument: `-logDeriv g = -deriv g / g`
(`g = LFunctionTrivChar₁ 1 = (s-1)ζ`) is continuous on `{s | s = 1 ∨ LFunctionTrivChar 1 s ≠ 0}`
(`continuousOn_neg_logDeriv_LFunctionTrivChar₁`), and `[1,2] ⊆` that set (`s = 1` at the left endpoint; for
`σ ∈ (1,2]`, `LFunctionTrivChar 1 σ = ζ(σ) ≠ 0`).  A compact interval + continuity ⟹ a uniform bound `A`.
So for every DLVP abscissa `σ = 1+1/(2K)` with `K ≥ 1/2` (`σ ≤ 2`), §1359's `A` is this absolute constant.
EXIT0 / std-axioms.
-/

open Complex ArithmeticFunction LSeries DirichletCharacter

theorem exists_logDeriv_LFunctionTrivChar₁_bound :
    ∃ A : ℝ, 0 ≤ A ∧ ∀ σ : ℝ, σ ∈ Set.Icc (1 : ℝ) 2 →
      ‖logDeriv (LFunctionTrivChar₁ 1) (σ : ℂ)‖ ≤ A := by
  have hmaps : Set.MapsTo (fun σ : ℝ => (σ : ℂ)) (Set.Icc (1 : ℝ) 2)
      {s | s = 1 ∨ LFunctionTrivChar 1 s ≠ 0} := by
    intro σ hσ
    rcases eq_or_lt_of_le hσ.1 with h1 | h1
    · left; rw [← h1]; simp
    · right
      have hσ1 : (σ : ℂ) ≠ 1 := by rw [Ne, Complex.ofReal_eq_one]; linarith
      rw [LFunctionTrivChar_eq_mul_riemannZeta hσ1, Nat.primeFactors_one, Finset.prod_empty, one_mul]
      exact riemannZeta_ne_zero_of_one_lt_re (by simpa using h1)
  have hcont : ContinuousOn
      (fun σ : ℝ => (fun s => -deriv (LFunctionTrivChar₁ 1) s / LFunctionTrivChar₁ 1 s) (σ : ℂ))
      (Set.Icc (1 : ℝ) 2) :=
    (continuousOn_neg_logDeriv_LFunctionTrivChar₁ (n := 1)).comp
      Complex.continuous_ofReal.continuousOn hmaps
  obtain ⟨C, hC⟩ := (isCompact_Icc).exists_bound_of_continuousOn hcont
  refine ⟨max C 0, le_max_right _ _, fun σ hσ => ?_⟩
  have hval := hC σ hσ
  have heq : logDeriv (LFunctionTrivChar₁ 1) (σ : ℂ)
      = -(-deriv (LFunctionTrivChar₁ 1) (σ : ℂ) / LFunctionTrivChar₁ 1 (σ : ℂ)) := by
    rw [logDeriv_apply, neg_div, neg_neg]
  rw [heq, norm_neg]
  exact le_trans hval (le_max_left _ _)

#print axioms exists_logDeriv_LFunctionTrivChar₁_bound
