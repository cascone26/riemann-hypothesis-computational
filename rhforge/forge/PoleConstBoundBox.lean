import Mathlib.NumberTheory.LSeries.DirichletContinuation
import Mathlib.NumberTheory.LSeries.Dirichlet
import Mathlib.NumberTheory.LSeries.Nonvanishing

set_option maxHeartbeats 1000000

/-!
# §1571 — uniform pole-constant on the COMPLEX rectangle `[1,2] × [-1,1]`

Complex-rectangle version of §1360.  `-deriv g / g` (`g = LFunctionTrivChar₁ 1 = (s-1)ζ`) is continuous on
`{s = 1 ∨ LFunctionTrivChar 1 s ≠ 0}`; the compact rectangle `S = {1 ≤ Re ≤ 2, |Im| ≤ 1}` sits inside that set
because `ζ(s) ≠ 0` for `Re s ≥ 1` (`riemannZeta_ne_zero_of_one_le_re` — the 1-LINE nonvanishing, junk-value-safe
at `s=1`).  Compactness (closed + bounded in the proper space ℂ) gives a uniform `A₀`.  This is the pole constant
for the SMALL-|t| region: unlike §1360 (real `[1,2]` only), it covers the whole strip near `s=1` that the
small-`γ` DLVP test points occupy.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem exists_logDeriv_LFunctionTrivChar₁_bound_box :
    ∃ A : ℝ, 0 ≤ A ∧ ∀ s : ℂ, 1 ≤ s.re → s.re ≤ 2 → |s.im| ≤ 1 →
      ‖logDeriv (LFunctionTrivChar₁ 1) s‖ ≤ A := by
  set S : Set ℂ := Complex.re ⁻¹' Set.Icc (1 : ℝ) 2 ∩ Complex.im ⁻¹' Set.Icc (-1 : ℝ) 1 with hSdef
  have hclosed : IsClosed S :=
    (isClosed_Icc.preimage Complex.continuous_re).inter (isClosed_Icc.preimage Complex.continuous_im)
  have hsubball : S ⊆ Metric.closedBall (0 : ℂ) 3 := by
    intro s hs
    obtain ⟨hre, him⟩ := hs
    rw [Set.mem_preimage, Set.mem_Icc] at hre him
    rw [Metric.mem_closedBall, dist_zero_right]
    calc ‖s‖ = ‖(s.re : ℂ) + (s.im : ℂ) * I‖ := by rw [Complex.re_add_im]
      _ ≤ ‖((s.re : ℝ) : ℂ)‖ + ‖((s.im : ℝ) : ℂ) * I‖ := norm_add_le _ _
      _ = |s.re| + |s.im| := by
          rw [Complex.norm_real, Real.norm_eq_abs, norm_mul, Complex.norm_I, mul_one,
            Complex.norm_real, Real.norm_eq_abs]
      _ ≤ 3 := by
          have h1 : |s.re| ≤ 2 := abs_le.mpr ⟨by linarith [hre.1], hre.2⟩
          have h2 : |s.im| ≤ 1 := abs_le.mpr ⟨him.1, him.2⟩
          linarith
  have hbounded : Bornology.IsBounded S := Metric.isBounded_closedBall.subset hsubball
  have hcompact : IsCompact S := Metric.isCompact_of_isClosed_isBounded hclosed hbounded
  have hsubset : S ⊆ {s : ℂ | s = 1 ∨ LFunctionTrivChar 1 s ≠ 0} := by
    intro s hs
    obtain ⟨hre, _⟩ := hs
    rw [Set.mem_preimage, Set.mem_Icc] at hre
    by_cases h1 : s = 1
    · left; exact h1
    · right
      rw [LFunctionTrivChar_eq_mul_riemannZeta h1, Nat.primeFactors_one, Finset.prod_empty, one_mul]
      exact riemannZeta_ne_zero_of_one_le_re hre.1
  have hcont : ContinuousOn
      (fun s => -deriv (LFunctionTrivChar₁ 1) s / LFunctionTrivChar₁ 1 s) S :=
    (continuousOn_neg_logDeriv_LFunctionTrivChar₁ (n := 1)).mono hsubset
  obtain ⟨C, hC⟩ := hcompact.exists_bound_of_continuousOn hcont
  refine ⟨max C 0, le_max_right _ _, fun s h1 h2 h3 => ?_⟩
  have hmem : s ∈ S := by
    refine ⟨?_, ?_⟩
    · rw [Set.mem_preimage, Set.mem_Icc]; exact ⟨h1, h2⟩
    · rw [Set.mem_preimage, Set.mem_Icc]; exact abs_le.mp h3
  have hval := hC s hmem
  have heq : logDeriv (LFunctionTrivChar₁ 1) s
      = -(-deriv (LFunctionTrivChar₁ 1) s / LFunctionTrivChar₁ 1 s) := by
    rw [logDeriv_apply, neg_div, neg_neg]
  rw [heq, norm_neg]
  exact le_trans hval (le_max_left _ _)

#print axioms exists_logDeriv_LFunctionTrivChar₁_bound_box
