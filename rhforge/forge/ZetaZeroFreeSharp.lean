import RegionR14Clean
import LogSharpClose
import BetaLtSigmaK
import PoleConstBound

/-!
# §1462 ★★★★★ — SHARP unconditional ζ zero-free region `β ≤ 1 - c/log|γ|`

Same tower as §1460 but with the OPTIMAL constant choice `B = C = 64·log|γ|` (so `K = 3A₀+320·log|γ| ≍ log|γ|`),
closing the two log bounds via §1461 (log-log absorption).  Result: an absolute `A₀ ≥ 0` with, for every ζ zero
`β+iγ` past an explicit threshold,

  `β ≤ 1 - 1/(14·(3A₀ + 4·(64 log|γ|) + 64 log|γ|))`   ( `= 1 - c/log|γ|`, the classical de la Vallée Poussin region ).

EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem riemannZeta_zero_free_region_sharp :
    ∃ A₀ : ℝ, 0 ≤ A₀ ∧ ∀ γ β : ℝ, Real.exp 16 ≤ |γ| →
      2 * Real.log (16 * (3 + 6 * A₀ + 640)) ≤ Real.log |γ| →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 →
      β ≤ 1 - 1 / (14 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|)) := by
  obtain ⟨A₀, hA₀nn, hA₀bd⟩ := exists_logDeriv_LFunctionTrivChar₁_bound
  refine ⟨A₀, hA₀nn, ?_⟩
  intro γ β hγe hγthr hζ
  have he1 : (1 : ℝ) < Real.exp 16 := by nlinarith [Real.add_one_le_exp (16 : ℝ)]
  have hγ : (1 : ℝ) < |γ| := by linarith
  have hL : 16 ≤ Real.log |γ| := by
    calc (16 : ℝ) = Real.log (Real.exp 16) := (Real.log_exp 16).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos 16) hγe
  have hLpos : 0 < Real.log |γ| := by linarith
  have hB : 0 ≤ 64 * Real.log |γ| := by linarith
  have hK0 : 0 < 3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ| := by nlinarith [hA₀nn, hLpos]
  have hKbig : 32 / 7 ≤ 3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ| := by
    nlinarith [hA₀nn, hL]
  have hσKmem : 1 + 1 / (2 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|))
      ∈ Set.Icc (1 : ℝ) 2 := by
    constructor
    · have : (0 : ℝ) < 1 / (2 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|)) := by positivity
      linarith
    · have : 1 / (2 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|)) ≤ 1 := by
        rw [div_le_one (by positivity)]; linarith
      linarith
  have hA_pole : ‖logDeriv (LFunctionTrivChar₁ 1)
      ((((1 + 1 / (2 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|)) : ℝ)) : ℂ))‖ ≤ A₀ :=
    hA₀bd _ hσKmem
  have hBclean : 32 * Real.log ((3 * |γ| + 10)
      * (3 + 2 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|))) ≤ 64 * Real.log |γ| := by
    have h := log_sharp_close (a := A₀) (x := |γ|) (d := 3) hA₀nn (by norm_num) (by norm_num) hγe hγthr
    have heq : (3 + 2 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|))
        = 3 + 6 * A₀ + 640 * Real.log |γ| := by ring
    rw [heq]; exact h
  have hCclean : 32 * Real.log ((3 * |2 * γ| + 10)
      * (3 + 2 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|))) ≤ 64 * Real.log |γ| := by
    have h := log_sharp_close (a := A₀) (x := |γ|) (d := 6) hA₀nn (by norm_num) (by norm_num) hγe hγthr
    have heq : (3 * |2 * γ| + 10)
          * (3 + 2 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|))
        = (6 * |γ| + 10) * (3 + 6 * A₀ + 640 * Real.log |γ|) := by
      rw [abs_mul]; simp only [abs_two]; ring
    rw [heq]; exact h
  exact region_riemannZeta_R14_clean hγ hKbig hA₀nn hB hB
    (beta_lt_sigmaK hK0 hζ) hζ hBclean hA_pole hCclean

#print axioms riemannZeta_zero_free_region_sharp
