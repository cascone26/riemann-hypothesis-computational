import LFunctionZeroFreeSharpQuad

set_option maxHeartbeats 1000000

/-!
# §1567 — clean citable form of the QUADRATIC (`χ²=1`) Dirichlet-L zero-free region

Quadratic analogue of §1552.  Repackages §1566 into the textbook shape: from `β ≤ 1-1/(14(3A₀+320·log|γ|))`,
`c₀ = 1/8960` gives `β ≤ 1 - c₀/log|γ|` once `log|γ| ≥ 3A₀/320`; all (modulus-dependent) thresholds fold into a
single `T₀`.  For every nontrivial `χ` with `χ²=1`:

  **∃ T₀ c₀ > 0, ∀ γ β, T₀ ≤ |γ| → LFunction χ(β+iγ) = 0 → β ≤ 1 - c₀/log|γ|.**

EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem LFunction_zero_free_region_clean_quad {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) (hχ2eq : χ ^ 2 = 1) :
    ∃ T₀ c₀ : ℝ, 0 < c₀ ∧ ∀ γ β : ℝ, T₀ ≤ |γ| →
      LFunction χ ((β : ℂ) + I * (γ : ℝ)) = 0 → β ≤ 1 - c₀ / Real.log |γ| := by
  obtain ⟨A₀, hA₀nn, hmain⟩ := LFunction_zero_free_region_sharp_quad hχ hχ2eq
  refine ⟨max (max (Real.exp 16)
      (Real.exp (2 * Real.log ((N : ℝ) * (2 + 6) * (3 + 6 * A₀ + 640)))))
      (max (Real.exp (2 * Real.log (16 * (3 + 6 * A₀ + 640)) + Real.log N / 16))
      (Real.exp (3 * A₀ / 320))), 1 / 8960, by norm_num, ?_⟩
  intro γ β hγ hL
  have h1 : Real.exp 16 ≤ |γ| := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hγ
  have hexp2 : Real.exp (2 * Real.log ((N : ℝ) * (2 + 6) * (3 + 6 * A₀ + 640))) ≤ |γ| :=
    le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hγ
  have hexp2b : Real.exp (2 * Real.log (16 * (3 + 6 * A₀ + 640)) + Real.log N / 16) ≤ |γ| :=
    le_trans (le_trans (le_max_left _ _) (le_max_right _ _)) hγ
  have hexp3 : Real.exp (3 * A₀ / 320) ≤ |γ| := le_trans (le_trans (le_max_right _ _) (le_max_right _ _)) hγ
  have h2 : 2 * Real.log ((N : ℝ) * (2 + 6) * (3 + 6 * A₀ + 640)) ≤ Real.log |γ| := by
    calc 2 * Real.log ((N : ℝ) * (2 + 6) * (3 + 6 * A₀ + 640))
        = Real.log (Real.exp (2 * Real.log ((N : ℝ) * (2 + 6) * (3 + 6 * A₀ + 640)))) :=
          (Real.log_exp _).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos _) hexp2
  have h2b : 2 * Real.log (16 * (3 + 6 * A₀ + 640)) + Real.log N / 16 ≤ Real.log |γ| := by
    calc 2 * Real.log (16 * (3 + 6 * A₀ + 640)) + Real.log N / 16
        = Real.log (Real.exp (2 * Real.log (16 * (3 + 6 * A₀ + 640)) + Real.log N / 16)) :=
          (Real.log_exp _).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos _) hexp2b
  have hL16 : 16 ≤ Real.log |γ| := by
    calc (16 : ℝ) = Real.log (Real.exp 16) := (Real.log_exp 16).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos 16) h1
  have hLpos : 0 < Real.log |γ| := by linarith
  have hL3A : 3 * A₀ / 320 ≤ Real.log |γ| := by
    calc 3 * A₀ / 320 = Real.log (Real.exp (3 * A₀ / 320)) := (Real.log_exp _).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos _) hexp3
  have hmain' := hmain γ β h1 h2 h2b hL
  have hKpos : 0 < 14 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|) := by
    nlinarith [hA₀nn, hLpos]
  have hcmp : 14 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|) ≤ 8960 * Real.log |γ| := by
    nlinarith [hL3A]
  have hdiv : 1 / (8960 * Real.log |γ|)
      ≤ 1 / (14 * (3 * A₀ + 4 * (64 * Real.log |γ|) + 64 * Real.log |γ|)) :=
    one_div_le_one_div_of_le hKpos hcmp
  have hc0 : (1 / 8960) / Real.log |γ| = 1 / (8960 * Real.log |γ|) := by rw [div_div]
  rw [hc0]
  linarith [hmain', hdiv]

#print axioms LFunction_zero_free_region_clean_quad
