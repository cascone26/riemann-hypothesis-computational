import Psi1CMellinForm
import Mathlib.Analysis.MellinTransform

set_option maxHeartbeats 1000000

/-!
# §1599 — `MellinConvergent psi1C (-σ₀-1)` for `σ₀ > 1`

`MellinConvergent f σ' = IntegrableOn (t ↦ t^{σ'-1}•f t) (Ioi 0)`.  Split `Ioi 0 = Ioc 0 1 ∪ Ioi 1`: on `Ioc 0 1`
the integrand vanishes (`psi1C = 0`, §1597), and on `Ioi 1` it matches `psi1C(t)·t^{-σ₀-2}` (§1596(b)).  This is
the first hypothesis of Mathlib's Mellin inversion for the smoothed Perron formula.  EXIT0 / std-axioms.
-/

open Complex Chebyshev MeasureTheory Set

theorem mellinConvergent_psi1C {σ₀ : ℝ} (hσ₀ : 1 < σ₀) :
    MellinConvergent psi1C ((-σ₀ - 1 : ℝ) : ℂ) := by
  rw [MellinConvergent]
  have hsplit : Set.Ioi (0 : ℝ) = Set.Ioc 0 1 ∪ Set.Ioi 1 :=
    (Set.Ioc_union_Ioi_eq_Ioi (by norm_num : (0:ℝ) ≤ 1)).symm
  rw [hsplit]
  apply MeasureTheory.IntegrableOn.union
  · -- Ioc 0 1: integrand ≡ 0
    apply (MeasureTheory.integrableOn_zero).congr_fun _ measurableSet_Ioc
    intro t ht
    rw [Set.mem_Ioc] at ht
    simp only []
    rw [psi1C_eq_zero_of_le_one ht.2, smul_zero]
  · -- Ioi 1: matches §1596(b)
    have h96 := psi1C_cpow_integrableOn_Ioi (s := (σ₀ : ℂ)) (by rwa [Complex.ofReal_re])
    apply h96.congr_fun _ measurableSet_Ioi
    intro t ht
    rw [Set.mem_Ioi] at ht
    have ht0 : t ≠ 0 := by intro h; rw [h] at ht; linarith
    simp only [smul_eq_mul]
    rw [mul_comm]
    congr 2
    push_cast
    ring

#print axioms mellinConvergent_psi1C
