import Mathlib.NumberTheory.Chebyshev
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic

set_option maxHeartbeats 1000000

/-!
# §1580 — ψ is constant on strips: `∫_n^{n+1} ψ(t)t^{-s-1} dt = ψ(n)·∫_n^{n+1} t^{-s-1} dt`

Third brick of the Mellin pillar.  `ψ(t) = ψ(⌊t⌋) = ψ(n)` for `t ∈ [n, n+1)` (`psi_eq_psi_coe_floor` +
`Nat.floor_eq_iff`); the endpoint `t = n+1` is a null set, so the interval integral sees the constant `ψ(n)`
(`integral_congr_ae` + `Real.volume_singleton`), which then factors out (`integral_const_mul`).  Combined with
§1579 this turns §1578's difference sum into a genuine integral over `[1, N]`.  EXIT0 / std-axioms.
-/

open Complex Chebyshev MeasureTheory

theorem psi_strip_integral {s : ℂ} {n : ℕ} (hn : 0 < n) :
    ∫ t : ℝ in (n : ℝ)..((n : ℝ) + 1), ((ψ t : ℝ) : ℂ) * (t : ℂ) ^ (-s - 1)
      = ((ψ (n : ℝ) : ℝ) : ℂ) * ∫ t : ℝ in (n : ℝ)..((n : ℝ) + 1), (t : ℂ) ^ (-s - 1) := by
  have hnR : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn
  have hsing : ∀ᵐ t : ℝ, t ≠ ((n : ℝ) + 1) := by
    refine MeasureTheory.ae_iff.mpr ?_
    have hset : {t : ℝ | ¬ t ≠ ((n : ℝ) + 1)} = {((n : ℝ) + 1)} := by
      ext t; simp
    rw [hset]
    exact Real.volume_singleton
  have hae : ∀ᵐ t : ℝ, t ∈ Set.uIoc (n : ℝ) ((n : ℝ) + 1) →
      ((ψ t : ℝ) : ℂ) * (t : ℂ) ^ (-s - 1) = ((ψ (n : ℝ) : ℝ) : ℂ) * (t : ℂ) ^ (-s - 1) := by
    filter_upwards [hsing] with t ht hmem
    rw [Set.uIoc_of_le (by linarith : (n : ℝ) ≤ (n : ℝ) + 1), Set.mem_Ioc] at hmem
    have hfl : ⌊t⌋₊ = n := by
      rw [Nat.floor_eq_iff (by linarith [hmem.1] : (0 : ℝ) ≤ t)]
      refine ⟨le_of_lt hmem.1, ?_⟩
      rcases lt_or_eq_of_le hmem.2 with h | h
      · exact_mod_cast h
      · exact absurd h ht
    have hψ : ψ t = ψ (n : ℝ) := by
      rw [psi_eq_psi_coe_floor t, hfl]
    rw [hψ]
  rw [intervalIntegral.integral_congr_ae hae, intervalIntegral.integral_const_mul]

#print axioms psi_strip_integral
