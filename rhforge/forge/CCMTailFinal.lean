import CCMSymHalf
import CCMSchurTest

/-!
# S3001 — sharpest 𝟙^⊥ far-tail bound: confluent Löwner + Taylor ½ + Schur (lane 2)

`confluent_perp_schur_half`: g C² on [−V,V], |g''| ≤ D₂, ν injective in [−V,V], Σx = 0 ⇒
  |xᵀ·lowner(−g'∘ν/π, g∘ν, ν)·x| ≤ (D₂/(2π))·S·Σx²,   S ≥ max_i Σ_k |ν_i − ν_k|.
Composes S2998 (`lowner_entry_sym_half`, δ_ik = D₂|ν_i−ν_k|/(2π)), S2994 (`schur_quad`) and S2990 (`quad_rowsum_zero`).
Expected: ~10× from the true top eigenvalue on 𝟙^⊥ (S2994: 20× with the D₂/π constant).
-/

open Set Finset Matrix

namespace CCMTailFinal

theorem confluent_perp_schur_half {ι : Type*} [Fintype ι] [DecidableEq ι] (g g' g'' : ℝ → ℝ)
    {V D₂ S : ℝ}
    (hg : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g (g' x) (Icc (-V) V) x)
    (hg' : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g' (g'' x) (Icc (-V) V) x)
    (hD : ∀ x ∈ Icc (-V) V, |g'' x| ≤ D₂) (ν : ι → ℝ) (hν : ∀ i, ν i ∈ Icc (-V) V)
    (hinj : Function.Injective ν) (hS : ∀ i, ∑ k, |ν i - ν k| ≤ S) (x : ι → ℝ)
    (hx : ∑ i, x i = 0) :
    |x ⬝ᵥ (CCMLowner.lowner (fun i => -(g' (ν i)) / Real.pi) (fun j => g (ν j)) ν *ᵥ x)|
      ≤ D₂ / (2 * Real.pi) * S * ∑ i, x i ^ 2 := by
  have hpi := Real.pi_pos
  set u : ι → ℝ := fun i => -(g' (ν i)) / (2 * Real.pi) with hu
  set P : Matrix ι ι ℝ := vecMulVec u (fun _ => (1 : ℝ)) + vecMulVec (fun _ => (1 : ℝ)) u with hP
  set δ : Matrix ι ι ℝ := fun i k => D₂ / (2 * Real.pi) * |ν i - ν k| with hδ
  have hPe : ∀ i k, P i k = u i + u k := by intro i k; simp [hP, vecMulVec_apply]
  have hent : ∀ i k, |CCMLowner.lowner (fun i => -(g' (ν i)) / Real.pi) (fun j => g (ν j)) ν i k
      - P i k| ≤ δ i k := by
    intro i k
    rw [hPe]
    by_cases hik : i = k
    · subst hik
      simp only [CCMLowner.lowner, if_true, hu, hδ, sub_self, abs_zero, mul_zero]
      have : -(g' (ν i)) / Real.pi - (-(g' (ν i)) / (2 * Real.pi) + -(g' (ν i)) / (2 * Real.pi)) = 0 := by
        field_simp; ring
      rw [this, abs_zero]
    · simp only [CCMLowner.lowner, if_neg hik]
      have hne : ν i ≠ ν k := fun e => hik (hinj e)
      have hs := CCMSymHalf.lowner_entry_sym_half g g' g'' hg hg' hD (hν i) (hν k) hne
      have e : (g (ν k) - g (ν i)) / (Real.pi * (ν i - ν k)) - (u i + u k)
          = (g (ν k) - g (ν i)) / (Real.pi * (ν i - ν k)) + (g' (ν i) + g' (ν k)) / (2 * Real.pi) := by
        simp only [hu]; field_simp; ring
      rw [e]
      calc _ ≤ D₂ * |ν k - ν i| / (2 * Real.pi) := hs
        _ = δ i k := by simp only [hδ]; rw [abs_sub_comm]; ring
  have hsym : ∀ i k, δ i k = δ k i := by intro i k; simp only [hδ]; rw [abs_sub_comm]
  have hrow : ∀ i, ∑ k, δ i k ≤ D₂ / (2 * Real.pi) * S := by
    intro i
    have hD0 : 0 ≤ D₂ := le_trans (abs_nonneg _) (hD (ν i) (hν i))
    simp only [hδ]; rw [← Finset.mul_sum]
    exact mul_le_mul_of_nonneg_left (hS i) (by positivity)
  have hq := CCMSchur.schur_quad _ P δ _ hent hsym hrow x
  rw [hP, CCMLownerSym.quad_rowsum_zero u x hx, sub_zero] at hq
  exact hq

end CCMTailFinal

#print axioms CCMTailFinal.confluent_perp_schur_half
