import CCMLownerSym

/-!
# S2992 — confluent Löwner matrices on 𝟙^⊥: |xᵀMx| ≤ (2V·sup|g''|/π)(Σ|x|)² (lane 2)

For the Cauchy far tail T = Σ_j w_j·(pair kernel), the diagonal is EXACTLY −g'_tail(ν_i)/π
(g_tail(ν) = πΣw(1/(ν−x) + 1/(ν+x)), g'_tail = −πΣw(1/(ν−x)² + 1/(ν+x)²)), i.e. T is the confluent Löwner (Pick-type)
matrix of g_tail.  For such matrices:
* `confluent_lowner_perp` — g C² on [−V,V] with |g''| ≤ D₂, ν injective in [−V,V], Σx = 0 ⇒
  |xᵀ·lowner(−g'∘ν/π, g∘ν, ν)·x| ≤ (2V·D₂/π)·(Σ|x|)².
This is the closed-form 𝟙^⊥ bound for the whole far tail (diagonal included).
-/

open Set Finset Matrix

namespace CCMConfluent

theorem confluent_lowner_perp {ι : Type*} [Fintype ι] [DecidableEq ι] (g g' g'' : ℝ → ℝ)
    {V D₂ : ℝ}
    (hg : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g (g' x) (Icc (-V) V) x)
    (hg' : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g' (g'' x) (Icc (-V) V) x)
    (hD : ∀ x ∈ Icc (-V) V, |g'' x| ≤ D₂) (ν : ι → ℝ) (hν : ∀ i, ν i ∈ Icc (-V) V)
    (hinj : Function.Injective ν) (x : ι → ℝ) (hx : ∑ i, x i = 0) :
    |x ⬝ᵥ (CCMLowner.lowner (fun i => -(g' (ν i)) / Real.pi) (fun j => g (ν j)) ν *ᵥ x)|
      ≤ 2 * V * D₂ / Real.pi * (∑ i, |x i|) ^ 2 := by
  have hpi := Real.pi_pos
  refine CCMLownerSym.near_rowsum_quad_perp _ (fun i => -(g' (ν i)) / (2 * Real.pi)) _
    (fun i k => ?_) x hx
  have hVi := hν i
  have hV0 : 0 ≤ V := by linarith [hVi.1, hVi.2]
  have hD0 : 0 ≤ D₂ := le_trans (abs_nonneg _) (hD (ν i) hVi)
  by_cases hik : i = k
  · subst hik
    simp only [CCMLowner.lowner, if_true]
    have : -(g' (ν i)) / Real.pi - (-(g' (ν i)) / (2 * Real.pi) + -(g' (ν i)) / (2 * Real.pi)) = 0 := by
      field_simp; ring
    rw [this, abs_zero]; positivity
  · simp only [CCMLowner.lowner, if_neg hik]
    have hne : ν i ≠ ν k := fun e => hik (hinj e)
    have hs := CCMLownerSym.lowner_entry_sym g g' g'' hg hg' hD (hν i) (hν k) hne
    have e : (g (ν k) - g (ν i)) / (Real.pi * (ν i - ν k))
        - (-(g' (ν i)) / (2 * Real.pi) + -(g' (ν k)) / (2 * Real.pi))
        = (g (ν k) - g (ν i)) / (Real.pi * (ν i - ν k)) + (g' (ν i) + g' (ν k)) / (2 * Real.pi) := by
      field_simp; ring
    rw [e]
    have hdist : |ν k - ν i| ≤ 2 * V := by
      rw [abs_le]; constructor <;> linarith [(hν i).1, (hν i).2, (hν k).1, (hν k).2]
    calc _ ≤ D₂ * |ν k - ν i| / Real.pi := hs
      _ ≤ D₂ * (2 * V) / Real.pi := by gcongr
      _ = 2 * V * D₂ / Real.pi := by ring

end CCMConfluent

#print axioms CCMConfluent.confluent_lowner_perp
