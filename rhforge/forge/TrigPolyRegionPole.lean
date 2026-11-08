import TrigPolyZetaRegion
import HpoleRiemannZetaBox

/-!
# S2932 — P-generic ζ zero-free region with the POLE input discharged (lane 1)

Feeds §1572 (`hpole_riemannZeta_box`: Re(−ζ'/ζ)(s) ≤ 1/(Re s − 1) + A₀ on 1 < Re s ≤ 2, |Im s| ≤ 1) into
S2926 (`zeta_region_P`).  Result: there is ONE absolute constant A₀ ≥ 0 such that for every nonnegative
cosine polynomial P (c₀,c₁ > 0, c_k ≥ 0), every ε > 0 with ε ≤ 𝒦 (so σ = 1+ε/𝒦 ≤ 2), and every zero β+iγ
with β < σ, the zero and growth inputs alone give
  β ≤ 1 − ε(c₁−c₀−ε)/((c₀+ε)𝒦),   𝒦 = c₀A₀ + c₁B + Σ c_k C_k.
The pole hypothesis is gone; remaining open inputs for 5-8-4-1: zero input at iγ and the growth constants
(S2928 + S2930 log absorption at 2γ, 3γ).
-/

open Complex

namespace TrigPolyRegion

theorem Tz_zero (σ γ : ℝ) : Tz σ γ 0 = (-logDeriv riemannZeta (σ : ℂ)).re := by
  unfold Tz; simp

theorem zeta_region_P_pole :
    ∃ A₀ : ℝ, 0 ≤ A₀ ∧ ∀ {γ β ε B : ℝ} {K : ℕ} (c : Fin (K + 2) → ℝ) (Cg : Fin K → ℝ),
      (∀ θ : ℝ, 0 ≤ ∑ k, c k * Real.cos (k * θ)) →
      0 < c 0 → 0 < c 1 → (∀ k, 0 ≤ c k) → 0 < ε → 0 ≤ B → (∀ k, 0 ≤ Cg k) →
      0 < c 0 * A₀ + c 1 * B + ∑ k : Fin K, c k.succ.succ * Cg k →
      ε ≤ c 0 * A₀ + c 1 * B + ∑ k : Fin K, c k.succ.succ * Cg k →
      β < 1 + ε / (c 0 * A₀ + c 1 * B + ∑ k : Fin K, c k.succ.succ * Cg k) →
      Tz (1 + ε / (c 0 * A₀ + c 1 * B + ∑ k : Fin K, c k.succ.succ * Cg k)) γ 1
        ≤ B - 1 / ((1 + ε / (c 0 * A₀ + c 1 * B + ∑ k : Fin K, c k.succ.succ * Cg k)) - β) →
      (∀ k : Fin K,
        Tz (1 + ε / (c 0 * A₀ + c 1 * B + ∑ k : Fin K, c k.succ.succ * Cg k)) γ (k + 2) ≤ Cg k) →
      β ≤ 1 - ε * (c 1 - c 0 - ε)
        / ((c 0 + ε) * (c 0 * A₀ + c 1 * B + ∑ k : Fin K, c k.succ.succ * Cg k)) := by
  obtain ⟨A₀, hA0, hA⟩ := hpole_riemannZeta_box
  refine ⟨A₀, hA0, ?_⟩
  intro γ β ε B K c Cg hP hc0 hc1 hck hε hB hCg hKpos hεK hβ hzero hgrowth
  apply zeta_region_P c Cg hP hc0 hc1 hck hε hA0 hB hCg hKpos hβ _ hzero hgrowth
  set 𝒦 := c 0 * A₀ + c 1 * B + ∑ k : Fin K, c k.succ.succ * Cg k
  have hq : 0 < ε / 𝒦 := div_pos hε hKpos
  have hq1 : ε / 𝒦 ≤ 1 := (div_le_one hKpos).mpr hεK
  rw [Tz_zero]
  have h := hA ((1 + ε / 𝒦 : ℝ) : ℂ) (by simp; linarith) (by simp; linarith) (by simp)
  simpa using h

end TrigPolyRegion

#print axioms TrigPolyRegion.zeta_region_P_pole
