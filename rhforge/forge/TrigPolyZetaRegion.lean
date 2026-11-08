/-
TrigPolyZetaRegion.lean (HP RH center, LEDGER S2926, 2026-11-08) — lane 1.
ζ-level composition of S2921 (general-P positivity) + S2922 (P-parametrized gap/region): for ANY nonnegative
cosine polynomial P = Σ_{k<K+2} c_k cos(kθ) with c₀, c₁ > 0, c_k ≥ 0, given the three analytic input bounds
at σ = 1 + ε/𝒦 (pole at σ, zero at σ+iγ, growth at σ+ikγ for k ≥ 2), a ζ zero β+iγ satisfies
  β ≤ 1 − ε(c₁ − c₀ − ε)/((c₀ + ε)𝒦),   𝒦 = c₀A + c₁B + Σ_{k≥2} c_k C_k.
This is the P-generic analogue of §1310 (which is P = 3-4-1); the analytic inputs are the existing tower's.
-/
import TrigPolyPositivityZeta
import TrigPolyZeroGap

open Complex

namespace TrigPolyRegion

/-- Re(−ζ'/ζ) at σ + ikγ. -/
noncomputable def Tz (σ γ : ℝ) (k : ℕ) : ℝ :=
  (-logDeriv riemannZeta ((σ : ℂ) + I * (((k : ℝ) * γ : ℝ) : ℂ))).re

theorem zeta_region_P {γ β ε A B : ℝ} {K : ℕ} (c : Fin (K + 2) → ℝ) (Cg : Fin K → ℝ)
    (hP : ∀ θ : ℝ, 0 ≤ ∑ k, c k * Real.cos (k * θ))
    (hc0 : 0 < c 0) (hc1 : 0 < c 1) (hck : ∀ k, 0 ≤ c k) (hε : 0 < ε)
    (hA : 0 ≤ A) (hB : 0 ≤ B) (hCg : ∀ k, 0 ≤ Cg k)
    (hKpos : 0 < c 0 * A + c 1 * B + ∑ k : Fin K, c k.succ.succ * Cg k)
    (hβ : β < 1 + ε / (c 0 * A + c 1 * B + ∑ k : Fin K, c k.succ.succ * Cg k))
    (hpole : Tz (1 + ε / (c 0 * A + c 1 * B + ∑ k : Fin K, c k.succ.succ * Cg k)) γ 0
      ≤ 1 / ((1 + ε / (c 0 * A + c 1 * B + ∑ k : Fin K, c k.succ.succ * Cg k)) - 1) + A)
    (hzero : Tz (1 + ε / (c 0 * A + c 1 * B + ∑ k : Fin K, c k.succ.succ * Cg k)) γ 1
      ≤ B - 1 / ((1 + ε / (c 0 * A + c 1 * B + ∑ k : Fin K, c k.succ.succ * Cg k)) - β))
    (hgrowth : ∀ k : Fin K,
      Tz (1 + ε / (c 0 * A + c 1 * B + ∑ k : Fin K, c k.succ.succ * Cg k)) γ (k + 2) ≤ Cg k) :
    β ≤ 1 - ε * (c 1 - c 0 - ε)
      / ((c 0 + ε) * (c 0 * A + c 1 * B + ∑ k : Fin K, c k.succ.succ * Cg k)) := by
  set 𝒦 := c 0 * A + c 1 * B + ∑ k : Fin K, c k.succ.succ * Cg k with h𝒦
  set σ := 1 + ε / 𝒦 with hσdef
  have hσ : 1 < σ := by have : 0 < ε / 𝒦 := div_pos hε hKpos; linarith
  have hpos := TrigPolyDLVP.trig_positivity_riemannZeta (γ := γ) hσ c hP
  have hsum : ∑ k : Fin (K + 2), c k *
      (-logDeriv riemannZeta ((σ : ℂ) + I * (((k : ℝ) * γ : ℝ) : ℂ))).re
      = c 0 * Tz σ γ 0 + c 1 * Tz σ γ 1 + ∑ k : Fin K, c k.succ.succ * Tz σ γ (k + 2) := by
    rw [Fin.sum_univ_succ, Fin.sum_univ_succ]
    simp only [Tz, Fin.val_zero, Fin.val_succ, Fin.succ_zero_eq_one, Fin.val_one, Nat.cast_add,
      Nat.cast_ofNat, Nat.cast_zero, Nat.cast_one]
    ring_nf
  rw [hsum] at hpos
  have hR : ∑ k : Fin K, c k.succ.succ * Tz σ γ (k + 2) ≤ ∑ k : Fin K, c k.succ.succ * Cg k :=
    Finset.sum_le_sum fun k _ => mul_le_mul_of_nonneg_left (hgrowth k) (hck _)
  have hD : 0 ≤ ∑ k : Fin K, c k.succ.succ * Cg k :=
    Finset.sum_nonneg fun k _ => mul_nonneg (hck _) (hCg k)
  exact TrigPolyGap.region_of_bounds hc0 hc1 hε hA hB hD hKpos hβ hpos hpole hzero hR

end TrigPolyRegion

#print axioms TrigPolyRegion.zeta_region_P
