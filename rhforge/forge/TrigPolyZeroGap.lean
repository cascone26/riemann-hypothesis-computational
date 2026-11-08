/-
TrigPolyZeroGap.lean (HP RH center, LEDGER S2922, 2026-11-08) — lane 1.
§1310 (`dlvp_gap_of_four_bounds`) parametrized by the trig polynomial: with the positivity
0 ≤ c₀T₀ + c₁T₁ + R (S2921 supplies it, R = Σ_{k≥2} c_k T_k), pole bound T₀ ≤ 1/(σ−1)+A, zero bound
T₁ ≤ B − 1/(σ−β), and R ≤ D, the gap is σ−β ≥ c₁(σ−1)/(c₀ + K(σ−1)), K = c₀A + c₁B + D.
Region: σ = 1 + ε/K gives β ≤ 1 − c(P,ε)/K with c(P,ε) = ε(c₁ − c₀ − ε)/(c₀ + ε) — explicit in P.
(3-4-1: c₀=3,c₁=4; 5-8-4-1: c₀=5,c₁=8 — larger c₁−c₀ widens the admissible ε range.)
-/
import Mathlib

namespace TrigPolyGap

theorem gap_of_bounds {c₀ c₁ σ β T₀ T₁ R A B D : ℝ}
    (hc₀ : 0 < c₀) (hc₁ : 0 < c₁) (hσ : 1 < σ) (hβ : β < σ)
    (hA : 0 ≤ A) (hB : 0 ≤ B) (hD : 0 ≤ D)
    (hpos : 0 ≤ c₀ * T₀ + c₁ * T₁ + R)
    (hpole : T₀ ≤ 1 / (σ - 1) + A)
    (hzero : T₁ ≤ B - 1 / (σ - β))
    (hR : R ≤ D) :
    c₁ * (σ - 1) / (c₀ + (c₀ * A + c₁ * B + D) * (σ - 1)) ≤ σ - β := by
  set K := c₀ * A + c₁ * B + D with hK
  have hs1 : 0 < σ - 1 := by linarith
  have hsb : 0 < σ - β := by linarith
  have hKn : 0 ≤ K := by positivity
  have p0 := mul_le_mul_of_nonneg_left hpole hc₀.le
  have p1 := mul_le_mul_of_nonneg_left hzero hc₁.le
  have h1 : c₁ * (1 / (σ - β)) ≤ c₀ * (1 / (σ - 1)) + K := by
    rw [hK]; nlinarith [p0, p1]
  have h2 := mul_le_mul_of_nonneg_right h1 (le_of_lt (mul_pos hsb hs1))
  have e1 : c₁ * (1 / (σ - β)) * ((σ - β) * (σ - 1)) = c₁ * (σ - 1) := by
    field_simp
  have e2 : (c₀ * (1 / (σ - 1)) + K) * ((σ - β) * (σ - 1))
      = c₀ * (σ - β) + K * (σ - 1) * (σ - β) := by
    field_simp
  rw [div_le_iff₀ (by positivity)]
  nlinarith [h2, e1, e2]

/-- **Explicit P-dependent zero-free region**: β ≤ 1 − ε(c₁ − c₀ − ε)/((c₀+ε)K) at σ = 1 + ε/K. -/
theorem region_of_gap {c₀ c₁ ε K β : ℝ} (hc₀ : 0 < c₀) (hε : 0 < ε) (hK : 0 < K)
    (hgap : c₁ * (ε / K) / (c₀ + K * (ε / K)) ≤ (1 + ε / K) - β) :
    β ≤ 1 - ε * (c₁ - c₀ - ε) / ((c₀ + ε) * K) := by
  have hKε : K * (ε / K) = ε := by field_simp
  rw [hKε] at hgap
  have : c₁ * (ε / K) / (c₀ + ε) = c₁ * ε / ((c₀ + ε) * K) := by
    field_simp
  rw [this] at hgap
  have hsplit : 1 - ε * (c₁ - c₀ - ε) / ((c₀ + ε) * K)
      = (1 + ε / K) - c₁ * ε / ((c₀ + ε) * K) := by
    field_simp; ring
  rw [hsplit]; linarith

/-- Composite: four bounds ⟹ region, with σ = 1 + ε/K. -/
theorem region_of_bounds {c₀ c₁ ε β T₀ T₁ R A B D : ℝ}
    (hc₀ : 0 < c₀) (hc₁ : 0 < c₁) (hε : 0 < ε)
    (hA : 0 ≤ A) (hB : 0 ≤ B) (hD : 0 ≤ D)
    (hKpos : 0 < c₀ * A + c₁ * B + D)
    (hβ : β < 1 + ε / (c₀ * A + c₁ * B + D))
    (hpos : 0 ≤ c₀ * T₀ + c₁ * T₁ + R)
    (hpole : T₀ ≤ 1 / ((1 + ε / (c₀ * A + c₁ * B + D)) - 1) + A)
    (hzero : T₁ ≤ B - 1 / ((1 + ε / (c₀ * A + c₁ * B + D)) - β))
    (hR : R ≤ D) :
    β ≤ 1 - ε * (c₁ - c₀ - ε) / ((c₀ + ε) * (c₀ * A + c₁ * B + D)) := by
  have hσ : 1 < 1 + ε / (c₀ * A + c₁ * B + D) := by
    have : 0 < ε / (c₀ * A + c₁ * B + D) := div_pos hε hKpos
    linarith
  have g := gap_of_bounds hc₀ hc₁ hσ hβ hA hB hD hpos hpole hzero hR
  refine region_of_gap hc₀ hε hKpos ?_
  simpa using g

end TrigPolyGap

#print axioms TrigPolyGap.gap_of_bounds
#print axioms TrigPolyGap.region_of_gap
#print axioms TrigPolyGap.region_of_bounds
