/-
TrigPolyConstCompare.lean (HP RH center, LEDGER S2924, 2026-11-08) — lane 1.
Normalization of S2922. In the classical inputs, B ≈ a·log t and D ≈ (Σ_{k≥2} c_k)·a·log t (A = O(1)), so
K ≈ a·log t·(S − c₀), S = Σ_k c_k, and the region reads β ≤ 1 − c*(P,ε)/(a log t) with
  c*(P,ε) = ε(c₁ − c₀ − ε) / ((c₀ + ε)(S − c₀)).
Certified comparison:
  3-4-1   (c₀=3, c₁=4, S=8):  c*(ε) ≤ 1/69 for EVERY ε (numerical max 0.014359 at ε≈0.4641);
  5-8-4-1 (c₀=5, c₁=8, S=18): c*(4/3) = 60/2223 ≥ 1/40   (numerical max 0.026991 at ε≈1.3246).
⟹ the 5-8-4-1 polynomial gives a strictly larger normalized zero-free constant (ratio ≥ 69/40 = 1.725).
-/
import Mathlib

namespace TrigPolyConst

/-- normalized constant -/
noncomputable def cstar (c₀ c₁ S ε : ℝ) : ℝ := ε * (c₁ - c₀ - ε) / ((c₀ + ε) * (S - c₀))

theorem cstar_341_le (ε : ℝ) (hε : 0 < ε) : cstar 3 4 8 ε ≤ 1 / 69 := by
  unfold cstar
  rw [div_le_div_iff₀ (by nlinarith) (by norm_num)]
  nlinarith [sq_nonneg (ε - 32 / 69)]

theorem cstar_5841_val : cstar 5 8 18 (4 / 3) = 60 / 2223 := by
  unfold cstar; norm_num

theorem cstar_5841_ge : 1 / 40 ≤ cstar 5 8 18 (4 / 3) := by
  rw [cstar_5841_val]; norm_num

/-- **5-8-4-1 strictly beats 3-4-1** (normalized constant, any ε for 3-4-1). -/
theorem p5841_beats_341 (ε : ℝ) (hε : 0 < ε) : cstar 3 4 8 ε < cstar 5 8 18 (4 / 3) := by
  have := cstar_341_le ε hε
  have h2 := cstar_5841_ge
  linarith [show (1:ℝ)/69 < 1/40 by norm_num]

end TrigPolyConst

#print axioms TrigPolyConst.cstar_341_le
#print axioms TrigPolyConst.cstar_5841_ge
#print axioms TrigPolyConst.p5841_beats_341
