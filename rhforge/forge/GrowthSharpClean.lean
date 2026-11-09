import GrowthSharpN

/-!
# S2963 — clean scalar growth with leading coefficient 8, uniform in R (lane 1)

Specialises S2961 to `N = ⌊|t|⌋₊`.  With `x = |t| ≥ 3`, `a = σ-R`, `p = 1-a = 1-σ+R`:
  `sharpMs + 1 ≤ x^p·(1/p + 3 + 4/a)`   (terms: x^p/p, 1, x^p/(x-R) ≤ x^p, (σ+x+R)N^{-a}/a ≤ 4x^p/a, +1),
and since `p ≤ R` the Borel–Carathéodory factor `4/(R/2)` turns `p·log x` into at most `8·log x`:

* `sharpMs_floor_add_one_le` — the scalar bound above.
* `growth_at_sharp_clean` — `(-ζ'/ζ(σ+it)).re ≤ 8·log|t| + (8/R)·(log(1/(1-σ+R)+3+4/(σ-R)) + log(3+1/(σ-1)))`
  for `1 < σ ≤ 2`, `0 < R < 1`, `σ - R < 1`, `3 ≤ |t|`.
The N = 1 route gave `(8/R)·log|t|` (32/3 at R=3/4, 128/15 at R=15/16).  Leading coefficient now 8 for
every admissible R; R only enters the O(1) term.  Next: substitute into the 5-8-4-1 region for the new c₀ rung.
-/

open Complex Metric

namespace GrowthSharpClean

theorem sharpMs_floor_add_one_le {σ t R : ℝ} (hσ : 1 < σ) (hσ2 : σ ≤ 2) (hR : 0 < R)
    (hR1 : R < 1) (haR : σ - R < 1) (ht : 3 ≤ |t|) :
    GrowthSharp.sharpMs σ t R ⌊|t|⌋₊ + 1
      ≤ |t| ^ (1 - (σ - R)) * (1 / (1 - (σ - R)) + 3 + 4 / (σ - R)) := by
  set x : ℝ := |t| with hxdef
  set a : ℝ := σ - R with hadef
  set p : ℝ := 1 - a with hpdef
  set N : ℕ := ⌊x⌋₊ with hNdef
  have hx0 : 0 ≤ x := by linarith
  have hapos : 0 < a := by rw [hadef]; linarith
  have hppos : 0 < p := by rw [hpdef]; linarith
  have hNx : (N : ℝ) ≤ x := Nat.floor_le hx0
  have hNlt : x < (N : ℝ) + 1 := Nat.lt_floor_add_one x
  have hNpos : (0 : ℝ) < (N : ℝ) := by linarith
  have hN2 : x ≤ 2 * (N : ℝ) := by linarith
  have hX1 : 1 ≤ x ^ p := Real.one_le_rpow (by linarith) hppos.le
  have hNp : (N : ℝ) ^ p ≤ x ^ p := Real.rpow_le_rpow hNpos.le hNx hppos.le
  have hNp1 : 1 ≤ (N : ℝ) ^ p := Real.one_le_rpow (by linarith) hppos.le
  -- T1
  have hT1 : ((N : ℝ) ^ (1 - a) - 1) / (1 - a) ≤ x ^ p / p := by
    rw [← hpdef]
    exact div_le_div_of_nonneg_right (by linarith) hppos.le
  -- T3
  have hT3 : (N : ℝ) ^ (1 - a) / (x - R) ≤ x ^ p := by
    rw [← hpdef]
    exact le_trans (div_le_self (by linarith) (by linarith)) hNp
  -- T4
  have hNa : (N : ℝ) * (N : ℝ) ^ (-a) = (N : ℝ) ^ p := by
    rw [hpdef, sub_eq_add_neg, Real.rpow_add hNpos, Real.rpow_one]
  have hNa0 : 0 ≤ (N : ℝ) ^ (-a) := Real.rpow_nonneg hNpos.le _
  have hT4 : (σ + x + R) * ((N : ℝ) ^ (-a) / a) ≤ 4 * x ^ p / a := by
    have h1 : (σ + x + R) * (N : ℝ) ^ (-a) ≤ 4 * (N : ℝ) ^ p := by
      rw [← hNa]
      have : σ + x + R ≤ 4 * (N : ℝ) := by linarith
      nlinarith
    have h2 : (σ + x + R) * ((N : ℝ) ^ (-a) / a) = (σ + x + R) * (N : ℝ) ^ (-a) / a := by ring
    rw [h2]
    exact div_le_div_of_nonneg_right (by linarith) hapos.le
  unfold GrowthSharp.sharpMs
  rw [← hadef, ← hxdef]
  have hexp : x ^ p * (1 / p + 3 + 4 / a) = x ^ p / p + 3 * x ^ p + 4 * x ^ p / a := by ring
  rw [hexp, ← hpdef]
  linarith

theorem growth_at_sharp_clean {σ t R : ℝ} (hσ : 1 < σ) (hσ2 : σ ≤ 2) (hR : 0 < R)
    (hR1 : R < 1) (haR : σ - R < 1) (ht : 3 ≤ |t|) :
    (-logDeriv riemannZeta ((σ : ℂ) + I * (t : ℂ))).re
      ≤ 8 * Real.log |t| + 8 / R * (Real.log (1 / (1 - (σ - R)) + 3 + 4 / (σ - R))
          + Real.log (3 + 1 / (σ - 1))) := by
  have hN : 1 ≤ ⌊|t|⌋₊ := Nat.le_floor (by norm_num; linarith)
  have g := GrowthSharp.growth_at_sharpN (R' := |t|) hN hσ hR hR1 haR (by linarith) le_rfl
  refine g.trans ?_
  set x : ℝ := |t| with hxdef
  set a : ℝ := σ - R with hadef
  set p : ℝ := 1 - a with hpdef
  set M : ℝ := GrowthSharp.sharpMs σ t R ⌊|t|⌋₊ with hMdef
  set K : ℝ := 1 / p + 3 + 4 / a with hKdef
  set D : ℝ := 3 + 1 / (σ - 1) with hDdef
  have hapos : 0 < a := by rw [hadef]; linarith
  have hppos : 0 < p := by rw [hpdef]; linarith
  have hpR : p ≤ R := by rw [hpdef, hadef]; linarith
  have hx0 : 0 < x := by linarith
  have hlogx : 0 ≤ Real.log x := Real.log_nonneg (by linarith)
  have hK : 0 < K := by rw [hKdef]; positivity
  have hD : 0 < D := by
    rw [hDdef]; have : 0 < 1 / (σ - 1) := by apply div_pos one_pos; linarith
    linarith
  have hX : 0 < x ^ p := Real.rpow_pos_of_pos hx0 p
  have hM1 : M + 1 ≤ x ^ p * K := by
    rw [hMdef, hKdef, hpdef, hadef, hxdef]
    exact sharpMs_floor_add_one_le hσ hσ2 hR hR1 haR ht
  have hMpos : 0 < M + 1 := by
    have hN0 : (1 : ℝ) ≤ (⌊|t|⌋₊ : ℝ) := by exact_mod_cast hN
    have h1 : 0 ≤ ((⌊|t|⌋₊ : ℝ) ^ (1 - (σ - R)) - 1) / (1 - (σ - R)) := by
      apply div_nonneg _ (by linarith)
      have := Real.one_le_rpow hN0 (show 0 ≤ 1 - (σ - R) by linarith)
      linarith
    have h3 : 0 ≤ (⌊|t|⌋₊ : ℝ) ^ (1 - (σ - R)) / (|t| - R) :=
      div_nonneg (Real.rpow_nonneg (by linarith) _) (by linarith)
    have h4 : 0 ≤ (σ + |t| + R) * ((⌊|t|⌋₊ : ℝ) ^ (-(σ - R)) / (σ - R)) :=
      mul_nonneg (by linarith [abs_nonneg t])
        (div_nonneg (Real.rpow_nonneg (by linarith) _) (by linarith))
    rw [hMdef]; unfold GrowthSharp.sharpMs; linarith
  have hlog : Real.log ((M + 1) * D) ≤ p * Real.log x + Real.log K + Real.log D := by
    calc Real.log ((M + 1) * D)
        ≤ Real.log ((x ^ p * K) * D) :=
          Real.log_le_log (mul_pos hMpos hD) (mul_le_mul_of_nonneg_right hM1 hD.le)
      _ = p * Real.log x + Real.log K + Real.log D := by
          rw [Real.log_mul (by positivity) hD.ne', Real.log_mul hX.ne' hK.ne',
            Real.log_rpow hx0]
  have hsplit : 4 * Real.log ((M + 1) * D) / (R / 2) = 8 / R * Real.log ((M + 1) * D) := by
    rw [div_div_eq_mul_div]; ring
  rw [hsplit]
  have h8R : 0 ≤ 8 / R := by positivity
  have step1 : 8 / R * Real.log ((M + 1) * D)
      ≤ 8 / R * (p * Real.log x + Real.log K + Real.log D) :=
    mul_le_mul_of_nonneg_left hlog h8R
  have step2 : 8 / R * (p * Real.log x) ≤ 8 * Real.log x := by
    have hpR' : p / R ≤ 1 := (div_le_one hR).mpr hpR
    have : 8 / R * (p * Real.log x) = 8 * Real.log x * (p / R) := by ring
    rw [this]
    nlinarith
  have hexp : 8 / R * (p * Real.log x + Real.log K + Real.log D)
      = 8 / R * (p * Real.log x) + 8 / R * (Real.log K + Real.log D) := by ring
  linarith

end GrowthSharpClean

#print axioms GrowthSharpClean.sharpMs_floor_add_one_le
#print axioms GrowthSharpClean.growth_at_sharp_clean
