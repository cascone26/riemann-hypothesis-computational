import RiemannLean.Bessel
open scoped Nat
/-- Upper bound `J₁(x) ≤ x/2` for `0 < x`, `x² < 8`: since `J₁(x)=(x/2)·S` with the alternating
series `S ≤ S₀ = 1`.  With `besselJ1_pos`, completes the J₁ characterization on `(0,2√2)`. -/
theorem besselJ1_le_half {x : ℝ} (hx0 : 0 < x) (hx : x ^ 2 < 8) : besselJ1 x ≤ x / 2 := by
  set y : ℝ := (x / 2) ^ 2 with hy
  have hy0 : 0 ≤ y := by rw [hy]; positivity
  have hy2 : y < 2 := by rw [hy]; nlinarith [hx]
  set f : ℕ → ℝ := fun k => y ^ k / ((k ! : ℝ) * ((k + 1)! : ℝ)) with hf
  have hfnn : ∀ k, 0 ≤ f k := fun k => by rw [hf]; positivity
  have hxhalf : (0 : ℝ) < x / 2 := by linarith
  have hfactor : besselJ1 x = (x / 2) * ∑' k, (-1 : ℝ) ^ k * f k := by
    unfold besselJ1; rw [← tsum_mul_left]; congr 1; ext k
    rw [hf, hy, pow_add, pow_mul, pow_one]; ring
  have hfa : Antitone f := by
    apply antitone_nat_of_succ_le; intro k
    have hfac1 : ((k + 1)! : ℝ) = ((k : ℝ) + 1) * (k ! : ℝ) := by
      rw [Nat.factorial_succ]; push_cast; ring
    have hfac2 : ((k + 2)! : ℝ) = ((k : ℝ) + 2) * ((k + 1)! : ℝ) := by
      rw [show k + 2 = (k + 1) + 1 from rfl, Nat.factorial_succ]; push_cast; ring
    have hkne : (k ! : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero k)
    have hk1ne : ((k + 1)! : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero (k + 1))
    have hstep : f (k + 1) = (y / (((k : ℝ) + 1) * ((k : ℝ) + 2))) * f k := by
      show y ^ (k + 1) / (((k + 1)! : ℝ) * ((k + 2)! : ℝ))
        = (y / (((k : ℝ) + 1) * ((k : ℝ) + 2))) * (y ^ k / ((k ! : ℝ) * ((k + 1)! : ℝ)))
      rw [hfac2, hfac1, pow_succ]; field_simp
    rw [hstep]
    have hc : y / (((k : ℝ) + 1) * ((k : ℝ) + 2)) ≤ 1 := by
      rw [div_le_one (by positivity)]; nlinarith [Nat.cast_nonneg (α := ℝ) k, hy2]
    calc (y / (((k : ℝ) + 1) * ((k : ℝ) + 2))) * f k ≤ 1 * f k :=
          mul_le_mul_of_nonneg_right hc (hfnn k)
      _ = f k := one_mul _
  have hfs : Summable f := by
    have hcomp := Real.summable_pow_div_factorial y
    apply hcomp.of_nonneg_of_le hfnn
    intro k
    have hk1 : (1 : ℝ) ≤ ((k + 1)! : ℝ) := by
      exact_mod_cast Nat.one_le_iff_ne_zero.mpr (Nat.factorial_ne_zero (k + 1))
    have hkpos : (0 : ℝ) < (k ! : ℝ) := by positivity
    show y ^ k / ((k ! : ℝ) * ((k + 1)! : ℝ)) ≤ y ^ k / (k ! : ℝ)
    gcongr; nlinarith [hk1, hkpos]
  have hb := alternating_series_error_bound f hfa hfs 0
  simp only [Finset.range_zero, Finset.sum_empty, sub_zero] at hb
  have hf0 : f 0 = 1 := by simp [hf]
  rw [abs_le] at hb
  have hS : (∑' k, (-1 : ℝ) ^ k * f k) ≤ 1 := by have h2 := hb.2; rwa [hf0] at h2
  rw [hfactor]
  calc (x / 2) * ∑' k, (-1 : ℝ) ^ k * f k ≤ (x / 2) * 1 :=
        mul_le_mul_of_nonneg_left hS (le_of_lt hxhalf)
    _ = x / 2 := mul_one _

#print axioms besselJ1_le_half
