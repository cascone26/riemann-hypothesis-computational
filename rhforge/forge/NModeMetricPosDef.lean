import NModeParseval

/-!
# §1738 — Full positive-definiteness: zero archimedean energy forces zero combination

The any-dimension conclusion from the N-mode Parseval §1737.  On a window of positive length, the archimedean
`L²` metric on the interval eigenbasis is positive definite: the only finite combination with zero self-energy
is the trivial one.

`0 < L` and `∫₀ᴸ (∑ᵢ cᵢ sin(ωᵢ u))² du = 0  ⟹  ∀ i ∈ s, cᵢ = 0`.

This completes the grounding (§1735/§1736 → any dimension) of the positive-definite inner-product hypothesis of
the Rosati critical-circle brick §1729.  EXIT0 / std-axioms.
-/

open Real intervalIntegral MeasureTheory

/-- **N-mode positive-definiteness.**  Zero self-energy forces all coefficients to vanish. -/
theorem n_mode_metric_pos_def {ι : Type*} [DecidableEq ι] (s : Finset ι) (ω c : ι → ℝ) (L : ℝ) (kk : ι → ℤ)
    (hL : 0 < L) (hω : ∀ i ∈ s, ω i ≠ 0) (hqu : ∀ i ∈ s, ω i * L = (kk i : ℝ) * π)
    (hdist : ∀ i ∈ s, ∀ j ∈ s, i ≠ j → ω i - ω j ≠ 0 ∧ ω i + ω j ≠ 0)
    (h0 : (∫ u in (0:ℝ)..L, (∑ i ∈ s, c i * Real.sin (ω i * u)) ^ 2) = 0) :
    ∀ i ∈ s, c i = 0 := by
  rw [n_mode_parseval s ω c L kk hω hqu hdist] at h0
  have hsum : ∑ i ∈ s, (c i) ^ 2 = 0 := by
    rcases mul_eq_zero.mp h0 with h | h
    · exact h
    · exact absurd h (ne_of_gt (by linarith : (0:ℝ) < L / 2))
  intro i hi
  have hi2 : (c i) ^ 2 = 0 :=
    (Finset.sum_eq_zero_iff_of_nonneg (fun j _ => sq_nonneg (c j))).mp hsum i hi
  exact sq_eq_zero_iff.mp hi2
