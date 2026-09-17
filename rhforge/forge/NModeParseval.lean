import EigenbasisOrthonormal
import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic

/-!
# §1737 — N-mode Parseval: full positive-definiteness of the archimedean metric

The two-mode identity §1735 generalises to any finite family of interval eigenmodes.  With frequencies
`ωᵢ·L = kᵢ·π`, pairwise `ωᵢ ≠ ±ωⱼ` for `i ≠ j`, the archimedean `L²` self-energy of an arbitrary finite
combination is the plain sum of squares (Parseval):

`∫₀ᴸ (∑ᵢ cᵢ·sin(ωᵢ u))² du = (∑ᵢ cᵢ²)·(L/2)`.

All cross terms vanish by orthogonality (§1734), each square contributes `L/2` (normalization).  Hence the
archimedean Gram matrix is `(L/2)·I` on the eigenbasis — **positive definite** for `L>0`.  This is the full
(any-dimension) grounding of the positive-definite inner product the Rosati critical-circle brick §1729 takes
as hypothesis.  EXIT0 / std-axioms.
-/

open Real intervalIntegral MeasureTheory

/-- **N-mode Parseval.**  `∫₀ᴸ (∑ᵢ cᵢ sin(ωᵢ u))² = (∑ᵢ cᵢ²)(L/2)`. -/
theorem n_mode_parseval {ι : Type*} [DecidableEq ι] (s : Finset ι) (ω c : ι → ℝ) (L : ℝ) (kk : ι → ℤ)
    (hω : ∀ i ∈ s, ω i ≠ 0) (hqu : ∀ i ∈ s, ω i * L = (kk i : ℝ) * π)
    (hdist : ∀ i ∈ s, ∀ j ∈ s, i ≠ j → ω i - ω j ≠ 0 ∧ ω i + ω j ≠ 0) :
    (∫ u in (0:ℝ)..L, (∑ i ∈ s, c i * Real.sin (ω i * u)) ^ 2)
      = (∑ i ∈ s, (c i) ^ 2) * (L / 2) := by
  have integ : ∀ i j : ι, IntervalIntegrable
      (fun u => (c i * Real.sin (ω i * u)) * (c j * Real.sin (ω j * u))) volume 0 L :=
    fun i j => (by fun_prop : Continuous _).intervalIntegrable 0 L
  have integ_inner : ∀ i ∈ s, IntervalIntegrable
      (fun u => ∑ j ∈ s, (c i * Real.sin (ω i * u)) * (c j * Real.sin (ω j * u))) volume 0 L :=
    fun i _ => (by fun_prop : Continuous _).intervalIntegrable 0 L
  -- per-pair integral value
  have key : ∀ i ∈ s, ∀ j ∈ s,
      (∫ u in (0:ℝ)..L, (c i * Real.sin (ω i * u)) * (c j * Real.sin (ω j * u)))
        = (if i = j then (c i) ^ 2 * (L / 2) else 0) := by
    intro i hi j hj
    by_cases hij : i = j
    · subst hij
      rw [intervalIntegral.integral_congr
            (g := fun u => (c i) ^ 2 * (Real.sin (ω i * u) * Real.sin (ω i * u)))
            (fun u _ => by ring),
          intervalIntegral.integral_const_mul, eigenmode_norm (ω i) L (kk i) (hω i hi) (hqu i hi)]
      simp
    · obtain ⟨hd1, hs1⟩ := hdist i hi j hj hij
      rw [intervalIntegral.integral_congr
            (g := fun u => (c i * c j) * (Real.sin (ω i * u) * Real.sin (ω j * u)))
            (fun u _ => by ring),
          intervalIntegral.integral_const_mul,
          eigenbasis_orthogonal (ω i) (ω j) L (kk i) (kk j) hd1 hs1 (hqu i hi) (hqu j hj)]
      simp [hij]
  calc (∫ u in (0:ℝ)..L, (∑ i ∈ s, c i * Real.sin (ω i * u)) ^ 2)
      = ∫ u in (0:ℝ)..L, ∑ i ∈ s, ∑ j ∈ s,
          (c i * Real.sin (ω i * u)) * (c j * Real.sin (ω j * u)) := by
        refine intervalIntegral.integral_congr (fun u _ => ?_)
        rw [sq, Finset.sum_mul_sum]
    _ = ∑ i ∈ s, ∫ u in (0:ℝ)..L, ∑ j ∈ s,
          (c i * Real.sin (ω i * u)) * (c j * Real.sin (ω j * u)) :=
        intervalIntegral.integral_finsetSum integ_inner
    _ = ∑ i ∈ s, ∑ j ∈ s, ∫ u in (0:ℝ)..L,
          (c i * Real.sin (ω i * u)) * (c j * Real.sin (ω j * u)) :=
        Finset.sum_congr rfl (fun i _ => intervalIntegral.integral_finsetSum (fun j _ => integ i j))
    _ = ∑ i ∈ s, ∑ j ∈ s, (if i = j then (c i) ^ 2 * (L / 2) else 0) :=
        Finset.sum_congr rfl (fun i hi => Finset.sum_congr rfl (fun j hj => key i hi j hj))
    _ = ∑ i ∈ s, (c i) ^ 2 * (L / 2) := by
        refine Finset.sum_congr rfl (fun i hi => ?_)
        rw [Finset.sum_ite_eq s i (fun _ => (c i) ^ 2 * (L / 2))]
        simp [hi]
    _ = (∑ i ∈ s, (c i) ^ 2) * (L / 2) := by rw [Finset.sum_mul]
