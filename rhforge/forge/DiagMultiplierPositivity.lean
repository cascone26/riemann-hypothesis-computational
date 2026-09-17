import Mathlib.Tactic

/-!
# §1759 — Diagonal-multiplier form positivity (the Fourier-diagonal criterion)

Cycle 14–15 recast the localized Weil form via Suzuki (4.6) + Wiener–Khinchin into a form that is DIAGONAL
in the Fourier variable: `q(w) = ∫ M(z) |ŵ(z)|² dz` with the localized Weil multiplier
`M(z;a) = log|z| + C₀ − log a − (2A+1) − 2 Σ_{n≤e^{2a}} Λ(n)/√n cos((log n/a) z)`.  For such a diagonal
form the positivity theory is elementary and exact — this brick states it in the discrete/finite model that
the numerics actually use (a multiplier `M : ι → ℝ` acting on mode-amplitudes `x : ι → ℝ`):

`q(x) = ∑ i, M i · (x i)²`.

Then: (1) `M ≥ 0` pointwise ⟹ `q` PSD; (2) the radical is exactly the modes with `M i = 0` (support of `x`
inside `{M = 0}`); (3) `M > 0` pointwise ⟹ `q` positive-definite.  This makes the Fourier-diagonal
positivity criterion formal: `λa ≥ 0` reduces to `M(·;a) ≥ 0` on the admissible frequency support — the
exact statement the symbol/uncertainty analysis (cycle 15) must establish.  Generalizes the Gram/radical
brick §1756 (there the radical was `ker L`; here it is `supp x ⊆ {M=0}`).  EXIT0 / std-axioms.
-/

open Finset

variable {ι : Type*} [Fintype ι]

/-- The diagonal multiplier form `q(x) = ∑ i, M i · (x i)²`. -/
def diagForm (M x : ι → ℝ) : ℝ := ∑ i, M i * (x i * x i)

/-- **PSD from a nonnegative multiplier.** -/
theorem diagForm_nonneg (M x : ι → ℝ) (hM : ∀ i, 0 ≤ M i) : 0 ≤ diagForm M x := by
  apply Finset.sum_nonneg
  intro i _
  exact mul_nonneg (hM i) (mul_self_nonneg (x i))

/-- **Radical = modes on `{M = 0}`.**  With a nonnegative multiplier the form vanishes iff every mode's
contribution vanishes (each amplitude sits where `M = 0` or is itself `0`). -/
theorem diagForm_eq_zero_iff (M x : ι → ℝ) (hM : ∀ i, 0 ≤ M i) :
    diagForm M x = 0 ↔ ∀ i, M i * (x i * x i) = 0 := by
  unfold diagForm
  rw [Finset.sum_eq_zero_iff_of_nonneg (fun i _ => mul_nonneg (hM i) (mul_self_nonneg (x i)))]
  simp

/-- **Positive-definiteness from a strictly positive multiplier.** -/
theorem diagForm_posDef (M x : ι → ℝ) (hM : ∀ i, 0 < M i)
    (h : diagForm M x = 0) : x = 0 := by
  have hnn : ∀ i, 0 ≤ M i := fun i => (hM i).le
  rw [diagForm_eq_zero_iff M x hnn] at h
  funext i
  simp only [Pi.zero_apply]
  have hi := h i
  have hxx : x i * x i = 0 := by
    rcases mul_eq_zero.mp hi with h1 | h2
    · exact absurd h1 (ne_of_gt (hM i))
    · exact h2
  rcases mul_eq_zero.mp hxx with h1 | h2
  · exact h1
  · exact h2
