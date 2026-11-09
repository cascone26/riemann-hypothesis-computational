import CCMSmallMargin

/-!
# S2960 — near-resonance LOCATES an ordinate on the lattice (2π/L)ℤ (lane 2)

S2958 turned a small CCM zero-side margin into `|sin(γ_j L/2)| < δ`.  This brick converts that into a metric
statement: via Jordan's inequality `|sin x| ≥ (2/π)|x|` on `|x| ≤ π/2` and reduction mod 2π,
* `abs_sin_ge_jordan` — `|x| ≤ π/2 → (2/π)|x| ≤ |sin x|`.
* `near_lattice_of_abs_sin_lt` — `|sin(y/2)| < δ → ∃ k:ℤ, |y − 2πk| < πδ`.
* `exists_near_lattice_of_margin` — margin `< μ` (relative to CCᵀ) ⇒ `∃ j k, |γ_j L − 2πk| < π√(μL)/2`.
So the observed margin law μ_λ ~ e^{−κλ²} pins some ordinate within `O(e^{−κλ²/2})` of the resonance lattice.
-/

open Matrix Finset Real

namespace CCMLattice

theorem abs_sin_ge_jordan {x : ℝ} (hx : |x| ≤ π / 2) : 2 / π * |x| ≤ |Real.sin x| := by
  rcases le_total 0 x with h | h
  · rw [abs_of_nonneg h] at hx ⊢
    exact (Real.mul_le_sin h hx).trans (le_abs_self _)
  · have h' : 0 ≤ -x := by linarith
    rw [abs_of_nonpos h] at hx ⊢
    have := Real.mul_le_sin h' hx
    rw [Real.sin_neg] at this
    exact this.trans (neg_le_abs _)

theorem near_lattice_of_abs_sin_lt {y δ : ℝ} (h : |Real.sin (y / 2)| < δ) :
    ∃ k : ℤ, |y - 2 * π * k| < π * δ := by
  refine ⟨round (y / (2 * π)), ?_⟩
  set k : ℤ := round (y / (2 * π)) with hk
  have hpi := Real.pi_pos
  have hround : |y / (2 * π) - k| ≤ 1 / 2 := abs_sub_round _
  have heq : y - 2 * π * k = (2 * π) * (y / (2 * π) - k) := by
    field_simp
  have hr' : |y - 2 * π * k| ≤ π := by
    rw [heq, abs_mul, abs_of_pos (by positivity)]
    nlinarith [abs_nonneg (y / (2 * π) - k)]
  set r : ℝ := y - 2 * π * k with hr
  have hy : y / 2 = r / 2 + k * π := by rw [hr]; ring
  have hs : |Real.sin (y / 2)| = |Real.sin (r / 2)| := by
    rw [hy, Real.sin_add_int_mul_pi, abs_mul]; simp
  have hr2 : |r / 2| ≤ π / 2 := by
    rw [abs_div, abs_two]; linarith
  have hj := abs_sin_ge_jordan hr2
  rw [abs_div, abs_two] at hj
  have hrw : 2 / π * (|r| / 2) = |r| / π := by field_simp
  have hlt : |r| / π < δ := by rw [← hrw]; linarith
  rw [div_lt_iff₀ hpi] at hlt
  linarith

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq n] [DecidableEq m]

theorem exists_near_lattice_of_margin {L μ : ℝ} (hL : 0 < L) (hμ : 0 ≤ μ) (C : Matrix n m ℝ)
    (γ : m → ℝ) (x : n → ℝ)
    (hsmall : x ⬝ᵥ ((C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ) *ᵥ x)
      < μ * (x ⬝ᵥ ((C * Cᵀ) *ᵥ x))) :
    ∃ j, ∃ k : ℤ, |γ j * L - 2 * π * k| < π * (Real.sqrt (μ * L) / 2) := by
  obtain ⟨j, hj⟩ := CCMSmall.exists_near_resonant_sqrt hL hμ C γ x hsmall
  obtain ⟨k, hk⟩ := near_lattice_of_abs_sin_lt hj
  exact ⟨j, k, hk⟩

end CCMLattice

#print axioms CCMLattice.abs_sin_ge_jordan
#print axioms CCMLattice.near_lattice_of_abs_sin_lt
#print axioms CCMLattice.exists_near_lattice_of_margin
