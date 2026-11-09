import CCMLatticeLocate

/-!
# S2962 — lattice distance ⇒ quantitative zero-side positivity (lane 2, converse of S2960)

S2960: small margin ⇒ some γ_j L near 2πℤ.  Converse here: if EVERY ordinate keeps lattice distance ≥ η,
the CCM zero-side Gram dominates a fixed multiple of the Cauchy Gram:
* `exists_k_abs_sin_half_ge` — ∀ y, ∃ k:ℤ, |y − 2πk| ≤ π ∧ |y − 2πk|/π ≤ |sin(y/2)|.
* `abs_sin_half_ge_of_lattice_gap` — (∀ k, η ≤ |y − 2πk|) ⇒ η/π ≤ |sin(y/2)|.
* `ccm_gram_ge_lattice_gap` / `ccm_rayleigh_ge_lattice_gap` — (∀ j k, η ≤ |γ_j L − 2πk|) ⇒
  G_ccm − (4η²/(π²L))·CCᵀ ⪰ 0, and the Rayleigh form.
Together with S2960: margin ≍ (lattice distance)² up to the constants π²/4 — the two-sided dictionary between
the observed μ_λ decay and resonance proximity (zero side, truncated, RH-conditional).
-/

open Matrix Finset Real

namespace CCMLatticeGap

theorem exists_k_abs_sin_half_ge (y : ℝ) :
    ∃ k : ℤ, |y - 2 * π * k| ≤ π ∧ |y - 2 * π * k| / π ≤ |Real.sin (y / 2)| := by
  refine ⟨round (y / (2 * π)), ?_⟩
  set k : ℤ := round (y / (2 * π)) with hk
  have hpi := Real.pi_pos
  have hround : |y / (2 * π) - k| ≤ 1 / 2 := abs_sub_round _
  have heq : y - 2 * π * k = (2 * π) * (y / (2 * π) - k) := by
    field_simp
  have hr' : |y - 2 * π * k| ≤ π := by
    rw [heq, abs_mul, abs_of_pos (by positivity)]
    nlinarith [abs_nonneg (y / (2 * π) - k)]
  refine ⟨hr', ?_⟩
  set r : ℝ := y - 2 * π * k with hr
  have hy : y / 2 = r / 2 + k * π := by rw [hr]; ring
  have hs : |Real.sin (y / 2)| = |Real.sin (r / 2)| := by
    rw [hy, Real.sin_add_int_mul_pi, abs_mul]; simp
  have hr2 : |r / 2| ≤ π / 2 := by
    rw [abs_div, abs_two]; linarith
  have hj := CCMLattice.abs_sin_ge_jordan hr2
  rw [abs_div, abs_two] at hj
  have hrw : 2 / π * (|r| / 2) = |r| / π := by field_simp
  rw [hs, ← hrw]; exact hj

theorem abs_sin_half_ge_of_lattice_gap {y η : ℝ} (h : ∀ k : ℤ, η ≤ |y - 2 * π * k|) :
    η / π ≤ |Real.sin (y / 2)| := by
  obtain ⟨k, -, hk⟩ := exists_k_abs_sin_half_ge y
  exact le_trans (div_le_div_of_nonneg_right (h k) Real.pi_pos.le) hk

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq n] [DecidableEq m]

theorem ccm_gram_ge_lattice_gap {L η : ℝ} (hL : 0 < L) (hη : 0 ≤ η) (C : Matrix n m ℝ)
    (γ : m → ℝ) (hgap : ∀ j, ∀ k : ℤ, η ≤ |γ j * L - 2 * π * k|) :
    (C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ
      - (4 * (η / π) ^ 2 / L) • (C * Cᵀ)).PosSemidef := by
  refine CCMGap.ccm_gram_ge_gap hL (div_nonneg hη Real.pi_pos.le) C γ (fun j => ?_)
  have := abs_sin_half_ge_of_lattice_gap (hgap j)
  rwa [show γ j * L / 2 = (γ j * L) / 2 by ring]

theorem ccm_rayleigh_ge_lattice_gap {L η : ℝ} (hL : 0 < L) (hη : 0 ≤ η) (C : Matrix n m ℝ)
    (γ : m → ℝ) (hgap : ∀ j, ∀ k : ℤ, η ≤ |γ j * L - 2 * π * k|) (x : n → ℝ) :
    4 * (η / π) ^ 2 / L * (x ⬝ᵥ ((C * Cᵀ) *ᵥ x))
      ≤ x ⬝ᵥ ((C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ) *ᵥ x) := by
  refine CCMGap.ccm_rayleigh_ge_gap hL (div_nonneg hη Real.pi_pos.le) C γ (fun j => ?_) x
  exact abs_sin_half_ge_of_lattice_gap (hgap j)

end CCMLatticeGap

#print axioms CCMLatticeGap.exists_k_abs_sin_half_ge
#print axioms CCMLatticeGap.abs_sin_half_ge_of_lattice_gap
#print axioms CCMLatticeGap.ccm_gram_ge_lattice_gap
#print axioms CCMLatticeGap.ccm_rayleigh_ge_lattice_gap
