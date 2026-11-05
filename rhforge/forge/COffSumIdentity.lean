import EOffInterpVector
import Mathlib.LinearAlgebra.Lagrange

open scoped MatrixOrder Matrix Polynomial
open Polynomial

/-!
# Off-critical-line lane: Lagrange interpolation identity `∑ c^{off}_n = 1`

Theorem: the interpolation coefficients `c_off K L z n` sum to 1 over all frequency nodes.
Proved via Mathlib's `Lagrange.leadingCoeff_eq_sum`: if a degree-exactly-`(#s - 1)` polynomial `P`
is evaluated at `#s` distinct nodes `v i`, its leading coefficient equals
`∑ i, P.eval (v i) / ∏ j≠i, (v i - v j)`. Taking `P := ` the `Polynomial ℝ` lift of `Q_off_K`
(monic of degree `4K`) at the `4K+1` frequency nodes gives exactly `∑ c_off[n] = Q_off_K.leadingCoeff = 1`.
-/

-- Abbreviation for the frequency node at index n.
noncomputable def freq_node (K : ℕ) (L : ℝ) (n : Fin (2 * (2 * K) + 1)) : ℝ :=
  2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L

-- Helper: frequency nodes are pairwise distinct when L ≠ 0.
lemma freq_nodes_distinct (K : ℕ) (L : ℝ) (hL : L ≠ 0) (m n : Fin (2 * (2 * K) + 1))
    (hmn : m ≠ n) : freq_node K L m ≠ freq_node K L n := by
  unfold freq_node
  intro h
  have h_eq : 2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L =
              2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L := h
  have pi_ne : Real.pi ≠ 0 := Real.pi_ne_zero
  have two_pi_ne : (2 : ℝ) * Real.pi ≠ 0 := mul_ne_zero (by norm_num) pi_ne
  have h1 : (2 : ℝ) * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) =
            (2 : ℝ) * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) := by
    have : 2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L * L =
            2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L * L := by
      rw [h_eq]
    simp only [div_mul_cancel₀ _ hL] at this
    exact this
  have h2 : ((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ) = ((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ) :=
    mul_left_cancel₀ two_pi_ne h1
  have nat_eq : (m : ℕ) = (n : ℕ) := by
    have h_cast : ((m : ℕ) : ℝ) = ((n : ℕ) : ℝ) := by linarith [h2]
    exact Nat.cast_injective h_cast
  exact hmn (Fin.ext nat_eq)

-- Helper: the product ∏_{j ≠ n} (ν_n - ν_j) is nonzero.
lemma lagrange_denom_nonzero (K : ℕ) (L : ℝ) (hL : L ≠ 0) (n : Fin (2 * (2 * K) + 1)) :
    (∏ j ∈ Finset.univ.erase n, (freq_node K L n - freq_node K L j)) ≠ 0 := by
  rw [Finset.prod_ne_zero_iff]
  intro j hj
  have hne : (j : Fin (2 * (2 * K) + 1)) ≠ n := by
    simp only [Finset.mem_erase] at hj
    exact hj.1
  have : freq_node K L n ≠ freq_node K L j :=
    freq_nodes_distinct K L hL n j (Ne.symm hne)
  exact sub_ne_zero.mpr this

-- The `Polynomial ℝ` lift of a single quartic factor of `Q_off_K`.
noncomputable def Q_off_K_factor_poly (a b : ℝ) : ℝ[X] :=
  X ^ 4 - C (2 * a) * X ^ 2 + C b

lemma Q_off_K_factor_poly_monic (a b : ℝ) : (Q_off_K_factor_poly a b).Monic := by
  unfold Q_off_K_factor_poly
  monicity!

lemma Q_off_K_factor_poly_natDegree (a b : ℝ) : (Q_off_K_factor_poly a b).natDegree = 4 := by
  unfold Q_off_K_factor_poly
  compute_degree!

lemma Q_off_K_factor_poly_eval (a b w : ℝ) :
    (Q_off_K_factor_poly a b).eval w = w ^ 4 - 2 * a * w ^ 2 + b := by
  unfold Q_off_K_factor_poly
  simp

-- The `Polynomial ℝ` lift of `Q_off_K`.
noncomputable def Q_off_K_poly (K : ℕ) (z : Fin K → ℂ) : ℝ[X] :=
  ∏ k : Fin K, Q_off_K_factor_poly ((z k).re ^ 2 - (z k).im ^ 2) (((z k).re ^ 2 + (z k).im ^ 2) ^ 2)

theorem Q_off_K_poly_eval (K : ℕ) (z : Fin K → ℂ) (w : ℝ) :
    (Q_off_K_poly K z).eval w = Q_off_K K z w := by
  unfold Q_off_K_poly Q_off_K
  rw [eval_prod]
  apply Finset.prod_congr rfl
  intro k _
  rw [Q_off_K_factor_poly_eval]

theorem Q_off_K_poly_monic (K : ℕ) (z : Fin K → ℂ) : (Q_off_K_poly K z).Monic := by
  unfold Q_off_K_poly
  exact Polynomial.monic_prod_of_monic _ _ (fun k _ => Q_off_K_factor_poly_monic _ _)

theorem Q_off_K_poly_natDegree (K : ℕ) (z : Fin K → ℂ) :
    (Q_off_K_poly K z).natDegree = 4 * K := by
  unfold Q_off_K_poly
  rw [Polynomial.natDegree_prod_of_monic _ _ (fun k _ => Q_off_K_factor_poly_monic _ _)]
  simp only [Q_off_K_factor_poly_natDegree]
  rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin]
  ring

-- Main theorem: the coefficients sum to 1.
theorem c_off_sum_eq_one (K : ℕ) (L : ℝ) (hL : L ≠ 0) (z : Fin K → ℂ) :
    ∑ n : Fin (2 * (2 * K) + 1), c_off K L z n = 1 := by
  have hmonic := Q_off_K_poly_monic K z
  have hdeg : (Q_off_K_poly K z).natDegree = 4 * K := Q_off_K_poly_natDegree K z
  have hInj : Set.InjOn (freq_node K L) (Finset.univ : Finset (Fin (2 * (2 * K) + 1))) := by
    intro m _ n _ h
    by_contra hmn
    exact freq_nodes_distinct K L hL m n hmn h
  have hdegree_eq : (Q_off_K_poly K z).degree = (4 * K : ℕ) := by
    rw [Polynomial.degree_eq_natDegree hmonic.ne_zero, hdeg]
  have hcard : (Finset.univ : Finset (Fin (2 * (2 * K) + 1))).card
      = (Q_off_K_poly K z).degree + 1 := by
    rw [Finset.card_univ, Fintype.card_fin, hdegree_eq]
    norm_cast
    ring
  have hlead := Lagrange.leadingCoeff_eq_sum hInj hcard
  rw [hmonic.leadingCoeff] at hlead
  have hrw : ∀ n : Fin (2 * (2 * K) + 1),
      (Q_off_K_poly K z).eval (freq_node K L n)
        / ∏ j ∈ Finset.univ.erase n, (freq_node K L n - freq_node K L j)
      = c_off K L z n := by
    intro n
    rw [Q_off_K_poly_eval]
    rw [c_off_eq_div_prod]
    rfl
  rw [Finset.sum_congr rfl (fun n _ => hrw n)] at hlead
  exact hlead.symm

#print axioms c_off_sum_eq_one
