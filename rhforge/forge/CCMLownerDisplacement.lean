import Mathlib

/-!
# S2968 — arithmetic side: W_R and W_P are Löwner matrices with rank-2 displacement (lane 2)

From `operators/ccm_parts.py`, both the archimedean block and the prime block have the divided-difference form
  `M[i,j] = d_i` (i = j),   `M[i,j] = (g_j − g_i)/(π(ν_i − ν_j))` (i ≠ j)
with ν = frequency n, and
* W_R : `g = Jsig` (odd sine-transform of ρ),       `d = archconst + D_diag`;
* W_P : `g_n = Σ_k Λ(k) k^{-1/2} sin(2πn log k / L)`, `d_n = Σ_k Λ(k) k^{-1/2} 2(L − log k) cos(2πn log k/L)/L`.
Theorems (generic in d, g, ν):
* `lowner_symm` — `Mᵀ = M`.
* `lowner_add` — `lowner d g ν + lowner d' g' ν = lowner (d + d') (g + g') ν`  (so W_R + W_P is one Löwner matrix).
* `lowner_displacement` — for injective ν:  `diag(ν)·M − M·diag(ν) = (1/π)•(𝟙 gᵀ − g 𝟙ᵀ)`  (rank ≤ 2),
  independent of the diagonal d.
Numeric cross-check (python parts(20,8)): commutators [diag(n), W_R], [diag(n), W_P], [diag(n), W02] all rank 2.
Reading (Report): the whole arithmetic matrix A = W02 − W_R − W_P has low displacement rank w.r.t. diag(n) —
the Cauchy/Löwner structure that made the zero-side bricks tractable is present on the arithmetic side too.
-/

open Matrix Finset Real

namespace CCMLowner

variable {ι : Type*} [DecidableEq ι]

noncomputable def lowner (d g ν : ι → ℝ) : Matrix ι ι ℝ := fun i j =>
  if i = j then d i else (g j - g i) / (π * (ν i - ν j))

theorem lowner_symm (d g ν : ι → ℝ) : (lowner d g ν)ᵀ = lowner d g ν := by
  ext i j
  simp only [transpose_apply, lowner]
  by_cases h : i = j
  · subst h; simp
  · have h' : j ≠ i := fun e => h e.symm
    rw [if_neg h', if_neg h]
    rw [show g i - g j = -(g j - g i) by ring, show ν j - ν i = -(ν i - ν j) by ring,
      mul_neg, neg_div_neg_eq]

theorem lowner_add (d g d' g' ν : ι → ℝ) :
    lowner d g ν + lowner d' g' ν = lowner (d + d') (g + g') ν := by
  ext i j
  simp only [Matrix.add_apply, lowner, Pi.add_apply]
  by_cases h : i = j
  · simp [h]
  · rw [if_neg h, if_neg h, if_neg h, ← add_div]; ring_nf

variable [Fintype ι]

theorem lowner_displacement (d g : ι → ℝ) {ν : ι → ℝ} (hν : Function.Injective ν) :
    diagonal ν * lowner d g ν - lowner d g ν * diagonal ν
      = (1 / π) • (vecMulVec (fun _ => (1 : ℝ)) g - vecMulVec g (fun _ => (1 : ℝ))) := by
  ext i j
  simp only [Matrix.sub_apply, diagonal_mul, mul_diagonal, Matrix.smul_apply, vecMulVec_apply,
    smul_eq_mul, lowner]
  by_cases h : i = j
  · subst h; simp; ring
  · rw [if_neg h]
    have hne : ν i - ν j ≠ 0 := sub_ne_zero.mpr (fun e => h (hν e))
    have hpi : (π : ℝ) ≠ 0 := Real.pi_ne_zero
    field_simp

end CCMLowner

#print axioms CCMLowner.lowner_symm
#print axioms CCMLowner.lowner_add
#print axioms CCMLowner.lowner_displacement
