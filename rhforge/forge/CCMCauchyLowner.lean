import CCMArithLowner

/-!
# S2972 — the zero side is Löwner too: Cauchy Gram G = C·diag(w)·Cᵀ has generator π·(C w) (lane 2)

The CCM zero-side vectors are Cauchy: ⟨e_n, ρ_γ⟩ ∝ (e^{−iγL} − 1)/(i(2πn/L − γ)), i.e. after moving
|e^{−iγL} − 1|² = 4 sin²(γL/2) into the weight w, `C_{nγ} = 1/(ν_n − x_γ)`.  For ANY Cauchy-type C
(hypothesis `(ν_i − x_j)·C_ij = 1`) and any weight w:
* `cauchy_gram_comm_entry` — `(ν_i − ν_k)·G_ik = (Cw)_k − (Cw)_i`.
* `cauchy_gram_displacement` — `[diag ν, G] = 𝟙(Cw)ᵀ − (Cw)𝟙ᵀ`  (rank ≤ 2).
* `cauchy_gram_eq_lowner` — injective ν ⇒ `G = lowner (diag G) (π•(Cw)) ν`.
* `ccm_full_lowner` — with S2970: `A − G = lowner (d_A − diag G) (g_A − π•(Cw)) ν`.
Reading (Report): the full CCM comparison matrix A − G (arithmetic minus zero side) is ONE Löwner matrix whose
generator is `g_A − π·Cw`: the explicit-formula balance is a statement about this single function of the frequency.
-/

open Matrix Finset Real

namespace CCMCauchy

variable {n m : Type*} [Fintype m] [DecidableEq n] [DecidableEq m]

omit [DecidableEq n] in
theorem cauchy_gram_comm_entry (ν : n → ℝ) (x : m → ℝ) (C : Matrix n m ℝ) (w : m → ℝ)
    (hC : ∀ i j, (ν i - x j) * C i j = 1) (i k : n) :
    (ν i - ν k) * (C * diagonal w * Cᵀ) i k = (C *ᵥ w) k - (C *ᵥ w) i := by
  rw [Matrix.mul_apply]
  simp only [Matrix.mul_diagonal, transpose_apply, mulVec, dotProduct]
  rw [Finset.mul_sum, ← Finset.sum_sub_distrib]
  apply Finset.sum_congr rfl
  intro j _
  linear_combination (w j * C k j) * hC i j - (w j * C i j) * hC k j

theorem cauchy_gram_displacement [Fintype n] (ν : n → ℝ) (x : m → ℝ) (C : Matrix n m ℝ)
    (w : m → ℝ) (hC : ∀ i j, (ν i - x j) * C i j = 1) :
    diagonal ν * (C * diagonal w * Cᵀ) - (C * diagonal w * Cᵀ) * diagonal ν
      = vecMulVec (fun _ => (1 : ℝ)) (C *ᵥ w) - vecMulVec (C *ᵥ w) (fun _ => (1 : ℝ)) := by
  ext i k
  simp only [Matrix.sub_apply, diagonal_mul, mul_diagonal, vecMulVec_apply, one_mul, mul_one]
  have h := cauchy_gram_comm_entry ν x C w hC i k
  linear_combination h

theorem cauchy_gram_eq_lowner (ν : n → ℝ) (x : m → ℝ) (C : Matrix n m ℝ) (w : m → ℝ)
    (hC : ∀ i j, (ν i - x j) * C i j = 1) (hν : Function.Injective ν) :
    C * diagonal w * Cᵀ
      = CCMLowner.lowner (fun i => (C * diagonal w * Cᵀ) i i) (π • (C *ᵥ w)) ν := by
  ext i k
  simp only [CCMLowner.lowner]
  by_cases h : i = k
  · subst h; simp
  · rw [if_neg h]
    have hne : ν i - ν k ≠ 0 := sub_ne_zero.mpr (fun e => h (hν e))
    have hpi : (π : ℝ) ≠ 0 := Real.pi_ne_zero
    have he := cauchy_gram_comm_entry ν x C w hC i k
    simp only [Pi.smul_apply, smul_eq_mul]
    rw [eq_div_iff (mul_ne_zero hpi hne)]
    linear_combination π * he

theorem ccm_full_lowner {L : ℝ} (hL : 0 < L) (ν : n → ℝ) (hν : Function.Injective ν)
    (x : m → ℝ) (C : Matrix n m ℝ) (w : m → ℝ) (hC : ∀ i j, (ν i - x j) * C i j = 1)
    (dR gR dP gP : n → ℝ) :
    (CCMW02.W02 L ν - CCMLowner.lowner dR gR ν - CCMLowner.lowner dP gP ν) - C * diagonal w * Cᵀ
      = CCMLowner.lowner
          ((fun i => CCMW02.W02 L ν i i) - dR - dP - fun i => (C * diagonal w * Cᵀ) i i)
          (-(CCMW02.cW L / 4) • CCMW02.vW L ν - gR - gP - π • (C *ᵥ w)) ν := by
  rw [CCMArith.ccm_arith_lowner hL hν]
  conv_lhs => rw [cauchy_gram_eq_lowner ν x C w hC hν]
  rw [CCMArith.lowner_sub]

end CCMCauchy

#print axioms CCMCauchy.cauchy_gram_comm_entry
#print axioms CCMCauchy.cauchy_gram_displacement
#print axioms CCMCauchy.cauchy_gram_eq_lowner
#print axioms CCMCauchy.ccm_full_lowner
