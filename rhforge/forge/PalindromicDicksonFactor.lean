import PowPairDickson
import PalindromicPairSum

/-!
# §1788 — Palindromic factorization P(z) = zᵐ·Q(z+z⁻¹)

The completing capstone of the Chebyshev degree-reduction.  Combining §1787 (pairing sum) with §1786
(pair-to-Dickson), a palindromic coefficient function `v` (`v i = v (2m-i)`) gives, for `z ≠ 0`,

    Σ_{i≤2m} vᵢ zⁱ  =  zᵐ · ( v m  +  Σ_{k<m} v(m+1+k) · D_{k+1}(z+z⁻¹) ).

So a palindromic polynomial of degree `2m` is `zᵐ` times a polynomial `Q` of degree `m` in `s = z+z⁻¹`.
Fed into §1785 (`root_on_circle_iff_chebyshev`), this makes the on-circle criterion apply to the actual
(self-inversive) kernel polynomial: its zeros lie on `|z|=1` ⟺ `Q`'s roots are real in `[-2,2]`.
EXIT0 / std-axioms.
-/

open Polynomial

variable {K : Type*} [Field K]

/-- **Palindromic → `zᵐ·Q(z+z⁻¹)` factorization.**  For `z ≠ 0` and palindromic `v`, the degree-`2m` sum
`Σ vᵢ zⁱ` factors as `zᵐ` times a Dickson combination in `z+z⁻¹`. -/
theorem palindromic_dickson_factor (v : ℕ → K) (z : K) (hz : z ≠ 0) (m : ℕ)
    (hpal : ∀ i, i ≤ 2 * m → v i = v (2 * m - i)) :
    ∑ i ∈ Finset.range (2 * m + 1), v i * z ^ i
      = z ^ m * (v m + ∑ k ∈ Finset.range m, v (m + 1 + k) * (dickson (k + 1)).eval (z + z⁻¹)) := by
  rw [palindromic_pair_sum (fun i => v i * z ^ i) m, mul_add, Finset.mul_sum]
  congr 1
  · ring
  · apply Finset.sum_congr rfl
    intro k hk
    have hkm : k < m := Finset.mem_range.mp hk
    have hpal2 : v (m - 1 - k) = v (m + 1 + k) := by
      rw [hpal (m - 1 - k) (by omega), show 2 * m - (m - 1 - k) = m + 1 + k by omega]
    have h86 : z ^ (m + 1 + k) + z ^ (m - 1 - k)
        = z ^ m * (dickson (k + 1)).eval (z + z⁻¹) := by
      rw [show m + 1 + k = m + (k + 1) by ring, show m - 1 - k = m - (k + 1) by omega]
      exact pow_pair_eq_dickson z hz m (k + 1) (by omega)
    rw [hpal2, ← mul_add, h86]
    ring
