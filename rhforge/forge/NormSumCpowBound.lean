import SumCpowNormBound
import SumRpowNegLe

/-!
# §1290 — Modulus bound on ζ's finite main sum: `‖∑_{k≤N} k^{-s}‖ ≤ (N^{1-σ}-1)/(1-σ)`, `σ = Re s`

Chaining the triangle/modulus bound §1287 (`‖∑ k^{-s}‖ ≤ ∑ k^{-σ}`) with the integral-comparison growth
bound §1289 (`∑_{k=2}^N k^{-σ} ≤ (N^{1-σ}-1)/(1-σ)`).  For `σ = Re s` in `[0,1)` the RHS is `≍ N^{1-σ}/(1-σ)`;
with `N ≍ |t|` this is the size control on the truncated Euler–Maclaurin main sum of §1286 → §1269's conditional
log-growth bound → sharp zero-free region.  EXIT0 / std-axioms.
-/

open Finset

theorem norm_sum_cpow_neg_le_bound {N : ℕ} (hN : 1 ≤ N) {s : ℂ}
    (hs0 : 0 ≤ s.re) (hs1 : s.re ≠ 1) :
    ‖∑ k ∈ Finset.Ioc 1 N, (k : ℂ) ^ (-s)‖ ≤ ((N : ℝ) ^ (1 - s.re) - 1) / (1 - s.re) :=
  le_trans norm_sum_cpow_neg_le (sum_rpow_neg_le hN hs0 hs1)

#print axioms norm_sum_cpow_neg_le_bound
