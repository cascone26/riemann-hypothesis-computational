import RiemannLean.Bessel
open scoped Nat
/-- `J₀` is even: `J₀(-x) = J₀(x)` (each series term depends on `x` only through `(x/2)^{2k}`). -/
theorem besselJ0_even (x : ℝ) : besselJ0 (-x) = besselJ0 x := by
  unfold besselJ0
  congr 1; ext k
  rw [show (-x) / 2 = -(x / 2) by ring, (even_two_mul k).neg_pow]

#print axioms besselJ0_even
