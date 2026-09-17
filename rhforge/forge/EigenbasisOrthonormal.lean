import ShiftOverlapDiagonal
import ShiftOverlapOffDiag
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic

/-!
# §1734 — The zero-shift limit recovers the `ℓ²` toy: the eigenbasis is orthonormal

At Arakelov distance `s = 0` the shift-overlap kernel `S_kl(s)` (§1730–§1733) reduces to the plain Gram matrix
of the interval eigenbasis `vₖ(u) = sin(ωₖ u)`, `ωₖ·L = k·π`:

* **orthogonality** `∫₀ᴸ sin(ωₖ u)·sin(ωₗ u) du = 0`  for `k ≠ l` (`a ≠ ±b`), and
* **normalization** `∫₀ᴸ sin(ωₖ u)² du = L/2`.

So the `δ_kl` (`ℓ²`) toy that the Lean Arakelov pairing used is EXACTLY the `s→0` limit of the true finite
intersection kernel — the diagonal `L/2` is the divisor self-energy at zero shift, and the off-diagonal
vanishing is the eigenbasis orthogonality.  This pins precisely where the toy is correct (`s=0`) and where the
genuine geometry lives (`s = log pᵏ > 0`, §1730–§1733).  EXIT0 / std-axioms.
-/

open Real intervalIntegral

/-- **Eigenbasis orthogonality** (zero-shift off-diagonal).  For `a·L = k·π`, `b·L = l·π`, `a ≠ ±b`. -/
theorem eigenbasis_orthogonal (a b L : ℝ) (k l : ℤ)
    (hd : a - b ≠ 0) (hs : a + b ≠ 0)
    (hLa : a * L = (k : ℝ) * π) (hLb : b * L = (l : ℝ) * π) :
    (∫ u in (0:ℝ)..L, Real.sin (a * u) * Real.sin (b * u)) = 0 := by
  have h := shift_overlap_offdiag a b 0 L hd hs
  simp only [add_zero, mul_zero, Real.sin_zero, sub_zero] at h
  rw [h]
  have hm : Real.sin ((a - b) * L) = 0 := by
    rw [show (a - b) * L = ((k - l : ℤ) : ℝ) * π by
      rw [show (a - b) * L = a * L - b * L by ring, hLa, hLb]; push_cast; ring]
    exact Real.sin_int_mul_pi _
  have hp : Real.sin ((a + b) * L) = 0 := by
    rw [show (a + b) * L = ((k + l : ℤ) : ℝ) * π by
      rw [show (a + b) * L = a * L + b * L by ring, hLa, hLb]; push_cast; ring]
    exact Real.sin_int_mul_pi _
  rw [hm, hp]; simp

/-- **Eigenmode normalization** (zero-shift diagonal).  For `a·L = k·π`, `a ≠ 0`. -/
theorem eigenmode_norm (a L : ℝ) (k : ℤ) (ha : a ≠ 0) (hLa : a * L = (k : ℝ) * π) :
    (∫ u in (0:ℝ)..L, Real.sin (a * u) * Real.sin (a * u)) = L / 2 := by
  have h := shift_overlap_diagonal a 0 L k ha hLa
  simp only [add_zero, mul_zero, sub_zero, Real.cos_zero, Real.sin_zero, mul_one, zero_div] at h
  rw [h]
