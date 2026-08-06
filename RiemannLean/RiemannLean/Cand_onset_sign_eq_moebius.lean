import Mathlib.NumberTheory.ArithmeticFunction.Moebius
import Mathlib.Data.Real.Sign
open ArithmeticFunction
/-- **Sign-theorem logical core** (LEDGER §690–698): for squarefree `n`, if the diagonal onset
factors as `μ(n) · P` with `P > 0` a positive Bessel product (the analytic input, §722/§694),
then its sign is exactly `μ(n)`. The Bessel positivity `P > 0` is proven in `Bessel.lean`
(`besselJ0_prod_pos`); this lemma is the elementary sign algebra `sgn(μ(n)·P) = μ(n)`. -/
theorem onset_sign_eq_moebius {n : ℕ} (hn : Squarefree n) {P : ℝ} (hP : 0 < P) :
    Real.sign ((moebius n : ℝ) * P) = (moebius n : ℝ) := by
  have hne : moebius n ≠ 0 := moebius_ne_zero_iff_squarefree.mpr hn
  rcases moebius_ne_zero_iff_eq_or.mp hne with h1 | h1
  · have : (moebius n : ℝ) = 1 := by rw [h1]; norm_num
    rw [this, one_mul, Real.sign_of_pos hP]
  · have : (moebius n : ℝ) = -1 := by rw [h1]; norm_num
    rw [this, neg_one_mul, Real.sign_of_neg (by linarith)]

#print axioms onset_sign_eq_moebius
