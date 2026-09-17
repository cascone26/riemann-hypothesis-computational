import HodgeFormModel

set_option maxHeartbeats 1000000

open IsDedekindDomain NumberField

noncomputable section

variable {K : Type*} [Field K] [NumberField K]

/-!
# §1711 — PRINCIPAL ARAKELOV DIVISORS FORM A HOMOMORPHIC IMAGE inside the degree-zero locus

`principalArakelov : Kˣ → ArakelovDivisor K` is multiplicative-to-additive (`x·y ↦ D_x + D_y`,
`1 ↦ 0`, `x⁻¹ ↦ −D_x`), so its image is a subgroup of the degree-zero locus (§1707:
`arakDeg ∘ principalArakelov = 0`).  This is the arithmetic Picard-group scaffold: the quotient
"degree-zero divisors / principal divisors" is the compactified class group the F₁ story pairs
with the Hodge form (§1710).
-/

theorem principalArakelov_mul (x y : Kˣ) :
    principalArakelov (K := K) (x * y)
      = principalArakelov (K := K) x + principalArakelov (K := K) y := by
  refine ArakelovDivisor.ext K ?_ ?_
  · show principalDiv (x * y) = principalDiv x + principalDiv y
    exact principalDiv_mul x y
  · funext w
    show Real.log (w ((x * y : Kˣ) : K))
      = Real.log (w (x : K)) + Real.log (w (y : K))
    have hx0 : w ((x : K)) ≠ 0 :=
      (InfinitePlace.pos_iff.mpr (Units.ne_zero x)).ne'
    have hy0 : w ((y : K)) ≠ 0 :=
      (InfinitePlace.pos_iff.mpr (Units.ne_zero y)).ne'
    have hmul : ((x * y : Kˣ) : K) = (x : K) * (y : K) := rfl
    rw [hmul, map_mul, Real.log_mul hx0 hy0]

theorem principalArakelov_one :
    principalArakelov (K := K) 1 = 0 := by
  refine ArakelovDivisor.ext K ?_ ?_
  · show principalDiv (1 : Kˣ) = 0
    exact principalDiv_one
  · funext w
    show Real.log (w ((1 : Kˣ) : K)) = 0
    have h1 : ((1 : Kˣ) : K) = 1 := rfl
    rw [h1, map_one, Real.log_one]

theorem principalArakelov_inv (x : Kˣ) :
    principalArakelov (K := K) x⁻¹ = -principalArakelov (K := K) x := by
  have h := principalArakelov_mul (K := K) x x⁻¹
  rw [mul_inv_cancel, principalArakelov_one] at h
  refine ArakelovDivisor.ext K ?_ ?_
  · show principalDiv x⁻¹ = -principalDiv x
    exact principalDiv_inv x
  · funext w
    have hcomp := congrArg (fun D : ArakelovDivisor K => D.infinite w) h
    have hadd : (principalArakelov (K := K) x + principalArakelov (K := K) x⁻¹).infinite w
        = (principalArakelov (K := K) x).infinite w
          + (principalArakelov (K := K) x⁻¹).infinite w := rfl
    have hzero : (0 : ArakelovDivisor K).infinite w = 0 := rfl
    rw [hadd, hzero] at hcomp
    show (principalArakelov (K := K) x⁻¹).infinite w
      = (-(principalArakelov (K := K) x)).infinite w
    have hneg : (-(principalArakelov (K := K) x)).infinite w
        = -((principalArakelov (K := K) x).infinite w) := rfl
    rw [hneg]
    linarith

/-- Every principal Arakelov divisor sits in the degree-zero locus and the Hodge form is
`≤ 0` there (strictly, unless the divisor is trivial): §1707 + §1710 combined. -/
theorem hodgeForm_principalArakelov_nonpos (x : Kˣ) :
    hodgeForm (principalArakelov (K := K) x) (principalArakelov (K := K) x) ≤ 0 := by
  by_cases h : principalArakelov (K := K) x = 0
  · rw [h]
    unfold hodgeForm
    rw [arakDeg_zero, intersectionPairing.zero_left K 0]
    ring_nf
    exact le_refl 0
  · exact le_of_lt (hodgeForm_self_neg_of_deg_zero (arakDeg_principalArakelov x) h)

#print axioms principalArakelov_mul
#print axioms principalArakelov_inv
#print axioms hodgeForm_principalArakelov_nonpos
