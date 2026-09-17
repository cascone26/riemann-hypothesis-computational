import IntersectionPairing

set_option maxHeartbeats 1000000

open IsDedekindDomain NumberField
noncomputable section

variable (K : Type*) [Field K] [NumberField K]

namespace intersectionPairing

/-- §1652 — the intersection pairing vanishes when the first argument is zero. -/
theorem zero_left (E : ArakelovDivisor K) : intersectionPairing K 0 E = 0 := by
  unfold intersectionPairing
  show (0 : HeightOneSpectrum (RingOfIntegers K) →₀ ℤ).sum (fun v n => (n : ℝ) * (E.finite v : ℝ)) +
      ∑ w : InfinitePlace K, (0 : ℝ) * E.infinite w * (w.mult : ℝ) = 0
  rw [Finsupp.sum_zero_index]
  simp

/-- §1652 — the intersection pairing vanishes when the second argument is zero (via `symmetric`). -/
theorem zero_right (D : ArakelovDivisor K) : intersectionPairing K D 0 = 0 := by
  rw [symmetric, zero_left]

/-- The local `ArakelovDivisor K` here has manual `Add`/`Neg`/`Sub`/`Zero` instances but no bundled
`AddGroup` (this is a self-contained chain, independent of the `Arakelov.lean`/`embFinite` chain which
does carry the bundled group structure) — so generic `AddGroup` lemmas like `add_neg_cancel` don't
apply directly to `ArakelovDivisor K` terms. Cancel componentwise via `ext` instead, reducing to the
`Finsupp`/`Pi` level where the ambient `AddCommGroup` instances DO apply. -/
private theorem add_neg_self (D : ArakelovDivisor K) : D + -D = (0 : ArakelovDivisor K) := by
  apply ArakelovDivisor.ext
  · show D.finite + -D.finite = 0
    exact add_neg_cancel D.finite
  · funext w
    show D.infinite w + -D.infinite w = 0
    exact add_neg_cancel (D.infinite w)

private theorem sub_eq_add_neg' (D D' : ArakelovDivisor K) : D - D' = D + -D' := by
  apply ArakelovDivisor.ext
  · show D.finite - D'.finite = D.finite + -D'.finite
    exact sub_eq_add_neg D.finite D'.finite
  · funext w
    show D.infinite w - D'.infinite w = D.infinite w + -D'.infinite w
    exact sub_eq_add_neg (D.infinite w) (D'.infinite w)

/-- §1652 — the pairing is additive-inverse-compatible in the first argument:
`⟨-D, E⟩ = -⟨D, E⟩`. From bilinearity: `⟨D,E⟩ + ⟨-D,E⟩ = ⟨D + (-D), E⟩ = ⟨0,E⟩ = 0`. -/
theorem neg_left (D E : ArakelovDivisor K) :
    intersectionPairing K (-D) E = -intersectionPairing K D E := by
  have h := bilinear_add_left K D E (-D)
  rw [add_neg_self, zero_left] at h
  linarith

/-- §1652 — companion: `⟨D, -E⟩ = -⟨D, E⟩` (via symmetry + `neg_left`). -/
theorem neg_right (D E : ArakelovDivisor K) :
    intersectionPairing K D (-E) = -intersectionPairing K D E := by
  rw [symmetric, neg_left, symmetric K E D]

/-- §1652 — subtraction in the first argument distributes over the pairing. -/
theorem bilinear_sub_left (D E D' : ArakelovDivisor K) :
    intersectionPairing K (D - D') E =
    intersectionPairing K D E - intersectionPairing K D' E := by
  rw [sub_eq_add_neg', bilinear_add_left K D E (-D'), neg_left]
  ring

/-- §1652 — subtraction in the second argument distributes over the pairing. -/
theorem bilinear_sub_right (D E E' : ArakelovDivisor K) :
    intersectionPairing K D (E - E') =
    intersectionPairing K D E - intersectionPairing K D E' := by
  rw [sub_eq_add_neg', bilinear_add_right K D E (-E'), neg_right]
  ring

end intersectionPairing

#print axioms intersectionPairing.zero_left
#print axioms intersectionPairing.zero_right
#print axioms intersectionPairing.neg_left
#print axioms intersectionPairing.neg_right
#print axioms intersectionPairing.bilinear_sub_left
#print axioms intersectionPairing.bilinear_sub_right
