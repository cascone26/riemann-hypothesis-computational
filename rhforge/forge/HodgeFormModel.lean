import ArakelovFullDegree

set_option maxHeartbeats 1000000

open IsDedekindDomain NumberField

noncomputable section

variable {K : Type*} [Field K] [NumberField K]

/-!
# §1710 — THE HODGE FORM ON THE CONCRETE MODEL: negdef on the degree-zero locus

The intersection pairing (§1652-family) is positive definite (it is an ℓ²-type form), so the
Hodge-shaped form on `ArakelovDivisor` is `hodgeForm D E := arakDeg D · arakDeg E − ⟨D,E⟩`:
symmetric, bilinear, and **strictly negative on nonzero degree-zero divisors** — exactly the
negative-definite-on-the-primitive-part hypothesis that the coordinate-free arithmetic
Hodge-index theorem (§1659 `hodge_signature_of_adjunction`) consumes.  This brick proves the
pairing's definiteness and the model's Hodge-shape.
-/

theorem intersectionPairing_self_nonneg (D : ArakelovDivisor K) :
    0 ≤ intersectionPairing K D D := by
  unfold intersectionPairing
  apply add_nonneg
  · rw [Finsupp.sum]
    apply Finset.sum_nonneg
    intro v _
    exact mul_self_nonneg _
  · apply Finset.sum_nonneg
    intro w _
    have h1 : (0:ℝ) ≤ D.infinite w * D.infinite w := mul_self_nonneg _
    have h2 : (0:ℝ) ≤ (w.mult : ℝ) := Nat.cast_nonneg _
    exact mul_nonneg h1 h2

theorem intersectionPairing_self_eq_zero_iff (D : ArakelovDivisor K) :
    intersectionPairing K D D = 0 ↔ D = 0 := by
  constructor
  · intro h
    have hfin0 : (0:ℝ) ≤ D.finite.sum (fun v n => (n : ℝ) * (D.finite v : ℝ)) := by
      rw [Finsupp.sum]
      exact Finset.sum_nonneg fun v _ => mul_self_nonneg _
    have hinf0 : (0:ℝ) ≤ ∑ w : InfinitePlace K,
        D.infinite w * D.infinite w * (w.mult : ℝ) :=
      Finset.sum_nonneg fun w _ =>
        mul_nonneg (mul_self_nonneg _) (Nat.cast_nonneg _)
    have hsum : D.finite.sum (fun v n => (n : ℝ) * (D.finite v : ℝ))
        + ∑ w : InfinitePlace K, D.infinite w * D.infinite w * (w.mult : ℝ) = 0 := h
    have hf : D.finite.sum (fun v n => (n : ℝ) * (D.finite v : ℝ)) = 0 := by linarith
    have hi : (∑ w : InfinitePlace K, D.infinite w * D.infinite w * (w.mult : ℝ)) = 0 := by
      linarith
    refine ArakelovDivisor.ext K ?_ ?_
    · -- finite part vanishes
      rw [Finsupp.sum] at hf
      have hall := (Finset.sum_eq_zero_iff_of_nonneg
        (fun v _ => mul_self_nonneg ((D.finite v : ℝ)))).mp hf
      have hsupp : D.finite.support = ∅ := by
        apply Finset.eq_empty_of_forall_notMem
        intro v hv
        have h1 := hall v hv
        have h2 : (D.finite v : ℝ) = 0 := by
          have := mul_self_eq_zero.mp h1
          exact this
        have h3 : D.finite v = 0 := by exact_mod_cast h2
        exact (Finsupp.mem_support_iff.mp hv) h3
      exact Finsupp.support_eq_empty.mp hsupp
    · -- infinite part vanishes
      have hall := (Finset.sum_eq_zero_iff_of_nonneg
        (fun w _ => mul_nonneg (mul_self_nonneg (D.infinite w)) (Nat.cast_nonneg w.mult))).mp hi
      funext w
      have h1 := hall w (Finset.mem_univ w)
      have hm : (0:ℝ) < (w.mult : ℝ) := by exact_mod_cast InfinitePlace.mult_pos
      have h2 : D.infinite w * D.infinite w = 0 := by
        by_contra hne
        have := mul_ne_zero hne hm.ne'
        exact this h1
      exact mul_self_eq_zero.mp h2
  · intro h
    rw [h]
    exact intersectionPairing.zero_left K 0

theorem intersectionPairing_self_pos_of_ne_zero {D : ArakelovDivisor K} (hD : D ≠ 0) :
    0 < intersectionPairing K D D :=
  lt_of_le_of_ne (intersectionPairing_self_nonneg D)
    (fun h => hD ((intersectionPairing_self_eq_zero_iff D).mp h.symm))

/-- The Hodge-shaped form on the concrete Arakelov model. -/
noncomputable def hodgeForm (D E : ArakelovDivisor K) : ℝ :=
  arakDeg D * arakDeg E - intersectionPairing K D E

theorem hodgeForm_symm (D E : ArakelovDivisor K) :
    hodgeForm D E = hodgeForm E D := by
  unfold hodgeForm
  rw [intersectionPairing.symmetric]
  ring

theorem hodgeForm_add_left (D D' E : ArakelovDivisor K) :
    hodgeForm (D + D') E = hodgeForm D E + hodgeForm D' E := by
  unfold hodgeForm
  rw [arakDeg_add, intersectionPairing.bilinear_add_left]
  ring

theorem hodgeForm_add_right (D E E' : ArakelovDivisor K) :
    hodgeForm D (E + E') = hodgeForm D E + hodgeForm D E' := by
  unfold hodgeForm
  rw [arakDeg_add, intersectionPairing.bilinear_add_right]
  ring

/-- **The Hodge-index shape on the model**: the Hodge form is strictly negative on every nonzero
degree-zero divisor — the negative-definite-primitive-part hypothesis of §1659, realized. -/
theorem hodgeForm_self_neg_of_deg_zero {D : ArakelovDivisor K}
    (hdeg : arakDeg D = 0) (hD : D ≠ 0) :
    hodgeForm D D < 0 := by
  unfold hodgeForm
  rw [hdeg]
  have := intersectionPairing_self_pos_of_ne_zero hD
  linarith

#print axioms intersectionPairing_self_eq_zero_iff
#print axioms hodgeForm_self_neg_of_deg_zero
