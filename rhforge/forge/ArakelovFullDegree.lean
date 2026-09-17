import PrincipalArakelovDegree
import IntersectionPairingZeroNeg

set_option maxHeartbeats 1000000

open IsDedekindDomain NumberField

noncomputable section

variable {K : Type*} [Field K] [NumberField K]

/-!
# §1707 — THE FULL ARAKELOV DEGREE on `ArakelovDivisor`: additive, principal ⇒ 0

Unites the two F₁ scaffolds: the degree theory (`deg`/`principalDiv`/product formula, §1003-family)
lives on bare finite divisors, while the intersection pairing (§1652-family) lives on
`ArakelovDivisor` (finite ⊕ archimedean).  This brick puts the FULL degree
`arakDeg D = deg D.finite + Σᶠ_w mult(w)·D.infinite(w)` on the pairing's type, proves it additive
(so the degree-zero locus is a subgroup of the pairing's world — where the arithmetic
Hodge-index-shape statements of §1659 live), and shows principal Arakelov divisors have degree
zero — the arithmetic analog of "principal divisors have degree 0 on a curve", direct from the
banked product formula (§1003 → arDeg_principalArakelov_eq_zero).
-/

/-- The full Arakelov degree: finite arithmetic degree + multiplicity-weighted archimedean part. -/
noncomputable def arakDeg (D : ArakelovDivisor K) : ℝ :=
  deg D.finite + ∑ᶠ w : InfinitePlace K, w.mult • D.infinite w

/-- The principal Arakelov divisor of a field unit: valuation divisor + archimedean logs. -/
noncomputable def principalArakelov (x : Kˣ) : ArakelovDivisor K :=
  ⟨principalDiv x, fun w => Real.log (w (x : K))⟩

theorem arakDeg_zero : arakDeg (0 : ArakelovDivisor K) = 0 := by
  have hfin : (0 : ArakelovDivisor K).finite = 0 := rfl
  have hinf : ∀ w : InfinitePlace K, (0 : ArakelovDivisor K).infinite w = 0 := fun _ => rfl
  unfold arakDeg
  rw [hfin]
  simp only [hinf, smul_zero, finsum_zero]
  unfold deg
  rw [Finsupp.sum_zero_index]
  ring

theorem arakDeg_add (D E : ArakelovDivisor K) :
    arakDeg (D + E) = arakDeg D + arakDeg E := by
  have hfin : (D + E).finite = D.finite + E.finite := rfl
  have hinf : ∀ w : InfinitePlace K, (D + E).infinite w = D.infinite w + E.infinite w :=
    fun _ => rfl
  unfold arakDeg
  rw [hfin]
  have h1 : deg (D.finite + E.finite) = deg D.finite + deg E.finite := by
    unfold deg
    rw [Finsupp.sum_add_index']
    · intro a
      simp
    · intro a b₁ b₂
      push_cast
      ring
  have h2 : (∑ᶠ w : InfinitePlace K, w.mult • (D + E).infinite w)
      = (∑ᶠ w : InfinitePlace K, w.mult • D.infinite w)
        + ∑ᶠ w : InfinitePlace K, w.mult • E.infinite w := by
    calc (∑ᶠ w : InfinitePlace K, w.mult • (D + E).infinite w)
        = ∑ᶠ w : InfinitePlace K,
            (w.mult • D.infinite w + w.mult • E.infinite w) := by
          simp only [hinf, smul_add]
      _ = _ := finsum_add_distrib (Set.toFinite _) (Set.toFinite _)
  rw [h1, h2]
  ring

theorem arakDeg_neg (D : ArakelovDivisor K) : arakDeg (-D) = -arakDeg D := by
  have hD : D + -D = (0 : ArakelovDivisor K) := by
    refine ArakelovDivisor.ext K ?_ ?_
    · show D.finite + -D.finite = 0
      simp
    · funext w
      show D.infinite w + -D.infinite w = 0
      ring
  have h := arakDeg_add D (-D)
  rw [hD, arakDeg_zero] at h
  linarith

theorem arakDeg_sub (D E : ArakelovDivisor K) :
    arakDeg (D - E) = arakDeg D - arakDeg E := by
  have h : D - E = D + -E := by
    refine ArakelovDivisor.ext K ?_ ?_
    · show D.finite - E.finite = D.finite + -E.finite
      exact sub_eq_add_neg _ _
    · funext w
      show D.infinite w - E.infinite w = D.infinite w + -E.infinite w
      ring
  rw [h, arakDeg_add, arakDeg_neg]
  ring

/-- **Principal Arakelov divisors have full degree zero** — the arithmetic product formula in the
intersection-pairing's world. -/
theorem arakDeg_principalArakelov (x : Kˣ) :
    arakDeg (principalArakelov (K := K) x) = 0 := by
  unfold arakDeg principalArakelov
  exact arDeg_principalArakelov_eq_zero x

#print axioms arakDeg_add
#print axioms arakDeg_principalArakelov
