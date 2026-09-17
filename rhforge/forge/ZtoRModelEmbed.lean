import ArakelovFullDegree
import ArakelovVecRHodge

set_option maxHeartbeats 1000000

open IsDedekindDomain NumberField

noncomputable section

variable {K : Type*} [Field K] [NumberField K]

/-!
# §1717 — THE ℤ-MODEL EMBEDS IN THE ℝ-MODEL, pairing- and degree-compatibly

`toVecR : ArakelovDivisor K → ArakelovVecR K` casts the finite ℤ-coefficients to ℝ and keeps the
archimedean part.  It is additive, matches the intersection pairing (`pairingFun ∘ toVecR = ⟨,⟩`)
and the full degree (`degFun ∘ toVecR = arakDeg`) — so the whole §1707–§1711 ℤ-story (principal
divisors, degree-zero subgroup, Hodge negativity) transports into the ℝ-model where §1716's
signature theorem lives: in particular every principal Arakelov divisor lands in `ker degLin`.
-/

/-- Cast the finite part ℤ → ℝ, keep the archimedean part. -/
noncomputable def toVecR (D : ArakelovDivisor K) : ArakelovVecR K :=
  (D.finite.mapRange (fun n : ℤ => (n : ℝ)) (by norm_num), D.infinite)

theorem toVecR_add (D E : ArakelovDivisor K) :
    toVecR (D + E) = toVecR D + toVecR E := by
  unfold toVecR
  have hfin : (D + E).finite = D.finite + E.finite := rfl
  have hinf : (D + E).infinite = fun w => D.infinite w + E.infinite w := rfl
  rw [hfin, hinf]
  refine Prod.ext ?_ rfl
  show (D.finite + E.finite).mapRange (fun n : ℤ => (n : ℝ)) (by norm_num)
      = D.finite.mapRange (fun n : ℤ => (n : ℝ)) (by norm_num)
        + E.finite.mapRange (fun n : ℤ => (n : ℝ)) (by norm_num)
  ext v
  simp [Finsupp.mapRange_apply]

theorem pairingFun_toVecR (D E : ArakelovDivisor K) :
    ArakelovVecR.pairingFun K (toVecR D) (toVecR E) = intersectionPairing K D E := by
  unfold ArakelovVecR.pairingFun intersectionPairing toVecR
  congr 1
  rw [Finsupp.sum_mapRange_index (fun v => zero_mul _)]
  apply Finsupp.sum_congr
  intro v _
  rw [Finsupp.mapRange_apply]

theorem degFun_toVecR (D : ArakelovDivisor K) :
    ArakelovVecR.degFun K (toVecR D) = arakDeg D := by
  unfold ArakelovVecR.degFun arakDeg toVecR deg
  congr 1
  · rw [Finsupp.sum_mapRange_index (fun v => zero_mul _)]
  · rw [finsum_eq_sum_of_fintype]
    apply Finset.sum_congr rfl
    intro w _
    rw [nsmul_eq_mul]

/-- Principal Arakelov divisors land in the kernel of the ℝ-model degree. -/
theorem degLin_toVecR_principal (x : Kˣ) :
    ArakelovVecR.degLin K (toVecR (principalArakelov (K := K) x)) = 0 := by
  rw [ArakelovVecR.degLin_apply, degFun_toVecR, arakDeg_principalArakelov]

#print axioms pairingFun_toVecR
#print axioms degLin_toVecR_principal
