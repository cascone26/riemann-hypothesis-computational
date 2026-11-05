/-
PartialFractionSignChange.lean  (HP RH center, LEDGER S2861, 2026-11-04)

Elementary parity lemma behind S2860's zero-counting sign-change law for the CCM/Shi ground state.
For R(y) = ∑ₖ uₖ/(y - xₖ) and two poles xᵢ < xⱼ with no other pole in [xᵢ, xⱼ]:
if uᵢ, uⱼ have the SAME sign then R has a zero strictly between xᵢ and xⱼ.
Contrapositive (`mul_neg_of_no_zero_between`): no zero between consecutive poles ⇒ the residues
change sign.  Summed over node intervals this gives
  #sign changes of u = (n-1) - #{node intervals containing an odd number of zeros of R}.
-/
import Mathlib

open Filter Topology Set

set_option linter.unusedSectionVars false

namespace PFSignChange

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- the partial-fraction function -/
noncomputable def R (x u : ι → ℝ) (y : ℝ) : ℝ := ∑ k, u k / (y - x k)

lemma R_split (x u : ι → ℝ) (i : ι) (y : ℝ) :
    R x u y = u i * (y - x i)⁻¹ + ∑ k ∈ Finset.univ.erase i, u k / (y - x k) := by
  unfold R
  rw [← Finset.add_sum_erase _ _ (Finset.mem_univ i), div_eq_mul_inv]

lemma R_neg (x u : ι → ℝ) (y : ℝ) : R x (fun k => - u k) y = - R x u y := by
  unfold R; rw [← Finset.sum_neg_distrib]; congr 1; ext k; ring

section
variable (x u : ι → ℝ) (i j : ι) (hij : x i < x j)
    (hout : ∀ k, k ≠ i → k ≠ j → x k ∉ Icc (x i) (x j))
include hij hout

lemma ne_of_mem_Ioo {y : ℝ} (hy : y ∈ Ioo (x i) (x j)) (k : ι) : y - x k ≠ 0 := by
  intro h
  have hyk : y = x k := by linarith
  by_cases hki : k = i
  · subst hki; linarith [hy.1]
  by_cases hkj : k = j
  · subst hkj; linarith [hy.2]
  exact hout k hki hkj ⟨by linarith [hy.1], by linarith [hy.2]⟩

lemma R_continuousOn : ContinuousOn (R x u) (Ioo (x i) (x j)) := by
  unfold R
  apply continuousOn_finsetSum
  intro k _
  exact continuousOn_const.div (continuousOn_id.sub continuousOn_const)
    (fun y hy => ne_of_mem_Ioo x i j hij hout hy k)

lemma rest_continuousAt_left :
    ContinuousAt (fun y => ∑ k ∈ Finset.univ.erase i, u k / (y - x k)) (x i) := by
  unfold ContinuousAt
  refine tendsto_finsetSum _ fun k hk => ?_
  have hki : k ≠ i := Finset.ne_of_mem_erase hk
  have hne : x i - x k ≠ 0 := by
    intro h
    have : x i = x k := by linarith
    by_cases hkj : k = j
    · subst hkj; linarith
    exact hout k hki hkj ⟨by linarith, by linarith⟩
  exact (continuousAt_const.div (continuousAt_id.sub continuousAt_const) hne).tendsto

lemma rest_continuousAt_right :
    ContinuousAt (fun y => ∑ k ∈ Finset.univ.erase j, u k / (y - x k)) (x j) := by
  unfold ContinuousAt
  refine tendsto_finsetSum _ fun k hk => ?_
  have hkj : k ≠ j := Finset.ne_of_mem_erase hk
  have hne : x j - x k ≠ 0 := by
    intro h
    have : x j = x k := by linarith
    by_cases hki : k = i
    · subst hki; linarith
    exact hout k hki hkj ⟨by linarith, by linarith⟩
  exact (continuousAt_const.div (continuousAt_id.sub continuousAt_const) hne).tendsto

lemma R_tendsto_left (hui : 0 < u i) : Tendsto (R x u) (𝓝[>] (x i)) atTop := by
  have h1 : Tendsto (fun y : ℝ => y - x i) (𝓝[>] (x i)) (𝓝[>] 0) := by
    apply tendsto_nhdsWithin_of_tendsto_nhds_of_eventually_within
    · have : Tendsto (fun y : ℝ => y - x i) (𝓝 (x i)) (𝓝 (x i - x i)) :=
        (continuous_id.sub continuous_const).tendsto (x i)
      simpa using this.mono_left nhdsWithin_le_nhds
    · filter_upwards [self_mem_nhdsWithin] with y hy
      simpa using hy
  have h2 : Tendsto (fun y : ℝ => u i * (y - x i)⁻¹) (𝓝[>] (x i)) atTop :=
    (tendsto_inv_nhdsGT_zero.comp h1).const_mul_atTop hui
  have h3 := (rest_continuousAt_left x u i j hij hout).tendsto.mono_left
    (nhdsWithin_le_nhds (s := Ioi (x i)))
  have := h2.atTop_add h3
  refine this.congr' (Eventually.of_forall fun y => ?_)
  exact (R_split x u i y).symm

lemma R_tendsto_right (huj : 0 < u j) : Tendsto (R x u) (𝓝[<] (x j)) atBot := by
  have h1 : Tendsto (fun y : ℝ => y - x j) (𝓝[<] (x j)) (𝓝[<] 0) := by
    apply tendsto_nhdsWithin_of_tendsto_nhds_of_eventually_within
    · have : Tendsto (fun y : ℝ => y - x j) (𝓝 (x j)) (𝓝 (x j - x j)) :=
        (continuous_id.sub continuous_const).tendsto (x j)
      simpa using this.mono_left nhdsWithin_le_nhds
    · filter_upwards [self_mem_nhdsWithin] with y hy
      simpa using hy
  have h2 : Tendsto (fun y : ℝ => u j * (y - x j)⁻¹) (𝓝[<] (x j)) atBot :=
    (tendsto_inv_nhdsLT_zero.comp h1).const_mul_atBot huj
  have h3 := (rest_continuousAt_right x u i j hij hout).tendsto.mono_left
    (nhdsWithin_le_nhds (s := Iio (x j)))
  have := h2.atBot_add h3
  refine this.congr' (Eventually.of_forall fun y => ?_)
  exact (R_split x u j y).symm

/-- **Main lemma.** Same-sign (positive) residues at consecutive poles force a zero between them. -/
theorem exists_zero_between_of_pos (hui : 0 < u i) (huj : 0 < u j) :
    ∃ y ∈ Ioo (x i) (x j), R x u y = 0 := by
  have hL : ∀ᶠ y in 𝓝[>] (x i), 0 < R x u y ∧ y ∈ Ioo (x i) (x j) :=
    ((R_tendsto_left x u i j hij hout hui).eventually_gt_atTop 0).and (Ioo_mem_nhdsGT hij)
  obtain ⟨a, ha, haI⟩ := hL.exists
  have hR : ∀ᶠ y in 𝓝[<] (x j), R x u y < 0 ∧ y ∈ Ioo a (x j) :=
    ((R_tendsto_right x u i j hij hout huj).eventually_lt_atBot 0).and (Ioo_mem_nhdsLT haI.2)
  obtain ⟨b, hb, hbI⟩ := hR.exists
  have hab : a ≤ b := le_of_lt hbI.1
  have hsub : Icc a b ⊆ Ioo (x i) (x j) := fun y hy =>
    ⟨lt_of_lt_of_le haI.1 hy.1, lt_of_le_of_lt hy.2 hbI.2⟩
  have hc : ContinuousOn (R x u) (Icc a b) := (R_continuousOn x u i j hij hout).mono hsub
  obtain ⟨c, hcI, hc0⟩ := intermediate_value_Icc' hab hc ⟨le_of_lt hb, le_of_lt ha⟩
  exact ⟨c, hsub hcI, hc0⟩

/-- Same sign (product positive) ⇒ a zero between the poles. -/
theorem exists_zero_between_of_mul_pos (h : 0 < u i * u j) :
    ∃ y ∈ Ioo (x i) (x j), R x u y = 0 := by
  rcases mul_pos_iff.mp h with ⟨hi, hj⟩ | ⟨hi, hj⟩
  · exact exists_zero_between_of_pos x u i j hij hout hi hj
  · obtain ⟨y, hy, h0⟩ := exists_zero_between_of_pos x (fun k => - u k) i j hij hout
      (by linarith) (by linarith)
    exact ⟨y, hy, by rw [R_neg] at h0; linarith⟩

/-- **Sign-change law (contrapositive).** No zero of R between consecutive poles, nonzero residues
⇒ the residues have opposite signs. -/
theorem mul_neg_of_no_zero_between (hi : u i ≠ 0) (hj : u j ≠ 0)
    (hno : ∀ y ∈ Ioo (x i) (x j), R x u y ≠ 0) : u i * u j < 0 := by
  rcases lt_trichotomy (u i * u j) 0 with h | h | h
  · exact h
  · exact absurd h (mul_ne_zero hi hj)
  · obtain ⟨y, hy, h0⟩ := exists_zero_between_of_mul_pos x u i j hij hout h
    exact absurd h0 (hno y hy)

end
end PFSignChange

#print axioms PFSignChange.exists_zero_between_of_mul_pos
#print axioms PFSignChange.mul_neg_of_no_zero_between
