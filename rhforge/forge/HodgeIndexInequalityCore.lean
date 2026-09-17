import HodgeIndexCore

set_option maxHeartbeats 2000000

open QuadraticMap Submodule

variable {F M : Type*} [Field F] [LinearOrder F] [IsStrictOrderedRing F]
  [AddCommGroup M] [Module F M] [FiniteDimensional F M]

private lemma Q_add_smul (Q : QuadraticForm F M) (a b : F) (u v : M) :
    Q (a • u + b • v) = a * a * Q u + b * b * Q v + a * b * polar Q u v := by
  rw [QuadraticMap.map_add Q (a • u) (b • v), QuadraticMap.map_smul, QuadraticMap.map_smul,
      QuadraticMap.polar_smul_left, QuadraticMap.polar_smul_right]
  simp only [smul_eq_mul]
  ring

/-- **Hodge index inequality (abstract core).** If `Q` has exactly one positive square and `Q x > 0`,
then `(polar Q x y)^2 ≥ 4 * Q x * Q y` for all `y` — the Lorentzian reverse-Cauchy–Schwarz that (in the
function-field case) yields RH. Over any ordered field. -/
theorem hodge_index_inequality
    (Q : QuadraticForm F M) (hsig : sigPos Q = 1)
    (x y : M) (hx : 0 < Q x) :
    4 * Q x * Q y ≤ polar Q x y ^ 2 := by
  set p : F := polar Q x y with hp
  set z : M := (2 * Q x) • y + (-p) • x with hzdef
  have hxz : polar Q x z = 0 := by
    simp only [hzdef, QuadraticMap.polar_add_right, QuadraticMap.polar_smul_right,
      QuadraticMap.polar_self, smul_eq_mul, ← hp]
    ring
  have hQz : Q z = Q x * (4 * Q x * Q y - p ^ 2) := by
    rw [hzdef, Q_add_smul, QuadraticMap.polar_comm Q y x, ← hp]; ring
  by_contra hcon
  push_neg at hcon
  have hQzpos : 0 < Q z := by rw [hQz]; apply mul_pos hx; linarith
  have hxz_indep : LinearIndependent F ![x, z] := by
    rw [LinearIndependent.pair_iff]
    intro s t hst
    have h1 : polar Q x (s • x + t • z) = 0 := by rw [hst]; exact QuadraticMap.polar_zero_right _ _
    rw [QuadraticMap.polar_add_right, QuadraticMap.polar_smul_right,
        QuadraticMap.polar_smul_right, QuadraticMap.polar_self, hxz] at h1
    simp only [smul_eq_mul, mul_zero, add_zero] at h1
    have hQx2 : (2 : F) * Q x ≠ 0 := by positivity
    have hs0 : s = 0 := by
      have hz2 : s * (2 * Q x) = 0 := by ring_nf; ring_nf at h1; linarith [h1]
      rcases mul_eq_zero.mp hz2 with h | h
      · exact h
      · exact absurd h hQx2
    refine ⟨hs0, ?_⟩
    rw [hs0, zero_smul, zero_add] at hst
    rcases smul_eq_zero.mp hst with h | h
    · exact h
    · exfalso; rw [h, QuadraticMap.map_zero] at hQzpos; exact lt_irrefl 0 hQzpos
  set W : Submodule F M := span F {x, z} with hWdef
  have hWpos : (Q.restrict W).PosDef := by
    intro w hw
    have hwmem : (w : M) ∈ span F {x, z} := w.2
    rw [Submodule.mem_span_pair] at hwmem
    obtain ⟨s, t, hst⟩ := hwmem
    have hval : (Q.restrict W) w = s * s * Q x + t * t * Q z := by
      have h1 : (Q.restrict W) w = Q (w : M) := restrict_apply _ _ _
      rw [h1, ← hst, Q_add_smul, hxz]; ring
    rw [hval]
    have hst0 : ¬ (s = 0 ∧ t = 0) := by
      rintro ⟨rfl, rfl⟩
      apply hw; apply Subtype.ext
      simp only [ZeroMemClass.coe_zero, ← hst, zero_smul, add_zero]
    rcases not_and_or.mp hst0 with hs | ht
    · have hss : 0 < s * s := (mul_self_nonneg s).lt_of_ne (Ne.symm (mul_ne_zero hs hs))
      nlinarith [mul_pos hss hx, mul_nonneg (mul_self_nonneg t) hQzpos.le]
    · have htt : 0 < t * t := (mul_self_nonneg t).lt_of_ne (Ne.symm (mul_ne_zero ht ht))
      nlinarith [mul_pos htt hQzpos, mul_nonneg (mul_self_nonneg s) hx.le]
  have hfr : Module.finrank F W = 2 := by
    have hW2 : W = span F (Set.range ![x, z]) := by
      rw [hWdef]; congr 1; ext q
      simp only [Set.mem_insert_iff, Set.mem_singleton_iff, Set.mem_range,
        Fin.exists_fin_two, Matrix.cons_val_zero, Matrix.cons_val_one]
      tauto
    rw [hW2, finrank_span_eq_card hxz_indep, Fintype.card_fin]
  have hle := le_sigPos_of_posDef Q hWpos
  rw [hfr, hsig] at hle
  omega

#print axioms hodge_index_inequality
