import HodgeAdjunction
import HodgeHyperbolic

set_option maxHeartbeats 1000000

open QuadraticMap Submodule

/-!
# §1663 — NON-VACUITY of the coordinate-free Hodge index: the hyperbolic-plane witness

§1659's hypotheses are satisfiable: on `F × F` with the hyperbolic form, `f = (1,1)` and
`s = (½,−½)` are an isotropic pair with `B f s = 1` — the adjunction relations of a fiber and a
section, realized concretely — and the primitive complement of `span{f,s} = F²` is trivial, so
negativity holds vacuously. Conclusion: signature `(1,1)`, matching `dim − 1 = 1`. Proves the
§1659 theorem has genuine content (its hypothesis set is consistent), the same discipline as the
original `HodgeHyperbolic` witness for the inequality form.
-/

variable {F : Type*} [Field F] [LinearOrder F] [IsStrictOrderedRing F]

theorem hyperbolic_adjunction_signature :
    sigPos (hyperbolic (F := F)).toQuadraticMap = 1 ∧
      sigNeg (hyperbolic (F := F)).toQuadraticMap = 1 := by
  have h2ne : (2 : F) ≠ 0 := two_ne_zero
  have hindep : LinearIndependent F ![((1 : F), (1 : F)), ((1 / 2 : F), (-(1 / 2) : F))] := by
    rw [LinearIndependent.pair_iff]
    intro a b hab
    have h1 := congrArg Prod.fst hab
    have h2 := congrArg Prod.snd hab
    simp only [Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul,
      Prod.fst_zero, Prod.snd_zero] at h1 h2
    constructor <;> nlinarith [h1, h2]
  have h := hodge_signature_of_adjunction (hyperbolic (F := F))
    (fun x y => by simp only [hyperbolic_apply]; ring)
    ((1 : F), (1 : F)) ((1 / 2 : F), (-(1 / 2) : F))
    (by norm_num [hyperbolic_apply])
    (by norm_num [hyperbolic_apply])
    (by norm_num [hyperbolic_apply])
    hindep
    (by
      intro w hw hw0
      exfalso
      have hf : hyperbolic (F := F) ((1 : F), (1 : F)) w = 0 :=
        hw ((1 : F), (1 : F)) (subset_span (by simp))
      have hs : hyperbolic (F := F) ((1 / 2 : F), (-(1 / 2) : F)) w = 0 :=
        hw ((1 / 2 : F), (-(1 / 2) : F)) (subset_span (by simp))
      simp only [hyperbolic_apply] at hf hs
      apply hw0
      have hw1 : w.1 = 0 := by nlinarith [hf, hs]
      have hw2 : w.2 = 0 := by nlinarith [hf, hs]
      exact Prod.ext hw1 hw2)
  have hrank : Module.finrank F (F × F) = 2 := by
    simp [Module.finrank_prod]
  rw [hrank] at h
  exact h

#print axioms hyperbolic_adjunction_signature
