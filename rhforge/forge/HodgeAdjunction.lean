import HodgeIndexCore
import Mathlib.LinearAlgebra.BilinearForm.Orthogonal
import Mathlib.LinearAlgebra.QuadraticForm.Basic

set_option maxHeartbeats 1000000

open QuadraticMap Submodule

/-!
# §1659 — COORDINATE-FREE ARITHMETIC HODGE INDEX from the adjunction relations

(Closes the long-scoped "§1026 target" from STATUS.md.) On a real Néron–Severi-type space with a
symmetric pairing `B`, the fiber/section adjunction relations — `B f f = 0` (fiber self-intersection),
`B s s = 0`, `B f s = 1` (a section meets a fiber once) — plus negative-definiteness on the primitive
complement `{f,s}^⊥`, force signature `(1, dim − 1)`: the Hodge-index shape. Strictly more general
than the coordinatized `(ℝ×ℝ)×V` model (§1023 arithB): applies directly to ANY real arithmetic-surface
NS-space presented by an actual intersection pairing. Positive line: `span{f+s}` (`Q = 2`);
negative complement: `span{f−s} ⊔ {f,s}^⊥` (`Q = −2a² + B w w`).
-/

variable {F M : Type*} [Field F] [LinearOrder F] [IsStrictOrderedRing F]
  [AddCommGroup M] [Module F M] [FiniteDimensional F M]

theorem hodge_signature_of_adjunction (B : LinearMap.BilinForm F M)
    (hsymm : ∀ x y, B x y = B y x)
    (f s : M) (hff : B f f = 0) (hss : B s s = 0) (hfs : B f s = 1)
    (hindep : LinearIndependent F ![f, s])
    (hneg : ∀ w ∈ B.orthogonal (Submodule.span F {f, s}), w ≠ 0 → B w w < 0) :
    sigPos B.toQuadraticMap = 1 ∧
      sigNeg B.toQuadraticMap = Module.finrank F M - 1 := by
  have hsf : B s f = 1 := by rw [hsymm]; exact hfs
  set W : Submodule F M := Submodule.span F {f, s} with hWdef
  have hfW : f ∈ W := subset_span (by simp)
  have hsW : s ∈ W := subset_span (by simp)
  have hrefl : B.IsRefl := fun x y h => by rw [hsymm y x]; exact h
  -- nonzero combinations from independence
  have hne_add : f + s ≠ 0 := by
    intro h
    have := (LinearIndependent.pair_iff.mp hindep 1 1 (by simpa using h)).1
    exact one_ne_zero this
  have hne_sub : f - s ≠ 0 := by
    intro h
    have := (LinearIndependent.pair_iff.mp hindep 1 (-1)
      (by rw [one_smul, neg_smul, one_smul, ← sub_eq_add_neg]; exact h)).1
    exact one_ne_zero this
  -- B restricted to W = span{f,s} is nondegenerate (hyperbolic plane)
  have hsepL : (B.restrict W).SeparatingLeft := by
    intro x hx
    obtain ⟨a, b, hab⟩ := Submodule.mem_span_pair.mp x.2
    have h1 := hx ⟨f, hfW⟩
    have h2 := hx ⟨s, hsW⟩
    rw [LinearMap.BilinForm.restrict_apply] at h1 h2
    have hb : b = 0 := by
      rw [← hab] at h1
      simpa [map_add, map_smul, LinearMap.add_apply, LinearMap.smul_apply, smul_eq_mul,
        hff, hsf] using h1
    have ha : a = 0 := by
      rw [← hab] at h2
      simpa [map_add, map_smul, LinearMap.add_apply, LinearMap.smul_apply, smul_eq_mul,
        hfs, hss] using h2
    apply Subtype.ext
    rw [← hab, ha, hb]
    simp
  have hrestrict : (B.restrict W).Nondegenerate := by
    refine ⟨hsepL, fun y hy => hsepL y fun z => ?_⟩
    have hz := hy z
    simp only [LinearMap.BilinForm.restrict_apply, LinearMap.domRestrict_apply] at hz ⊢
    rw [hsymm]
    exact hz
  -- complement structure and dimensions
  have hcompl : IsCompl W (B.orthogonal W) :=
    LinearMap.BilinForm.isCompl_orthogonal_of_restrict_nondegenerate hrefl hrestrict
  have hWrank : Module.finrank F W = 2 := by
    rw [hWdef, show ({f, s} : Set M) = Set.range ![f, s] from
      (Matrix.range_cons_cons_empty f s ![]).symm]
    simpa using finrank_span_eq_card hindep
  have hsumrank := Submodule.finrank_add_eq_of_isCompl hcompl
  have horthrank : Module.finrank F (B.orthogonal W) = Module.finrank F M - 2 := by omega
  have hMge : 2 ≤ Module.finrank F M := by omega
  -- the positive line and the negative hyperplane
  set P : Submodule F M := Submodule.span F {f + s} with hPdef
  set N : Submodule F M := Submodule.span F {f - s} ⊔ B.orthogonal W with hNdef
  have hPrank : Module.finrank F P = 1 := finrank_span_singleton hne_add
  have hdisj : Disjoint (Submodule.span F {f - s}) (B.orthogonal W) := by
    have hsub : Submodule.span F {f - s} ≤ W :=
      span_le.mpr (Set.singleton_subset_iff.mpr (W.sub_mem hfW hsW))
    exact hcompl.disjoint.mono_left hsub
  have hNrank : Module.finrank F N = Module.finrank F M - 1 := by
    have h := Submodule.finrank_sup_add_finrank_inf_eq (Submodule.span F {f - s})
      (B.orthogonal W)
    rw [disjoint_iff.mp hdisj, finrank_bot, finrank_span_singleton hne_sub, ← hNdef] at h
    omega
  -- positivity on P
  have hPpos : ((B.toQuadraticMap).restrict P).PosDef := by
    intro x hx0
    obtain ⟨a, ha⟩ := Submodule.mem_span_singleton.mp x.2
    have ha0 : a ≠ 0 := by
      intro h
      apply hx0
      apply Subtype.ext
      rw [← ha, h, zero_smul]
      simp
    have hQ : B x.1 x.1 = 2 * (a * a) := by
      rw [← ha]
      simp only [map_add, map_smul, LinearMap.add_apply, LinearMap.smul_apply, smul_eq_mul,
        hff, hss, hfs, hsf]
      ring
    rw [restrict_apply, LinearMap.BilinMap.toQuadraticMap_apply, hQ]
    have : 0 < a * a := mul_self_pos.mpr ha0
    linarith
  -- negativity on N
  have hNneg : (((-B.toQuadraticMap)).restrict N).PosDef := by
    intro x hx0
    obtain ⟨y, hy, w, hw, hyw⟩ := Submodule.mem_sup.mp x.2
    obtain ⟨a, ha⟩ := Submodule.mem_span_singleton.mp hy
    have hBfw : B f w = 0 := hw f hfW
    have hBsw : B s w = 0 := hw s hsW
    have hBwf : B w f = 0 := by rw [hsymm]; exact hBfw
    have hBws : B w s = 0 := by rw [hsymm]; exact hBsw
    have hQ : B x.1 x.1 = -(2 * (a * a)) + B w w := by
      rw [← hyw, ← ha]
      simp only [map_add, map_sub, map_smul, LinearMap.add_apply, LinearMap.sub_apply,
        LinearMap.smul_apply, smul_eq_mul, hff, hss, hfs, hsf, hBfw, hBsw, hBwf, hBws]
      ring
    rw [restrict_apply, neg_apply, LinearMap.BilinMap.toQuadraticMap_apply, hQ]
    by_cases hw0 : w = 0
    · have ha0 : a ≠ 0 := by
        intro h
        apply hx0
        apply Subtype.ext
        rw [← hyw, ← ha, h, zero_smul, hw0, add_zero]
        simp
      have hww : B w w = 0 := by rw [hw0]; simp
      have : 0 < a * a := mul_self_pos.mpr ha0
      rw [hww]
      linarith
    · have hww : B w w < 0 := hneg w hw hw0
      have : 0 ≤ a * a := mul_self_nonneg a
      linarith
  exact hodgeIndex_signature B.toQuadraticMap P N hPrank hPpos hNneg hNrank (by omega)

#print axioms hodge_signature_of_adjunction
