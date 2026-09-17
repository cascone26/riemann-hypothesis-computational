import HodgeAdjunction
import HyperbolicPairComplexPlace

set_option maxHeartbeats 2000000

open IsDedekindDomain NumberField QuadraticMap Submodule

noncomputable section

variable (K : Type*) [Field K] [NumberField K]

namespace ArakelovVecR

/-!
# §1716 ★★★ — §1659 INSTANTIATED: the Hodge signature (1,1) on the model's hyperbolic plane

The wiring proof-of-life: restrict `hodgeBilin` to `W = span{F,S}` for any hyperbolic pair.  On
`W` the pair GENERATES, so the orthogonal complement of `span{f,s}` is forced to zero
(`B(f,w) = b`, `B(s,w) = a` for `w = af+bs`) and §1659's negativity hypothesis holds vacuously;
`finrank W = 2` from independence; `hodge_signature_of_adjunction` returns **signature (1,1)** —
the arithmetic Hodge index shape, machine-checked end-to-end on the concrete Arakelov model of
every number field (§1715 supplies the pair).
-/

theorem hodge_signature_on_hyperbolic_plane
    {F S : ArakelovVecR K}
    (hFF : hodgeBilin K F F = 0) (hSS : hodgeBilin K S S = 0)
    (hFS : hodgeBilin K F S = 1) (hind : LinearIndependent ℝ ![F, S]) :
    sigPos (R := ℝ) (M := span ℝ ({F, S} : Set (ArakelovVecR K)))
      (((hodgeBilin K).restrict
        (span ℝ ({F, S} : Set (ArakelovVecR K)))).toQuadraticMap) = 1 ∧
    sigNeg (R := ℝ) (M := span ℝ ({F, S} : Set (ArakelovVecR K)))
      (((hodgeBilin K).restrict
        (span ℝ ({F, S} : Set (ArakelovVecR K)))).toQuadraticMap) = 1 := by
  classical
  set W : Submodule ℝ (ArakelovVecR K) := span ℝ ({F, S} : Set (ArakelovVecR K)) with hWdef
  haveI : FiniteDimensional ℝ W := FiniteDimensional.span_of_finite ℝ (Set.toFinite _)
  set B' := (hodgeBilin K).restrict W with hB'def
  have hFW : F ∈ W := subset_span (by simp)
  have hSW : S ∈ W := subset_span (by simp)
  set f : W := ⟨F, hFW⟩ with hfdef
  set s : W := ⟨S, hSW⟩ with hsdef
  have hres : ∀ x y : W, B' x y = hodgeBilin K (x : ArakelovVecR K) (y : ArakelovVecR K) :=
    fun x y => rfl
  have hsymm' : ∀ x y : W, B' x y = B' y x := by
    intro x y
    rw [hres, hres, hodgeBilin_symm]
  have hff' : B' f f = 0 := by rw [hres]; exact hFF
  have hss' : B' s s = 0 := by rw [hres]; exact hSS
  have hfs' : B' f s = 1 := by rw [hres]; exact hFS
  have hsf' : B' s f = 1 := by rw [hsymm']; exact hfs'
  -- independence upstairs
  have hind' : LinearIndependent ℝ ![f, s] := by
    apply LinearIndependent.of_comp W.subtype
    have hcomp : (W.subtype ∘ ![f, s]) = ![F, S] := by
      funext i
      fin_cases i <;> rfl
    rw [hcomp]
    exact hind
  -- the pair generates W
  have hspan_top : span ℝ ({f, s} : Set W) = ⊤ := by
    rw [Submodule.eq_top_iff']
    intro w
    have hw : (w : ArakelovVecR K) ∈ span ℝ ({F, S} : Set (ArakelovVecR K)) := w.2
    rw [Submodule.mem_span_pair] at hw
    obtain ⟨a, b, hab⟩ := hw
    rw [Submodule.mem_span_pair]
    refine ⟨a, b, Subtype.ext ?_⟩
    simpa using hab
  -- the orthogonal complement of the generating pair is zero → hneg vacuous
  have hneg : ∀ w ∈ B'.orthogonal (span ℝ ({f, s} : Set W)), w ≠ 0 → B' w w < 0 := by
    intro w hw hwne
    exfalso
    rw [LinearMap.BilinForm.mem_orthogonal_iff] at hw
    have h1 : B' f w = 0 := hw f (subset_span (by simp))
    have h2 : B' s w = 0 := hw s (subset_span (by simp))
    have hwmem : w ∈ span ℝ ({f, s} : Set W) := by rw [hspan_top]; trivial
    rw [Submodule.mem_span_pair] at hwmem
    obtain ⟨a, b, hab⟩ := hwmem
    have hBfw : B' f w = b := by
      rw [← hab, map_add, map_smul, map_smul, smul_eq_mul, smul_eq_mul,
        hff', hfs']
      ring
    have hBsw : B' s w = a := by
      rw [← hab, map_add, map_smul, map_smul, smul_eq_mul, smul_eq_mul,
        hsf', hss']
      ring
    have hb0 : b = 0 := by rw [← hBfw]; exact h1
    have ha0 : a = 0 := by rw [← hBsw]; exact h2
    apply hwne
    rw [← hab, ha0, hb0]
    simp
  -- finrank W = 2
  have hrange : ({F, S} : Set (ArakelovVecR K)) = Set.range ![F, S] := by
    ext x
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff, Set.mem_range,
      Fin.exists_fin_two, Matrix.cons_val_zero, Matrix.cons_val_one]
    tauto
  have hfr : Module.finrank ℝ W = 2 := by
    have h2 : W = span ℝ (Set.range ![F, S]) := by rw [hWdef, hrange]
    rw [h2, finrank_span_eq_card hind]
    simp
  obtain ⟨hpos, hnegdim⟩ :=
    hodge_signature_of_adjunction B' hsymm' f s hff' hss' hfs' hind' hneg
  refine ⟨hpos, ?_⟩
  rw [hnegdim, hfr]

/-- **The arithmetic Hodge signature (1,1) lives in the Arakelov model of every number field**:
combine any prime with §1715's hyperbolic pair. -/
theorem exists_hodge_signature_plane (v : HeightOneSpectrum (RingOfIntegers K)) :
    ∃ F S : ArakelovVecR K,
      sigPos (R := ℝ) (M := span ℝ ({F, S} : Set (ArakelovVecR K)))
        (((hodgeBilin K).restrict
          (span ℝ ({F, S} : Set (ArakelovVecR K)))).toQuadraticMap) = 1 ∧
      sigNeg (R := ℝ) (M := span ℝ ({F, S} : Set (ArakelovVecR K)))
        (((hodgeBilin K).restrict
          (span ℝ ({F, S} : Set (ArakelovVecR K)))).toQuadraticMap) = 1 := by
  obtain ⟨F, S, hFF, hSS, hFS, hind⟩ := exists_hyperbolic_pair_any K v
  exact ⟨F, S, hodge_signature_on_hyperbolic_plane K hFF hSS hFS hind⟩

#print axioms hodge_signature_on_hyperbolic_plane
#print axioms exists_hodge_signature_plane
