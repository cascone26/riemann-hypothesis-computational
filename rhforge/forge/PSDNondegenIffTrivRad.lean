import PSDNullRadical

/-!
# §1749 — Nondegeneracy ⟺ trivial radical (PSD forms)

The general-principle mirror of `NeronTateRadical` (§1089, which proves the same statement concretely for the
Néron–Tate height: nondegenerate exactly modulo torsion).  For a positive-semidefinite symmetric bilinear form
the two natural notions of "no degeneracy" coincide:

`(∀ u, (∀ w, B u w = 0) → u = 0)  ⟺  (∀ u, B u u = 0 → u = 0)`.

Left = nondegenerate (trivial full radical); right = anisotropic (no nonzero null vector).  For a general form
left ⟸ right needs nothing, and left ⟹ right needs §1744 (a null vector pairs to zero with everything).  With
§1748 this closes the triple `nondegenerate ⟺ trivial radical ⟺ positive definite` for the archimedean-type
PSD forms — the abstract skeleton of "the Weil radical is trivial ⟺ RH".  EXIT0 / std-axioms.
-/

/-- **Nondegenerate ⟺ anisotropic** for a PSD symmetric bilinear form. -/
theorem psd_nondegenerate_iff_trivial_radical {M : Type*} [AddCommGroup M] [Module ℝ M]
    (B : LinearMap.BilinForm ℝ M) (hsymm : ∀ x y, B x y = B y x) (hpsd : ∀ x, 0 ≤ B x x) :
    (∀ u, (∀ w, B u w = 0) → u = 0) ↔ (∀ u, B u u = 0 → u = 0) := by
  constructor
  · intro hnd u hu
    exact hnd u (fun w => psd_null_in_radical B hsymm hpsd u hu w)
  · intro htr u hu
    exact htr u (hu u)
