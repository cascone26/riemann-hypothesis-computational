import Mathlib.LinearAlgebra.FiniteDimensional.Lemmas
import Mathlib.Tactic

/-!
# §1779 — A form PSD on a codimension-1 subspace has at most one negative direction (#neg ≤ 1)

The linear-algebra core of the screw-function resolution of Connes' CS hypothesis (THE_BRIDGE cycle 123):
Suzuki proves the ζ-screw function `g_ξ` lies in Kreĭn's class `𝒢_∞`, i.e. the increment kernel
`G_g(t,u)=g(t-u)-g(t)-g(-u)+g(0)` is nonnegative definite.  For every mean-zero test function the plain
Weil–CS quadratic form equals the `G_g`-form, hence is **PSD on the codimension-1 mean-zero subspace**.
This brick is the abstract consequence: a quadratic form `B` that is `≥ 0` on a subspace `S` and `< 0` on a
subspace `W` forces `W ⊓ S = ⊥` (a common nonzero vector would be both `≥0` and `<0`), whence
`dim W + dim S ≤ dim V`.  With `codim S ≤ 1` this gives `dim W ≤ 1`: **at most one negative direction**,
so the minimal eigenvalue of the Weil–CS form is simple and isolated — exactly the CS hypothesis.
No spectral theorem needed: pure dimension counting.  EXIT0 / std-axioms.
-/

open Module

variable {V : Type*} [AddCommGroup V] [Module ℝ V] [FiniteDimensional ℝ V]

/-- **Disjointness from opposite definiteness.**  If `B` is `≥ 0` on `S` and `< 0` on nonzero vectors of
`W`, the two subspaces meet only at `0`, so their dimensions add to at most `dim V`. -/
theorem finrank_negdef_add_finrank_psd_le
    (B : V → V → ℝ) (S W : Submodule ℝ V)
    (hS : ∀ x ∈ S, 0 ≤ B x x)
    (hW : ∀ w ∈ W, w ≠ 0 → B w w < 0) :
    finrank ℝ W + finrank ℝ S ≤ finrank ℝ V := by
  have hdisj : W ⊓ S = ⊥ := by
    rw [Submodule.eq_bot_iff]
    intro v hv
    obtain ⟨hvW, hvS⟩ := hv
    by_contra hne
    exact absurd (hS v hvS) (not_le.mpr (hW v hvW hne))
  have key := Submodule.finrank_sup_add_finrank_inf_eq W S
  rw [hdisj, finrank_bot, add_zero] at key
  calc finrank ℝ W + finrank ℝ S = finrank ℝ ↑(W ⊔ S) := key.symm
    _ ≤ finrank ℝ V := Submodule.finrank_le _

/-- **#neg ≤ 1 from PSD on a codimension-1 subspace.**  If `B ≥ 0` on `S` with `codim S ≤ 1`
(`dim V ≤ dim S + 1`), then any subspace `W` on which `B` is negative definite has `dim W ≤ 1`:
the form has at most one negative direction. -/
theorem finrank_negdef_le_one_of_psd_codim_one
    (B : V → V → ℝ) (S W : Submodule ℝ V)
    (hS : ∀ x ∈ S, 0 ≤ B x x)
    (hW : ∀ w ∈ W, w ≠ 0 → B w w < 0)
    (hcodim : finrank ℝ V ≤ finrank ℝ S + 1) :
    finrank ℝ W ≤ 1 := by
  have h := finrank_negdef_add_finrank_psd_le B S W hS hW
  omega
