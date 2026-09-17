import QuadSectorSplit

/-!
# §1751 — The radical of a parity-decoupled PSD form splits by parity

Formalises tonight's verified numerical finding (the Λ=0 saturation is balanced across parity sectors): for a
positive-semidefinite parity-decoupled form (`W vₑ vₒ = 0`), a sum is null iff each parity component is null.

`W (vₑ + vₒ) (vₑ + vₒ) = 0  ⟺  W vₑ vₑ = 0 ∧ W vₒ vₒ = 0`.

Proof: `§1740` gives `W (vₑ+vₒ)(vₑ+vₒ) = W vₑ vₑ + W vₒ vₒ`, a sum of two nonnegative terms — zero iff both are.
So `rad(W)` respects the even/odd decomposition; the saturation locus decomposes as `rad(W_even) ⊕ rad(W_odd)`,
exactly the balanced-parity null space measured on the exact Weil form.  EXIT0 / std-axioms.
-/

/-- **Parity-split of the radical.**  A parity-decoupled PSD form is null on a sum iff null on each part. -/
theorem weil_radical_parity_split {M : Type*} [AddCommGroup M] [Module ℝ M]
    (W : LinearMap.BilinForm ℝ M) (hsymm : ∀ x y, W x y = W y x) (hpsd : ∀ x, 0 ≤ W x x)
    (ve vo : M) (hcross : W ve vo = 0) :
    W (ve + vo) (ve + vo) = 0 ↔ (W ve ve = 0 ∧ W vo vo = 0) := by
  rw [quad_sector_split W hsymm ve vo hcross]
  constructor
  · intro h
    exact ⟨le_antisymm (by linarith [hpsd vo]) (hpsd ve),
           le_antisymm (by linarith [hpsd ve]) (hpsd vo)⟩
  · intro ⟨h1, h2⟩
    rw [h1, h2]; ring
