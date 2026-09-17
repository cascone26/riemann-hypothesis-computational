import QuadSectorSplit

/-!
# §1742 — Weil positivity from sector positivity

Combining §1740 (sector split) directly: if the Weil form `W` is parity-decoupled (`W vₑ vₒ = 0`) and is
nonnegative on each parity sector, then it is nonnegative on their sum:

`W vₑ vₑ ≥ 0` and `W vₒ vₒ ≥ 0` (and `W vₑ vₒ = 0`) ⟹ `W (vₑ + vₒ) (vₑ + vₒ) ≥ 0`.

With §1739 (the Weil form's cross term genuinely vanishes for opposite parity), this reduces the full
Weil-positivity problem to two independent sector problems — the even sector and the odd sector — exactly the
block structure the numerics exhibit.  EXIT0 / std-axioms.
-/

/-- **Sector-wise sufficiency.**  Nonnegativity on each parity sector gives nonnegativity on the sum. -/
theorem weil_psd_of_sectors {M : Type*} [AddCommGroup M] [Module ℝ M]
    (W : LinearMap.BilinForm ℝ M) (hsymm : ∀ x y, W x y = W y x)
    (ve vo : M) (hcross : W ve vo = 0) (he : 0 ≤ W ve ve) (ho : 0 ≤ W vo vo) :
    0 ≤ W (ve + vo) (ve + vo) := by
  rw [quad_sector_split W hsymm ve vo hcross]
  linarith
