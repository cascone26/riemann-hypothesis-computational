import EigenbasisOrthonormal
import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic

/-!
# §1735 — Two-mode Parseval: the archimedean metric is diagonal and positive on the eigenbasis

Combining §1734's orthogonality and normalization, the `L²` inner product of any two-mode combination of
interval eigenmodes is `⟨c vₖ + d vₗ, c vₖ + d vₗ⟩ = (c² + d²)·(L/2)`:

`∫₀ᴸ (c·sin(a u) + d·sin(b u))² du = (c² + d²)·(L/2)`   (`a·L=k·π`, `b·L=l·π`, `a ≠ 0,±b`).

The cross term drops (orthogonality) and each square contributes `L/2` (normalization).  So on the eigenbasis
the archimedean metric is diagonal with positive entries `L/2 > 0` — it is **positive definite** (the Parseval
sum `(c²+d²)L/2` is `0` iff `c=d=0` when `L>0`).  This is exactly the positive-definite inner-product
hypothesis the Rosati critical-circle brick §1729 (`⟪F x, F y⟫ = q⟪x,y⟫ ⟹ |λ|=√q`) requires — now grounded on
the concrete archimedean L² pairing rather than assumed.  EXIT0 / std-axioms.
-/

open Real intervalIntegral MeasureTheory

/-- **Two-mode Parseval.**  `∫₀ᴸ (c sin(a u) + d sin(b u))² = (c²+d²)(L/2)`. -/
theorem two_mode_parseval (a b L c d : ℝ) (k l : ℤ)
    (ha : a ≠ 0) (hb : b ≠ 0) (hd : a - b ≠ 0) (hs : a + b ≠ 0)
    (hLa : a * L = (k : ℝ) * π) (hLb : b * L = (l : ℝ) * π) :
    (∫ u in (0:ℝ)..L, (c * Real.sin (a * u) + d * Real.sin (b * u)) ^ 2)
      = (c ^ 2 + d ^ 2) * (L / 2) := by
  have Iaa := eigenmode_norm a L k ha hLa
  have Ibb := eigenmode_norm b L l hb hLb
  have Iab := eigenbasis_orthogonal a b L k l hd hs hLa hLb
  have iAB : IntervalIntegrable
      (fun u => c ^ 2 * (Real.sin (a * u) * Real.sin (a * u))
                  + 2 * c * d * (Real.sin (a * u) * Real.sin (b * u))) volume 0 L :=
    (by fun_prop : Continuous _).intervalIntegrable 0 L
  have iC : IntervalIntegrable
      (fun u => d ^ 2 * (Real.sin (b * u) * Real.sin (b * u))) volume 0 L :=
    (by fun_prop : Continuous _).intervalIntegrable 0 L
  have iA : IntervalIntegrable
      (fun u => c ^ 2 * (Real.sin (a * u) * Real.sin (a * u))) volume 0 L :=
    (by fun_prop : Continuous _).intervalIntegrable 0 L
  have iB : IntervalIntegrable
      (fun u => 2 * c * d * (Real.sin (a * u) * Real.sin (b * u))) volume 0 L :=
    (by fun_prop : Continuous _).intervalIntegrable 0 L
  rw [intervalIntegral.integral_congr
        (g := fun u => c ^ 2 * (Real.sin (a * u) * Real.sin (a * u))
                        + 2 * c * d * (Real.sin (a * u) * Real.sin (b * u))
                        + d ^ 2 * (Real.sin (b * u) * Real.sin (b * u)))
        (fun u _ => by ring)]
  rw [intervalIntegral.integral_add iAB iC, intervalIntegral.integral_add iA iB,
      intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
      intervalIntegral.integral_const_mul, Iaa, Iab, Ibb]
  ring
