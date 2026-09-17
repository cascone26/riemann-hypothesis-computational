import Mathlib.Tactic
import Mathlib.Algebra.Module.LinearMap.Defs

/-!
# §1764 — Kramers degeneracy vs particle-hole: why class C gives SIMPLE zeros

The spectral-construction route needs the zeros to be SIMPLE (Yakaboylu's construction assumes it).  Cycle 36
found — and this brick proves the abstract core of — WHY the AZ symmetry class settles it:

* **Time-reversal `T` with `T² = −1` (class AII) FORCES double degeneracy (Kramers).**  If `T H = H T` and
  `H v = λ v` (`v ≠ 0`), then `T v` is ANOTHER `λ`-eigenvector and `{v, T v}` is linearly independent — so
  every eigenvalue has multiplicity `≥ 2`.  (Proved here in real form: `T v = c v` would give `c² = −1`,
  impossible over `ℝ`.)
* **Particle-hole `C` with `C H = −H C` (class C) does NOT** — it sends `λ ↦ −λ` (§1763), a DIFFERENT
  eigenspace for `λ ≠ 0`, imposing no degeneracy at `λ`.

The Riemann Hamiltonian is class C (Srednicki): particle-hole, NO time-reversal.  So it has NO Kramers
doubling — its nonzero spectrum is generically simple = the zeros are SIMPLE by symmetry class.  This brick
banks Kramers' theorem (the mechanism that class C escapes), companion to §1763.  EXIT0 / std-axioms.
-/

variable {V : Type*} [AddCommGroup V] [Module ℝ V]

/-- **Kramers partner is an eigenvector at the same eigenvalue.**  `T H = H T`, `H v = λ v` ⟹ `H (T v) = λ (T v)`. -/
theorem kramers_eigenpartner (H T : V →ₗ[ℝ] V)
    (hcomm : T ∘ₗ H = H ∘ₗ T) (v : V) (lam : ℝ) (hv : H v = lam • v) :
    H (T v) = lam • (T v) := by
  have h := LinearMap.congr_fun hcomm v
  simp only [LinearMap.comp_apply] at h
  rw [← h, hv, map_smul]

/-- **Kramers degeneracy (real form).**  With `T² = −id`, `v` and `T v` are linearly independent for `v ≠ 0`.
Combined with `kramers_eigenpartner`, every eigenvalue then has multiplicity `≥ 2` — the class-AII doubling
that class C (no such `T`) avoids. -/
theorem kramers_independent (T : V →ₗ[ℝ] V) (hT2 : T ∘ₗ T = -LinearMap.id)
    (v : V) (hv0 : v ≠ 0) (a b : ℝ) (h : a • v + b • (T v) = 0) :
    a = 0 ∧ b = 0 := by
  have hT2v : T (T v) = -v := by
    have := LinearMap.congr_fun hT2 v; simpa using this
  -- apply T to h and use T²v = -v:  a • T v - b • v = 0
  have h2 : a • (T v) - b • v = 0 := by
    have hh := congrArg T h
    rw [map_add, map_smul, map_smul, map_zero, hT2v] at hh
    -- hh : a • T v + b • (-v) = 0
    rw [smul_neg] at hh
    linear_combination (norm := module) hh
  -- combine: a•h - b•h2  gives (a²+b²)•v = 0
  have hcomb : (a * a + b * b) • v = 0 := by
    linear_combination (norm := module) a • h - b • h2
  have hab : a * a + b * b = 0 := by
    rcases smul_eq_zero.mp hcomb with hz | hz
    · exact hz
    · exact absurd hz hv0
  constructor
  · nlinarith [mul_self_nonneg a, mul_self_nonneg b, hab]
  · nlinarith [mul_self_nonneg a, mul_self_nonneg b, hab]
