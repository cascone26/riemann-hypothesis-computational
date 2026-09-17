import Mathlib.Analysis.InnerProductSpace.Basic
import SelfInversiveEigenpoly

/-!
# §1777 — Self-inversive + no zeros inside the disk ⟹ all zeros on |z|=1 (on-circle core)

The on-circle half of the Connes–van Suijlekom / Carathéodory-Fejér engine, reduced to an ELEMENTARY core that
avoids the full trigonometric-moment theory.  By §1776 the (palindromic/even) eigenvector polynomial is
self-inversive: every zero `ρ ≠ 0` is paired with `1/ρ`.  The remaining Carathéodory content — that the singular
PSD Toeplitz structure puts NO zero strictly inside the open unit disk — then forces every zero onto the unit
circle: if `‖ρ‖ ≥ 1` and its partner `‖1/ρ‖ ≥ 1`, then `‖ρ‖ = 1`.

So the CS engine's on-circle conclusion factors as: [§1776 reciprocal pairing] + [PSD ⟹ no interior zeros]
⟹ zeros on `|z|=1` (= critical line).  This brick discharges the first-plus-geometry step; the owed analytic
input is exactly "no zeros inside the disk" (the Schur/PSD condition).  EXIT0 / std-axioms.
-/

open Finset

variable {n : ℕ}

/-- **On-circle from self-inversive + no interior zeros.**  If the coefficient vector is palindromic (§1776) and
the polynomial `∑ vᵢ zⁱ` has no zero strictly inside the unit disk (`‖z‖ ≥ 1` at every nonzero root), then every
nonzero root lies exactly on the unit circle. -/
theorem self_inversive_on_circle (v : Fin n → ℂ) (hpal : ∀ i, v (Fin.rev i) = v i)
    (hno_inside : ∀ z : ℂ, z ≠ 0 → (∑ i : Fin n, v i * z ^ (i : ℕ)) = 0 → 1 ≤ ‖z‖)
    (ρ : ℂ) (hρ : ρ ≠ 0) (hroot : (∑ i : Fin n, v i * ρ ^ (i : ℕ)) = 0) :
    ‖ρ‖ = 1 := by
  -- 1/ρ is also a root (self-inversive, §1776)
  have hrec : (∑ i : Fin n, v i * (ρ⁻¹) ^ (i : ℕ)) = 0 :=
    self_inversive_reciprocal_zero v hpal ρ hρ hroot
  have hinvne : ρ⁻¹ ≠ 0 := inv_ne_zero hρ
  -- both ρ and 1/ρ have norm ≥ 1
  have h2 : 1 ≤ ‖ρ‖ := hno_inside ρ hρ hroot
  have h1 : 1 ≤ ‖ρ⁻¹‖ := hno_inside ρ⁻¹ hinvne hrec
  rw [norm_inv] at h1
  -- ‖ρ‖ ≥ 1 and ‖ρ‖⁻¹ ≥ 1 ⟹ ‖ρ‖ = 1
  have hpos : (0 : ℝ) < ‖ρ‖ := norm_pos_iff.mpr hρ
  have hmul : ‖ρ‖ * ‖ρ‖⁻¹ = 1 := mul_inv_cancel₀ (ne_of_gt hpos)
  have hle : ‖ρ‖ ≤ 1 := by
    nlinarith [h1, hpos, hmul, mul_nonneg (le_of_lt hpos) (sub_nonneg.mpr h1)]
  linarith
