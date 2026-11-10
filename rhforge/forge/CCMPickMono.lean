import Mathlib

/-!
# S3009 — Pick generators are monotone between poles (lane 2: the necessary condition behind S3007's numerics)

For g(ν) = Σ_j w_j/(ν − x_j) with w ≥ 0:
* `cauchy_term_antitone` — a ≤ b on the same side of x ((a−x)(b−x) > 0) ⇒ 1/(b−x) ≤ 1/(a−x).
* `pick_generator_antitone` — a ≤ b with every pole x_j outside [a,b] (same side) and w ≥ 0 ⇒ g(b) ≤ g(a).
So if all far zeros are real (w_γ ≥ 0), the far-tail generator g_tail = πΣ w(1/(ν−γ)+1/(ν+γ)) is DECREASING on any window
free of zeros — exactly what S3007 observed for the ARITHMETIC generator (−g_tail/π increasing, min step 2.5e-2–3.2e-2).
A violation (an increase of g_tail on a zero-free window) would certify a non-real zero above T (Report: contrapositive).
-/

open Finset

namespace CCMPickMono

theorem cauchy_term_antitone {a b x : ℝ} (hab : a ≤ b) (hside : 0 < (a - x) * (b - x)) :
    1 / (b - x) ≤ 1 / (a - x) := by
  have ha : a - x ≠ 0 := by intro h; rw [h, zero_mul] at hside; exact lt_irrefl _ hside
  have hb : b - x ≠ 0 := by intro h; rw [h, mul_zero] at hside; exact lt_irrefl _ hside
  have key : 1 / (a - x) - 1 / (b - x) = (b - a) / ((a - x) * (b - x)) := by
    field_simp; ring
  have : 0 ≤ (b - a) / ((a - x) * (b - x)) := div_nonneg (by linarith) hside.le
  linarith

theorem pick_generator_antitone {m : Type*} [Fintype m] (w xz : m → ℝ) (hw : ∀ j, 0 ≤ w j)
    {a b : ℝ} (hab : a ≤ b) (hside : ∀ j, 0 < (a - xz j) * (b - xz j)) :
    ∑ j, w j * (1 / (b - xz j)) ≤ ∑ j, w j * (1 / (a - xz j)) :=
  Finset.sum_le_sum (fun j _ => mul_le_mul_of_nonneg_left (cauchy_term_antitone hab (hside j)) (hw j))

end CCMPickMono

#print axioms CCMPickMono.cauchy_term_antitone
#print axioms CCMPickMono.pick_generator_antitone
