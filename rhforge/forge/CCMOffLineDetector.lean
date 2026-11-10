import CCMPickMono

/-!
# S3011 — an off-line-zero DETECTOR from monotonicity of the generator (lane 2)

Contrapositive of S3009: for a pole representation g(ν) = Σ_j w_j/(ν − x_j) (weights w_j = ccmW-type, ≥ 0 iff the
corresponding zero is real), an INCREASE of g on a pole-free interval forces a negative weight:
* `exists_neg_weight_of_increase` — a ≤ b, every pole outside [a,b], Σ w/(b−x) > Σ w/(a−x) ⇒ ∃ j, w_j < 0.
* `pick_iff_detector_quiet` — equivalently: all w_j ≥ 0 ⇒ no increase on any pole-free [a,b].
Reading (Report): applied to the arithmetic tail generator (prime side, S2970/S3007), a single certified increase on a
zero-free window would certify an off-line zero above T. S3007 numerics (λ² = 20, 28) show strict decrease — the detector
is quiet, as RH predicts.
-/

open Finset

namespace CCMDetector

theorem exists_neg_weight_of_increase {m : Type*} [Fintype m] (w xz : m → ℝ) {a b : ℝ} (hab : a ≤ b)
    (hside : ∀ j, 0 < (a - xz j) * (b - xz j))
    (hinc : ∑ j, w j * (1 / (a - xz j)) < ∑ j, w j * (1 / (b - xz j))) :
    ∃ j, w j < 0 := by
  by_contra h
  push Not at h
  have := CCMPickMono.pick_generator_antitone w xz h hab hside
  linarith

theorem pick_iff_detector_quiet {m : Type*} [Fintype m] (w xz : m → ℝ) (hw : ∀ j, 0 ≤ w j) {a b : ℝ}
    (hab : a ≤ b) (hside : ∀ j, 0 < (a - xz j) * (b - xz j)) :
    ¬ (∑ j, w j * (1 / (a - xz j)) < ∑ j, w j * (1 / (b - xz j))) := by
  have := CCMPickMono.pick_generator_antitone w xz hw hab hside
  linarith

end CCMDetector

#print axioms CCMDetector.exists_neg_weight_of_increase
#print axioms CCMDetector.pick_iff_detector_quiet
