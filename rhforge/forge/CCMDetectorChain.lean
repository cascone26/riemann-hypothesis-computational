import Mathlib

/-!
# S3017 — bridge from the Arb step certificate (S3015) to a robust grid-quiet statement (lane 2)

The Arb certifier (operators/arb_offline_detector_cert.py) outputs certified lower bounds s_k ≥ s > 0 on the grid steps
(g̃(ν_k) − g̃(ν_{k+1}))/π of the zero-subtracted generator.  This brick turns that certificate format into statements:
* `strictAnti_of_steps` — positive steps on the grid ⇒ g̃ strictly decreasing along the whole grid (chain).
* `drop_ge_of_steps` — steps ≥ s ⇒ g̃(ν_i) − g̃(ν_j) ≥ (j − i)·s for i ≤ j (quantitative margin).
* `quiet_robust` — ROBUSTNESS: if the true generator differs from the certified one by at most η at every grid node
  and 2η < s, the true generator is still strictly decreasing on the grid.
Reading (Report): S3015 certifies s ≥ π·2.1e-2 at N = 24 (λ² ≤ 36), so the quiet verdict survives any error in the
τ formula (CCM Lemma 5.1) below ~3e-2 in g̃ — eight orders above S2859's formula-vs-zero-side agreement (2e-10).
-/

namespace CCMDetectorChain

theorem drop_ge_of_steps (f : ℕ → ℝ) (s : ℝ) (n : ℕ) (hstep : ∀ k, k < n → s ≤ f k - f (k + 1)) :
    ∀ i d, i + d ≤ n → (d : ℝ) * s ≤ f i - f (i + d) := by
  intro i d
  induction d with
  | zero => intro _; simp
  | succ d ih =>
    intro hid
    have h1 := ih (by omega)
    have h2 := hstep (i + d) (by omega)
    rw [show i + (d + 1) = i + d + 1 by omega]
    push_cast
    linarith

theorem strictAnti_of_steps (f : ℕ → ℝ) (n : ℕ) (hstep : ∀ k, k < n → f (k + 1) < f k) :
    ∀ i j, i < j → j ≤ n → f j < f i := by
  intro i j hij hjn
  obtain ⟨d, rfl⟩ : ∃ d, j = i + (d + 1) := ⟨j - i - 1, by omega⟩
  induction d with
  | zero => simpa using hstep i (by omega)
  | succ d ih =>
    have h1 := ih (by omega) (by omega)
    have h2 := hstep (i + (d + 1)) (by omega)
    rw [show i + (d + 1 + 1) = i + (d + 1) + 1 by omega]
    linarith

theorem quiet_robust (f g : ℕ → ℝ) (s η : ℝ) (n : ℕ) (hstep : ∀ k, k < n → s ≤ f k - f (k + 1))
    (herr : ∀ k, k ≤ n → |g k - f k| ≤ η) (hη : 2 * η < s) :
    ∀ i j, i < j → j ≤ n → g j < g i := by
  apply strictAnti_of_steps g n
  intro k hk
  have h1 := hstep k hk
  have h2 := herr k (by omega)
  have h3 := herr (k + 1) (by omega)
  rw [abs_le] at h2 h3
  linarith [h2.1, h2.2, h3.1, h3.2]

end CCMDetectorChain

#print axioms CCMDetectorChain.drop_ge_of_steps
#print axioms CCMDetectorChain.strictAnti_of_steps
#print axioms CCMDetectorChain.quiet_robust
