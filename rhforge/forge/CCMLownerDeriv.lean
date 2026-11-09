import CCMTailMoment2
import CCMLownerDisplacement

/-!
# S2985 — Löwner entries are difference quotients: derivative control (lane 2)

S2983's honest reading: the CCM margin lives in the far tail's FINE structure; entrywise expansions converge slowly.
The exact tool: every CCM block is Löwner (S2968/S2970/S2972), with off-diagonal entries
`(g(ν_k) − g(ν_i))/(π(ν_i − ν_k))` = −(1/π)·(difference quotient of g).  So derivative bounds on the generator g
control the whole off-diagonal part at once (mean value theorem):
* `lowner_entry_le_of_deriv` — g with derivative g' on [−V, V], |g'| ≤ D there, ν_i ≠ ν_k in [−V, V] ⇒
  |(g ν_k − g ν_i)/(π(ν_i − ν_k))| ≤ D/π.
* `lowner_offdiag_le_of_deriv` — the same for every off-diagonal entry of `CCMLowner.lowner d (g ∘ ν) ν`.
Next: feed the exact tail generator g_tail(ν) = π Σ_j w_j·2ν/(ν² − x_j²) (S2974 numerics) — its k-th derivative is
O(Σ w/x^{k+1}) and the k-th divided differences give the moment hierarchy with sharp constants.
-/

open Set

namespace CCMLownerDeriv

theorem lowner_entry_le_of_deriv (g g' : ℝ → ℝ) {V D : ℝ}
    (hg : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g (g' x) (Icc (-V) V) x)
    (hD : ∀ x ∈ Icc (-V) V, |g' x| ≤ D) {a b : ℝ} (ha : a ∈ Icc (-V) V) (hb : b ∈ Icc (-V) V)
    (hab : a ≠ b) :
    |(g b - g a) / (Real.pi * (a - b))| ≤ D / Real.pi := by
  have hmv := (convex_Icc (-V) V).norm_image_sub_le_of_norm_hasDerivWithin_le hg
    (fun x hx => by simpa [Real.norm_eq_abs] using hD x hx) ha hb
  simp only [Real.norm_eq_abs] at hmv
  have hne : a - b ≠ 0 := sub_ne_zero.mpr hab
  have hpos : 0 < |a - b| := abs_pos.mpr hne
  have hpi := Real.pi_pos
  rw [abs_div, abs_mul, abs_of_pos hpi, div_le_div_iff₀ (by positivity) hpi]
  rw [abs_sub_comm b a] at hmv
  nlinarith [abs_nonneg (g b - g a)]

theorem lowner_offdiag_le_of_deriv {ι : Type*} [DecidableEq ι] (d : ι → ℝ) (ν : ι → ℝ)
    (g g' : ℝ → ℝ) {V D : ℝ}
    (hg : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g (g' x) (Icc (-V) V) x)
    (hD : ∀ x ∈ Icc (-V) V, |g' x| ≤ D) (hν : ∀ i, ν i ∈ Icc (-V) V)
    (hinj : Function.Injective ν) {i k : ι} (hik : i ≠ k) :
    |CCMLowner.lowner d (fun j => g (ν j)) ν i k| ≤ D / Real.pi := by
  simp only [CCMLowner.lowner, if_neg hik]
  exact lowner_entry_le_of_deriv g g' hg hD (hν i) (hν k) (fun e => hik (hinj e))

end CCMLownerDeriv

#print axioms CCMLownerDeriv.lowner_entry_le_of_deriv
#print axioms CCMLownerDeriv.lowner_offdiag_le_of_deriv
