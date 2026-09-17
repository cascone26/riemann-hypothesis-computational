import ProdNormSubGe
import SphereSubHalfBall

/-!
# §1376 — `‖∏ (s-ρ)^{m ρ}‖ ≥ (R/2)^{Σ m ρ}` on the sphere (half-disk zeros)

Combining §1374 (half-disk zeros stay `≥ R/2` from the `R`-sphere) with §1375 (product lower bound): for `s` on
`sphere c R` and a finite set `S` of half-disk points (`S ⊆ closedBall c (R/2)`), the partial product
`P_half(s) = ∏_{ρ∈S} (s - ρ)^{m ρ}` obeys `‖P_half(s)‖ ≥ (R/2)^{Σ m ρ}`.  This is the `Plow` for §1373, giving the
sphere bound `‖ζ/P_half‖ ≤ Mζ/(R/2)^n` — the crux of the two-radii factor bound.  EXIT0 / std-axioms.
-/

open Complex

theorem norm_prod_sub_half_ge {S : Finset ℂ} {s c : ℂ} {R : ℝ} (hR : 0 ≤ R) (m : ℂ → ℕ)
    (hs : s ∈ Metric.sphere c R) (hS : ∀ ρ ∈ S, ρ ∈ Metric.closedBall c (R / 2)) :
    (R / 2) ^ (∑ ρ ∈ S, m ρ) ≤ ‖∏ ρ ∈ S, (s - ρ) ^ (m ρ)‖ := by
  rw [norm_prod]
  simp_rw [norm_pow]
  exact prod_norm_sub_ge (by linarith) m (fun ρ hρ => sphere_sub_halfball_ge hs (hS ρ hρ))

#print axioms norm_prod_sub_half_ge
