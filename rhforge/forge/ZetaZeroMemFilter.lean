import ZetaZeroMemSupport

/-!
# §1407 — the DLVP zero lies in the **filtered** (half-disk) support (§1406's `hρ₀`)

§1406's `hρ₀` needs the kept zero `ρ = β+iγ` in `S = support.filter p` (the half-disk zeros).  From §1353 it is in
the divisor support, and if `p ρ` holds (i.e. `ρ ∈ closedBall (σK+iγ) (R/2)`, which is `σK-β ≤ R/2`) then it is in
the filter.  EXIT0 / std-axioms.
-/

open Metric Classical

theorem zeta_zero_mem_filter {c : ℂ} {R' : ℝ} (hball : (1 : ℂ) ∉ ball c R') {ρ z₀ : ℂ}
    (hρ : ρ ∈ ball c R') (hζρ : riemannZeta ρ = 0) (hz₀ : z₀ ∈ ball c R')
    (hz₀ne : riemannZeta z₀ ≠ 0) (p : ℂ → Prop) (hp : p ρ) :
    ρ ∈ (divisor_riemannZeta_ball_support_finite c R').toFinset.filter p :=
  Finset.mem_filter.mpr ⟨zeta_zero_mem_divisor_support hball hρ hζρ hz₀ hz₀ne, hp⟩

#print axioms zeta_zero_mem_filter
