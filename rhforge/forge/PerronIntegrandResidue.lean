import NegLogDerivZetaResidue
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv

set_option maxHeartbeats 1000000

/-!
# §1608 — the residue of the Perron integrand at `s=1` is `x²/2` (the PNT MAIN TERM)

`F(s) = -ζ'/ζ(s)·x^{s+1}/(s(s+1))` has a simple pole at `s=1` with residue exactly `x²/2`:
`(s-1)·F(s) = [(s-1)·(-ζ'/ζ(s))]·[x^{s+1}/(s(s+1))] → 1 · x²/2 = x²/2`, combining the `-ζ'/ζ` residue `1` (§1603)
with the continuity of the smoothed kernel `x^{s+1}/(s(s+1))` at `s=1` (value `x²/2`).  This is the numerical value
the contour shift deposits as the main term: `Ψ(x) = x²/2 + O(error)`.  EXIT0 / std-axioms.
-/

open Complex Filter Topology

theorem perron_integrand_residue_one {x : ℝ} (hx : 0 < x) :
    Tendsto (fun s : ℂ => (s - 1) *
        (-logDeriv riemannZeta s * (x : ℂ) ^ (s + 1) / (s * (s + 1))))
      (𝓝[≠] (1 : ℂ)) (𝓝 ((x : ℂ) ^ (2 : ℂ) / 2)) := by
  have hxne : (x : ℂ) ≠ 0 := by exact_mod_cast ne_of_gt hx
  -- residue of `-ζ'/ζ` at 1 is 1
  have hres := neg_logDeriv_zeta_residue_one
  -- the kernel `x^{s+1}/(s(s+1))` is continuous at 1, value `x^{2}/2`
  have hden1 : (1 : ℂ) * ((1 : ℂ) + 1) ≠ 0 := by norm_num
  have hKcont : ContinuousAt (fun s : ℂ => (x : ℂ) ^ (s + 1) / (s * (s + 1))) 1 := by
    apply ContinuousAt.div
    · exact (((differentiableAt_id (𝕜 := ℂ)).add_const 1).const_cpow (Or.inl hxne)).continuousAt
    · exact ((differentiableAt_id (𝕜 := ℂ)).mul
        ((differentiableAt_id (𝕜 := ℂ)).add_const 1)).continuousAt
    · exact hden1
  have hK : Tendsto (fun s : ℂ => (x : ℂ) ^ (s + 1) / (s * (s + 1)))
      (𝓝[≠] (1 : ℂ)) (𝓝 ((x : ℂ) ^ ((1 : ℂ) + 1) / ((1 : ℂ) * ((1 : ℂ) + 1)))) :=
    hKcont.continuousWithinAt
  have hmul := hres.mul hK
  -- rewrite the target function into the factored product form
  have hfun : (fun s : ℂ => (s - 1) *
        (-logDeriv riemannZeta s * (x : ℂ) ^ (s + 1) / (s * (s + 1))))
      = (fun s : ℂ => ((s - 1) * (-logDeriv riemannZeta s)) *
        ((x : ℂ) ^ (s + 1) / (s * (s + 1)))) := by
    funext s; ring
  rw [hfun]
  -- the limit value `1 · x^{1+1}/(1·(1+1)) = x^2/2`
  have h2 : (1 : ℂ) + 1 = 2 := by norm_num
  have hval : (1 : ℂ) * ((x : ℂ) ^ ((1 : ℂ) + 1) / ((1 : ℂ) * ((1 : ℂ) + 1)))
      = (x : ℂ) ^ (2 : ℂ) / 2 := by rw [h2]; ring
  rw [← hval]
  exact hmul

#print axioms perron_integrand_residue_one
