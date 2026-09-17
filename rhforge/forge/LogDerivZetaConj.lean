import Mathlib.NumberTheory.Harmonic.ZetaAsymp
import Mathlib.Analysis.Calculus.LogDeriv
import Mathlib.Analysis.Calculus.Deriv.Slope
import Mathlib.Tactic

set_option maxHeartbeats 1000000

open Complex Filter Topology

/-!
# §1690 — SCHWARZ-REFLECTION for `ζ'/ζ`: `‖logDeriv ζ (conj s)‖ = ‖logDeriv ζ s‖`

The bottom-connector (`Im = -T`) bound comes free from the top one by conjugation symmetry.
Built from scratch (Mathlib has `riemannZeta_conj` but NO derivative-conjugation lemma):
the reflection derivative `HasDerivAt (conj∘f∘conj) (conj d)` via slope conjugation, then
`deriv ζ (conj z) = conj (deriv ζ z)` (off `z = 1`), then the `logDeriv` and norm forms.
-/

/-- The Schwarz-reflection derivative: if `f` has derivative `d` at `conj z`, the reflected
function `conj ∘ f ∘ conj` has derivative `conj d` at `z`. -/
theorem HasDerivAt.conj_comp_conj {f : ℂ → ℂ} {z d : ℂ}
    (hf : HasDerivAt f d ((starRingEnd ℂ) z)) :
    HasDerivAt (fun s => (starRingEnd ℂ) (f ((starRingEnd ℂ) s))) ((starRingEnd ℂ) d) z := by
  rw [hasDerivAt_iff_tendsto_slope] at hf ⊢
  have hmap : Tendsto (starRingEnd ℂ) (𝓝[≠] z) (𝓝[≠] ((starRingEnd ℂ) z)) := by
    apply tendsto_nhdsWithin_of_tendsto_nhds_of_eventually_within
    · exact (Complex.continuous_conj.tendsto z).mono_left nhdsWithin_le_nhds
    · filter_upwards [self_mem_nhdsWithin] with s hs
      simp only [Set.mem_compl_iff, Set.mem_singleton_iff] at hs ⊢
      intro h
      exact hs ((starRingEnd ℂ).injective h)
  have hconj : Tendsto (starRingEnd ℂ) (𝓝 d) (𝓝 ((starRingEnd ℂ) d)) :=
    (Complex.continuous_conj.tendsto d)
  have hcomp := hconj.comp (hf.comp hmap)
  apply hcomp.congr'
  filter_upwards [self_mem_nhdsWithin] with s hs
  simp only [Set.mem_compl_iff, Set.mem_singleton_iff] at hs
  simp only [Function.comp_apply, slope_def_field]
  simp only [map_div₀, map_sub, Complex.conj_conj]

/-- `ζ' (conj z) = conj (ζ' z)` away from the pole. -/
theorem deriv_riemannZeta_conj {z : ℂ} (hz : z ≠ 1) :
    deriv riemannZeta ((starRingEnd ℂ) z) = (starRingEnd ℂ) (deriv riemannZeta z) := by
  have hζrefl : riemannZeta = fun s => (starRingEnd ℂ) (riemannZeta ((starRingEnd ℂ) s)) := by
    funext s
    rw [riemannZeta_conj, RingHomInvPair.comp_apply_eq]
  have hd : HasDerivAt riemannZeta (deriv riemannZeta z) z :=
    (differentiableAt_riemannZeta hz).hasDerivAt
  have hd' : HasDerivAt riemannZeta (deriv riemannZeta z)
      ((starRingEnd ℂ) ((starRingEnd ℂ) z)) := by
    rwa [RingHomInvPair.comp_apply_eq]
  have hrefl := hd'.conj_comp_conj
  rw [← hζrefl] at hrefl
  exact hrefl.deriv

/-- `logDeriv ζ (conj z) = conj (logDeriv ζ z)` away from the pole. -/
theorem logDeriv_riemannZeta_conj {z : ℂ} (hz : z ≠ 1) :
    logDeriv riemannZeta ((starRingEnd ℂ) z) = (starRingEnd ℂ) (logDeriv riemannZeta z) := by
  rw [logDeriv_apply, logDeriv_apply, deriv_riemannZeta_conj hz, riemannZeta_conj, ← map_div₀]

/-- The norm symmetry: the bottom edge inherits the top edge's bound. -/
theorem norm_logDeriv_riemannZeta_conj {z : ℂ} (hz : z ≠ 1) :
    ‖logDeriv riemannZeta ((starRingEnd ℂ) z)‖ = ‖logDeriv riemannZeta z‖ := by
  rw [logDeriv_riemannZeta_conj hz]
  exact RCLike.norm_conj _

#print axioms HasDerivAt.conj_comp_conj
#print axioms logDeriv_riemannZeta_conj
#print axioms norm_logDeriv_riemannZeta_conj
