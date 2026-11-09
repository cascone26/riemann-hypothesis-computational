import CCMWeightShape
import CCMResonance

/-!
# S2946 — two-regime Loewner envelope for the CCM zero side (lane 2)

From S2943 `ccmW_le_min` (w_γ ≤ min(4/L, Lγ²)) and S2937 Loewner monotonicity:
* `ccm_gram_le_envelope` : G_ccm ⪯ C·diag(min(4/L, Lγ_j²))·Cᵀ;
* `ccm_gram_le_quad` : G_ccm ⪯ L·C·diag(γ_j²)·Cᵀ  (low-lying zeros contribute O(Lγ²));
* trace/Rayleigh forms of both.
The envelope is the sharp combination of the high-zero (4/L) and low-zero (Lγ²) brackets.
-/

open Matrix Finset

namespace CCMEnvelope

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq n] [DecidableEq m]

theorem ccm_gram_le_envelope {L : ℝ} (hL : 0 < L) (C : Matrix n m ℝ) (γ : m → ℝ) :
    (C * diagonal (fun j => min (4 / L) (L * γ j ^ 2)) * Cᵀ
      - C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ).PosSemidef :=
  CCMLoewner.gram_mono C (fun j => CCMWeight.ccmW_le_min hL (γ j))

theorem ccm_gram_le_quad {L : ℝ} (hL : 0 < L) (C : Matrix n m ℝ) (γ : m → ℝ) :
    (C * diagonal (fun j => L * γ j ^ 2) * Cᵀ
      - C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ).PosSemidef :=
  CCMLoewner.gram_mono C (fun j => CCMWeight.ccmW_le_sq hL (γ j))

theorem ccm_rayleigh_le_envelope {L : ℝ} (hL : 0 < L) (C : Matrix n m ℝ) (γ : m → ℝ) (x : n → ℝ) :
    x ⬝ᵥ ((C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ) *ᵥ x)
      ≤ x ⬝ᵥ ((C * diagonal (fun j => min (4 / L) (L * γ j ^ 2)) * Cᵀ) *ᵥ x) :=
  CCMResonance.rayleigh_mono C (fun j => CCMWeight.ccmW_le_min hL (γ j)) x

theorem ccm_trace_le_envelope {L : ℝ} (hL : 0 < L) (C : Matrix n m ℝ) (γ : m → ℝ) :
    (C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ).trace
      ≤ (C * diagonal (fun j => min (4 / L) (L * γ j ^ 2)) * Cᵀ).trace :=
  CCMResonance.trace_mono C (fun j => CCMWeight.ccmW_le_min hL (γ j))

theorem ccm_trace_le_quad {L : ℝ} (hL : 0 < L) (C : Matrix n m ℝ) (γ : m → ℝ) :
    (C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ).trace
      ≤ (C * diagonal (fun j => L * γ j ^ 2) * Cᵀ).trace :=
  CCMResonance.trace_mono C (fun j => CCMWeight.ccmW_le_sq hL (γ j))

end CCMEnvelope

#print axioms CCMEnvelope.ccm_gram_le_envelope
#print axioms CCMEnvelope.ccm_gram_le_quad
#print axioms CCMEnvelope.ccm_rayleigh_le_envelope
#print axioms CCMEnvelope.ccm_trace_le_envelope
#print axioms CCMEnvelope.ccm_trace_le_quad
