import DerivBoundBCGen
import NegLeReOfNormLe

/-!
# §1318 — Borel–Carathéodory lower bound on `Re f'(0)`: `Re(f'(0)) ≥ -(4M+6‖f 0‖)/R`

Composing the forge's Borel–Carathéodory derivative bound §1237 (`‖f'(0)‖ ≤ (4M+6‖f 0‖)/R` for `f` analytic on
`ball 0 R` with `Re f ≤ M`) with §1317 (`‖g‖≤B ⟹ -B≤Re g`).  Applied to `f = log(ζ/P)` (holomorphic branch of the
log of the nonvanishing quotient), `f'(0) = logDeriv(ζ/P)` at the disk centre, so this is exactly the estimate
`Re(logDeriv(ζ/P)) ≥ -B` that §1314 needs (with `B = (4M+6‖f 0‖)/R ≍ log|t|`).  EXIT0 / std-axioms.
-/

open Complex

theorem re_deriv_ge_of_re_le {f : ℂ → ℂ} {R M : ℝ} (hR : 0 < R) (hM : 0 < M)
    (hf : DifferentiableOn ℂ f (Metric.ball 0 R))
    (hre : ∀ z ∈ Metric.ball (0 : ℂ) R, (f z).re ≤ M) :
    -((4 * M + 6 * ‖f 0‖) / R) ≤ (deriv f 0).re :=
  neg_le_re_of_norm_le (deriv_bound_of_re_le_general hR hM hf hre)

#print axioms re_deriv_ge_of_re_le
