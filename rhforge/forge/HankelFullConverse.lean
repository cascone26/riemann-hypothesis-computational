import HankelConverseCore
import RealLinearSurj
import ConjClosedProdReal
import ProdEvalNeZero
import RealPolyCoeffVector

/-!
# §1803 — Full Hermite converse: a non-real root forces the Hankel form negative (CAPSTONE)

Assembles the general-case pieces into the complete converse.  Given a non-real `r₀`, a conjugation-closed
"rest" of the other roots with `r₀, r̄₀ ∉ rest`, build the isolating real polynomial `P = Q·L`:
`Q = ∏_{s∈rest}(X−s)` (real by §1800, nonzero at `r₀` by §1801), `L` the real linear factor from §1798 with
`L(r₀) = i/Q(r₀)`.  Then `P` is real, `P(r₀)=i`, and `P` vanishes on `rest`; extracting its real coefficient
vector (§1802) and feeding §1797 gives complex Hankel form `= -2 < 0`.  Hence a conj-closed root multiset with
ANY non-real element has a NON-PSD Hankel form — the Hermite criterion's hard direction, unconditional.
EXIT0 / std-axioms.
-/

open Polynomial Finset

/-- **Full Hermite converse.**  A conjugation-closed root multiset containing a non-real element has a
non-positive-semidefinite complex Hankel form (a real vector drives it to `-2`). -/
theorem hankel_not_psd_of_nonreal (rest : Multiset ℂ) (r₀ : ℂ)
    (hr : r₀.im ≠ 0) (hcc : rest.map (starRingEnd ℂ) = rest)
    (hr0 : r₀ ∉ rest) (hr0bar : (starRingEnd ℂ r₀) ∉ rest) :
    ∃ (n : ℕ) (c : Fin (n + 1) → ℝ),
      ((r₀ ::ₘ (starRingEnd ℂ r₀) ::ₘ rest).map
        (fun r => (∑ i : Fin (n + 1), (c i : ℂ) * r ^ (i : ℕ)) ^ 2)).sum = -2 := by
  set Q : Polynomial ℂ := (rest.map (fun s => X - C s)).prod with hQdef
  have hQreal : Polynomial.map (starRingEnd ℂ) Q = Q := conj_closed_prod_real rest hcc
  have hQr0 : Q.eval r₀ ≠ 0 := prod_eval_ne_zero rest r₀ hr0
  obtain ⟨a, b, hab⟩ := real_linear_surj_at_nonreal r₀ hr (Complex.I / Q.eval r₀)
  set L : Polynomial ℂ := C (a : ℂ) + C (b : ℂ) * X with hLdef
  have hLreal : Polynomial.map (starRingEnd ℂ) L = L := by
    simp [hLdef, Polynomial.map_add, Polynomial.map_mul, Polynomial.map_C, Polynomial.map_X,
          Complex.conj_ofReal]
  have hLeval : ∀ r, L.eval r = (a : ℂ) + (b : ℂ) * r := by intro r; simp [hLdef]
  set P : Polynomial ℂ := Q * L with hPdef
  have hPreal : Polynomial.map (starRingEnd ℂ) P = P := by
    rw [hPdef, Polynomial.map_mul, hQreal, hLreal]
  have hPr0 : P.eval r₀ = Complex.I := by
    rw [hPdef, Polynomial.eval_mul, hLeval, hab]
    field_simp
  have hPrest : ∀ s ∈ rest, P.eval s = 0 := by
    intro s hs
    rw [hPdef, Polynomial.eval_mul]
    have hQs : Q.eval s = 0 := by
      rw [hQdef]
      simp only [Polynomial.eval_multiset_prod, Multiset.map_map]
      apply Multiset.prod_eq_zero
      rw [Multiset.mem_map]
      exact ⟨s, hs, by simp⟩
    rw [hQs, zero_mul]
  obtain ⟨c, hc⟩ := real_poly_coeff_vector P hPreal
  refine ⟨P.natDegree, c, ?_⟩
  apply hankel_isolated_nonreal_neg rest c r₀
  · rw [hc r₀]; exact hPr0
  · intro s hs; rw [hc s]; exact hPrest s hs
