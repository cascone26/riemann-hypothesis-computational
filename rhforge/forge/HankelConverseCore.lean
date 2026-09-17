import ComplexHankel

/-!
# §1797 — Hermite converse core: an isolated non-real root forces the Hankel form negative

The logical heart of the converse (PSD ⟹ all roots real).  Suppose a real coefficient vector `c` *isolates* a
non-real root `r₀` of the (conjugation-closed) root multiset `r₀ ::ₘ conj r₀ ::ₘ rest`: namely `P_c(r₀) = i`
and `P_c(s) = 0` for every `s ∈ rest`.  Because `c` is real, `P_c(conj r₀) = conj(P_c r₀) = -i` automatically,
so the complex Hankel form (= `∑_r P_c(r)²`, §1795) evaluates to `i² + (-i)² + 0 = -2 < 0`.  Hence the form is
NOT positive semidefinite — the contradiction that, once the isolating `c` is produced by interpolation, proves
"PSD ⟹ no non-real roots".  EXIT0 / std-axioms.
-/

open Finset

variable {d : ℕ}

/-- **Converse core.**  A real vector isolating a non-real root drives the complex Hankel form to `-2`. -/
theorem hankel_isolated_nonreal_neg (rest : Multiset ℂ) (c : Fin (d + 1) → ℝ) (r₀ : ℂ)
    (hr₀ : (∑ i : Fin (d + 1), (c i : ℂ) * r₀ ^ (i : ℕ)) = Complex.I)
    (hrest : ∀ s ∈ rest, (∑ i : Fin (d + 1), (c i : ℂ) * s ^ (i : ℕ)) = 0) :
    ((r₀ ::ₘ (starRingEnd ℂ r₀) ::ₘ rest).map
        (fun r => (∑ i : Fin (d + 1), (c i : ℂ) * r ^ (i : ℕ)) ^ 2)).sum = -2 := by
  have hconjP : (∑ i : Fin (d + 1), (c i : ℂ) * (starRingEnd ℂ r₀) ^ (i : ℕ))
      = starRingEnd ℂ (∑ i : Fin (d + 1), (c i : ℂ) * r₀ ^ (i : ℕ)) := by
    rw [map_sum]
    refine Finset.sum_congr rfl (fun i _ => ?_)
    rw [map_mul, map_pow, Complex.conj_ofReal]
  have hzero : rest.map (fun r => (∑ i : Fin (d + 1), (c i : ℂ) * r ^ (i : ℕ)) ^ 2)
      = rest.map (fun _ => (0 : ℂ)) := by
    refine Multiset.map_congr rfl (fun s hs => ?_)
    rw [hrest s hs]; ring
  simp only [Multiset.map_cons, Multiset.sum_cons]
  rw [hr₀, hconjP, hr₀, Complex.conj_I, hzero]
  simp only [Multiset.map_const', Multiset.sum_replicate, smul_zero, neg_sq, Complex.I_sq]
  ring
