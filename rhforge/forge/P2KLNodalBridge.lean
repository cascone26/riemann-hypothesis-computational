import QOffKRationalIdentity

open scoped MatrixOrder Matrix
open Polynomial

/-!
# Off-critical-line lane: `P_{2K,L}` IS the Lagrange nodal polynomial (bridging notations)

Connects §2366-2367's hand-built `P_2KL` (a plain real-valued function, matching the paper's own
Eq (7.10) notation `P_{2K,L}(w):=∏_n(w−ν_{n,L})`) to the abstract `Lagrange.nodal` object that
§2389's `Q_off_K_rational_identity` states its conclusion in terms of — letting future work restate
that identity using the paper's own `P_{2K,L}` notation directly, rather than the unfamiliar
`Lagrange.nodal`. The two are literally the same product, differing only in which "package" (a bare
function vs. a genuine `Polynomial` term via `Lagrange.eval_nodal`) computes it. Search-before-
build: grepped for `P_2KL_eq_nodal`/`nodal_eq_P_2KL` — fresh. Fully general, no zeta/RH dependence.
-/
theorem P_2KL_eq_nodal_eval (K : ℕ) (L : ℝ) (w : ℝ) :
    ((P_2KL K L w : ℝ) : ℂ)
      = (Lagrange.nodal (Finset.univ : Finset (Fin (2 * (2 * K) + 1)))
            (fun n : Fin (2 * (2 * K) + 1) =>
              (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))).eval (w : ℂ) := by
  unfold P_2KL
  rw [Lagrange.eval_nodal]
  push_cast
  apply Finset.prod_congr rfl
  intro n _
  push_cast
  ring

#print axioms P_2KL_eq_nodal_eval
