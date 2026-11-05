import Eq7p17DetIdentity
import DLNKOff0

open scoped MatrixOrder Matrix
open Polynomial

/-!
# Eq (7.17) — Polynomial Uniqueness Approach

Theorem 7.2's determinant identity via polynomial point-agreement:
- det_param(w) and det(G_off0)*Q_off_K(w) are both polynomials in w, degree ≤ 4*K
- If they agree at 4*K+1 distinct points, they are equal (Lagrange/Vandermonde)
- The frequency nodes ν_{n,L} provide natural evaluation points

KEY OBSERVATION:
- S_T_off_real is (4*K+1) × (4*K+1)
- K_off0 = Cᵀ * (S_T_off_real * D_LN) * C, where D_LN is diagonal with frequency nodes
- The null-vector fact S_T_off_real *ᵥ e_off = 0 constrains the structure
- Question: does this null-vector relationship imply that det_param vanishes at certain points?
-/

-- Step 1: Understand the S_T_off_real / det_param_matrix relationship
-- K_off0 is defined as: Cᵀ * (S_T_off_real * D_LN) * C
-- So det_param(w) = det(w·G_off0 - Cᵀ*(S_T_off_real*D_LN)*C)

-- The null-vector fact: S_T_off_real *ᵥ e_off = 0
-- This is proven in STOffFinalAssembly.lean as S_T_off_real_mulVec_e_off_eq_zero

-- Step 2: Can we relate the null-vector e_off to the evaluation of det_param at frequency nodes?
-- The vector e_off has components e_off(n) = √L * Q_off_K(ν_{n,L}) / P'_{2K,L}(ν_{n,L})
-- So e_off is defined precisely at the frequency nodes!

-- Observation: if e_off is a null vector of S_T_off_real, and e_off is constructed from
-- Q_off_K evaluated at the frequency nodes, maybe there's a connection to det_param?

-- Step 3: Verify the S_T_off_real ↔ K_off0 relationship is sound
-- Note: K_off0 unfolds to exactly this formula by definition
theorem K_off0_relation (K : ℕ) (L : ℝ) (z : Fin K → ℂ) :
    K_off0 (2 * K) L K z =
      let N := 2 * K
      (C_matrix (2 * N + 1) (by omega))ᵀ *
      (S_T_off_real N L K z * D_LN N L) *
      (C_matrix (2 * N + 1) (by omega)) := by
  unfold K_off0
  simp only []

-- Step 4: What we need to prove for polynomial uniqueness
-- We need: det_param(2*K, L, K, z, w) = det(G_off0(2*K, L, K, z)) * Q_off_K(K, z, w)
-- This is true for all w if it's true at 4*K+1 distinct points (by polynomial uniqueness)

-- HONEST ASSESSMENT:
-- The connection between the null-vector fact (S_T_off_real *ᵥ e_off = 0) and
-- point evaluations of det_param is NOT immediately clear from the definitions.
--
-- Possible interpretations:
-- 1. Maybe det_param vanishes at the frequency nodes (needs proof)
-- 2. Maybe det_param takes specific values at the frequency nodes (needs identification)
-- 3. Maybe the null-vector fact constrains det_param through the K_off0 structure
--
-- WITHOUT clarity on point (1), (2), or (3), the "natural candidate points" suggested
-- in the strategy cannot be formalized. The strategy as stated assumes there IS a
-- connection, but that connection is not yet explicitly proven in the corpus.

#print axioms K_off0_relation
