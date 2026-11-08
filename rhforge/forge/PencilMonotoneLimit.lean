/-
PencilMonotoneLimit.lean (HP RH center, LEDGER S2906, 2026-11-07)
S2884 COROLLARY in Lean: monotone limits. Let e M k be the k-th LARGEST squared pencil eigenvalue at level M
(decreasing order, junk outside range). The capstone S2884.pencil_interlacing gives e (M+1) (k+1) ≤ e M k.
Reindex by "j-th SMALLEST": g j n := e (n + j + 1) n. Then each g j is antitone in the level, and if the values
are ≥ 0 (positivity of the pencil, S ≻ 0) it CONVERGES:  θ²_j(∞) := lim g j exists, and θ²_j(∞) ≤ θ²_j(N) ∀ N.
-/
import Mathlib

open Filter Topology

namespace PencilLimit

/-- the j-th smallest eigenvalue along levels -/
def g (e : ℕ → ℕ → ℝ) (j n : ℕ) : ℝ := e (n + j + 1) n

theorem g_antitone (e : ℕ → ℕ → ℝ) (h : ∀ M k, k < M → e (M + 1) (k + 1) ≤ e M k) (j : ℕ) :
    Antitone (g e j) := by
  apply antitone_nat_of_succ_le
  intro n
  have := h (n + j + 1) n (by omega)
  simp only [g]
  rwa [show n + j + 1 + 1 = n + 1 + j + 1 by omega] at this

theorem g_tendsto (e : ℕ → ℕ → ℝ) (h : ∀ M k, k < M → e (M + 1) (k + 1) ≤ e M k)
    (hpos : ∀ M k, 0 ≤ e M k) (j : ℕ) :
    ∃ l : ℝ, Tendsto (g e j) atTop (𝓝 l) ∧ 0 ≤ l ∧ ∀ n, l ≤ g e j n := by
  have hanti := g_antitone e h j
  have hbdd : BddBelow (Set.range (g e j)) := ⟨0, by rintro _ ⟨n, rfl⟩; exact hpos _ _⟩
  refine ⟨⨅ n, g e j n, tendsto_atTop_ciInf hanti hbdd, ?_, fun n => ciInf_le hbdd n⟩
  exact le_ciInf fun n => hpos _ _

end PencilLimit

#print axioms PencilLimit.g_antitone
#print axioms PencilLimit.g_tendsto
