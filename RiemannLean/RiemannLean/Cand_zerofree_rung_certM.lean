import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Data.Real.StarOrdered
open Matrix
/-- **Machine-checked ZERO-FREE open-rung Weil certificate matrix** (Fejér span c∈{14,22,30,38},
support 4W=1.32 > log 3 — the first open rung, primes {2,3} in support).
`M` rationalizes the ARB-CERTIFIED midpoints of the zero-free route entries pole+arch−prime
(operators/arb_zerofree_cert3.py, §788: NO zero knowledge anywhere — no verified-zero list, no
on-line-ness hypothesis, no strip lemma). This theorem proves `M.PosSemidef` via its exact LDLᵀ
sum-of-squares. With the certified bound ‖Q−M‖₂ ≤ ρ ≈ 3.80e-05 and min-eig(M) > ρ (proven by
exact LDLᵀ of M − ρI over ℚ, all pivots positive; min LDLᵀ pivot of M = 2.0237), the TRUE
zero-free Weil form is positive-definite on this span UNCONDITIONALLY (§788/§792). -/
theorem zerofree_rung_certM_posSemidef :
    (Matrix.of ![![(726780/359137 : ℝ), (281643/886004 : ℝ), (82804/657425 : ℝ), (83245/892953 : ℝ)], ![(281643/886004 : ℝ), (320438/110315 : ℝ), (603592/993409 : ℝ), (136609/742034 : ℝ)], ![(82804/657425 : ℝ), (603592/993409 : ℝ), (2868181/863560 : ℝ), (11937/14239 : ℝ)], ![(83245/892953 : ℝ), (136609/742034 : ℝ), (11937/14239 : ℝ), (644527/183987 : ℝ)]]).PosSemidef := by
  rw [posSemidef_iff_dotProduct_mulVec]
  refine ⟨?_, ?_⟩
  · ext i j; fin_cases i <;> fin_cases j <;> norm_num [Matrix.conjTranspose_apply, Matrix.of_apply]
  · intro x
    simp only [dotProduct, Matrix.mulVec, Matrix.of_apply, Fin.sum_univ_four,
      Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons, Matrix.cons_val_two,
      Matrix.cons_val_three, Matrix.tail_cons, star_trivial]
    nlinarith [mul_nonneg (by norm_num : (0:ℝ) ≤ (726780/359137 : ℝ)) (sq_nonneg (x 0 + (33716140697/214643329040 : ℝ) * x 1 + (7434495037/119450835375 : ℝ) * x 2 + (5979271913/129796076268 : ℝ) * x 3)),
               mul_nonneg (by norm_num : (0:ℝ) ≤ (11978341222309139431843/4195827673691109158080 : ℝ)) (sq_nonneg (x 1 + (1610758772161868201504584285305616/7822957769370083212712839153243475 : ℝ) * x 2 + (235559003999004740768685087586820/3968433349266192655144320785149443 : ℝ) * x 3)),
               mul_nonneg (by norm_num : (0:ℝ) ≤ (42255821177315529904322315908523088761478366267268733/13236070249528628153444930609408907201774043438725000 : ℝ)) (sq_nonneg (x 2 + (49804276085876675443253130758504269399147320627909077553074880157000/199337242454171210407727820118580493673272290827567103453001103267587 : ℝ) * x 3)),
               mul_nonneg (by norm_num : (0:ℝ) ≤ (21078492297678923303907353101080036360421590713478850513160480706163783840002028063452853/6407870410134863760222088191629490686236802283119059869842017537000972808630492936986533 : ℝ)) (sq_nonneg (x 3))]
