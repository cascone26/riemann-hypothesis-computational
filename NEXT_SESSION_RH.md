# Next RH session — the analytic coercivity target (start here)

This session (2026-07-19/20) mapped the Weil operator τ exhaustively and produced two unconditional theorems
(RH_RESULTS.md). Do NOT re-map τ — every angle (margin scaling, near-null modes, √τ, Loewner, FE sectors, the
DBN toy) lands on the same sign-indefinite core |ω|<6.3. Read FIRST: LEDGER.md §10–§34d, OBSTRUCTION_MAP.md,
memories `connes-decisive-diagnostic`, `rh-dbn-criticality-crack` (the DBN proof-route was CIRCULAR — do not
retry the toy τ+t·diag(r²)).

**The one real proof target** (unchanged, from Connes/Track-2 analysis): prove UNCONDITIONAL uniform coercivity
  λ_min(τ_N) ≥ 0  for ALL N (all bandwidths), i.e. A(g)+B(g) ≥ P(g) bandwidth-uniformly.
This IS RH (Bombieri/Weil). The content is ENTIRELY analytic — numerics can't close it (finitely many N).
Two non-circular handles only:
  (i) explicit coercivity / SOS from the digamma symbol Ω(r)=Re ψ(¼+ir/2)−log π's Hadamard structure + an
      unconditional prime-sum bound. Ω is sign-indefinite (Ω<0 for |r|<6.29) — the negative core is the whole
      fight (§12). The exact off-diagonal Pick form F_n=... and the closed form F_n=4−4cos(nθ)cosh(nu) are tools.
  (ii) Hermite–Biehler for ξ̂ not using any zero location.
**Existing-knowledge route to criticality:** RH=Λ=0; Rodgers–Tao proved Λ≥0, Polymath15 proved Λ≤0.22. The
real (not toy) de Bruijn–Newman flow is the object; pushing 0.22→0 is a heavy but genuine program.

Come in fresh, pick (i) or the Polymath15 route, and do HARD ANALYSIS, not more τ probing.
