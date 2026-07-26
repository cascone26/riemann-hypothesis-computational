# SUPERPRISM 2026-07-24 (Fable, Sam) — the metric question resolved into the a-flow

Context: Lane B (§234–§240) compressed RH to "does the τ-Krein / scaling-isometric metric exist"
(AKM template). This brainstorm collided that with the repo's own Suzuki asset (§41–45, §59b) and
produced one new reformulation + two structural lenses. Grep-first discipline applied: the naive
"finite HB + Hurwitz" idea was §59b already (Suzuki Thm 1.5 = finite HB PROVEN; §60 anti-rediscovery
flag on Cor 1.6). What survives is new.

## 1. Hurwitz sharpening (why the θ(a)-limit is the hard part, precisely)
If e^φ W(a,θ(a);z) → z²ξ(½−iz)/ξ′(½−iz) locally uniformly on compacts, Hurwitz forces every ξ-zero
to be a limit of W's (proven-real) zeros ⟹ RH follows from convergence ALONE. So tuned locally-uniform
convergence is itself RH-hard — the route cannot be closed by soft analysis; it must WEAKEN the
convergence mode. Trace-level weakening collapses back to the Weil criterion (checked: spectral-measure
convergence against PW test class = explicit-formula positivity = Weil). The only mode carrying NEW
information is PER-BRANCH convergence — exactly what §45b measured (0.073→0.025) but never framed.

## 2. THE NEW FRAME — RH as gap-persistence of the Suzuki a-flow
Given branch convergence + strip-counting convergence (RvM N(T) counts the strip regardless of Re ρ),
an off-line zero pair can hide in the real spectra ONLY as two adjacent branches colliding IN THE LIMIT
at the pair's common ordinate:
    **RH + simple zeros ⟺ ∀ fixed k: liminf_{a→∞} gap_k(a) > 0**, gap_k = z_{k+1}(a) − z_k(a) of D_{a,θ(a)}.
Structure of the flow (all previously-proven pieces, newly assembled):
  - between thresholds a = ½log q: smooth scaling flow (no prime data changes);
  - at each threshold: shift-operator coupling S_{log q} enters CUBICALLY SOFTLY (= Onset Lemma §219,
    reappearing as the insertion events of this flow);
  - Lehmer pairs = near-misses of the branch dynamics (dynamical, not statistical);
  - the dBN t-flow is CLOSED (§47e); this is a DIFFERENT flow (in a, on the operator side, unrun).
Honest limits: finite-a simplicity is likely provable (interlacing) but insufficient — the content is
the LIMINF. Per-k (not uniform-in-k) is what's needed, much weaker than level repulsion. The danger
mechanism is resonant accumulation of insertion couplings at a pair frequency — the ζ′/ζ circularity
shadow appears HERE as a dynamical resonance question, which is a genuinely different handle.

## 3. One geometry under both lanes (arc lens)
Prime q's tail-basis kernel has Toeplitz symbol EXACTLY (π/2)·𝟙_arc, arc measure 2b·log q (verified
5 digits, banking lane). Criticality a = ½log q ⟺ q's arc reaches half-measure π ⟺ q's shift enters
the a-flow. Banking-lane wall positions = a-flow insertion events = arc half-wraps: three languages,
one threshold set. As a→∞ the arc-interference sum is Σ Λ(q)q^{−1/2+imb} = ζ′/ζ on the line — the
circularity, geometrically visible. Banking certifies the flow's positivity windows; RH needs its gap
dynamics.

## Experiment fired tonight
operators/suzuki_gap_flow.py — branch trajectories z_k(a) + gap_k(a), a ∈ [0.30, 1.01] dense with
refinement at thresholds {0.347, 0.549, 0.693, 0.805, 0.973}, θ ∈ {0, π/2, π, 3π/2}. Measures:
(i) branch smoothness across insertions (Onset-softness check on the OPERATOR side), (ii) min-gap
collision detector, (iii) drift of gaps from momentum-lattice regularity toward the irregular zeta
gaps (6.89, 3.99, 5.41, 2.51, ...). Next phase if the frame survives contact: Hellmann–Feynman
per-prime decomposition of d(gap_k)/da — the "resonance ledger" of which primes move which gaps.
Results → LEDGER §241.

## OUTCOME (same night): experiment returned a NEW measured phenomenon — see LEDGER §241
Resonance staircase, θ-ROBUST: branches of D_{a,θ} transit fast (up to 5× free-scaling) between zeros and
PARK at zeta ordinates (velocity ~0, ≥10× sub-free) — plateau positions agree across all four extensions
θ∈{0,π/2,π,3π/2} and sit at γ₁..γ₆ with offsets ≤0.1 shrinking in a. Handoffs are avoided crossings; the
Cor 1.6 tuned θ(a) = choosing WHICH branch is parked; the parking is universal. New RH observable: plateau
minimum velocity v_k(a) — RH ⟺ v_k(a)→0 ∀k (an off-line zero of width δ floors its plateau velocity at ~δ,
Fano/finite-size-scaling logic). Next: v_k(a) scaling to larger a (needs §42 singularity-subtraction
discretization), then per-prime Hellmann–Feynman decomposition of the parking force.
