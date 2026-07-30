# Task 3 — What S(T)-shape can a self-adjoint extension inject? (structural companion to §466)

**Roundtable handoff Task 3, 2026-07-30.** Turns the empirical NULL of §466 (no U(2) extension of the
CCM archimedean prolate operator correlates with arithmetic S(T)) into a *structural* statement about why.

## Setup
- D = CCM archimedean prolate operator (λ=√2 build), deficiency indices **(2,2)** (§424, Arb-certified) —
  both singular endpoints x=√2 and x=∞ limit-circle. So D is not essentially self-adjoint; its self-adjoint
  extensions form a 4-real-parameter **U(2)** manifold.
- N_θ(T) = counting function of extension A_θ; ⟨N⟩(T) = (T/2π)log(T/2πe)+7/8 = the RvM main term, which is
  **θ-invariant and exactly correct** (§425, Arb-certified).
- F_θ(T) = N_θ(T) − ⟨N⟩(T)  (the operator's counting *fluctuation*).
- S(T) = (1/π) arg ζ(½+iT)  (the arithmetic fluctuation; N_ζ(T) = ⟨N⟩(T) + S(T)).

## The rigorous core (unconditional)
**Krein / finite-rank resolvent bound.** Any two self-adjoint extensions of a symmetric operator with
deficiency indices (n,n) have resolvent difference of rank ≤ n. For deficiency (2,2), n=2:
```
        (A_θ − z)^{-1} − (A_θ' − z)^{-1}   has rank ≤ 2.
```
A rank-≤2 resolvent perturbation moves the eigenvalue counting function by at most 2 (eigenvalue
interlacing / Weyl inequalities):
```
        | N_θ(T) − N_θ'(T) | ≤ 2      for ALL T, ALL θ, θ' ∈ U(2).
```
Since ⟨N⟩ is θ-invariant, F_θ − F_θ' = N_θ − N_θ', hence
```
   ‖ F_θ − F_θ' ‖_∞  ≤ 2     ⇒     diam_{L∞} { F_θ : θ ∈ U(2) }  ≤ 4.
```
**The entire extension family of fluctuations lives in a uniform sup-norm tube of half-width 2 around the
bare fluctuation F_0.** This is exactly what the extension knob "can inject": a *uniformly bounded* (≤2
levels, anywhere, forever) additive correction to the counting function. It can locally slide the level
ladder, but it cannot change the fluctuation by more than 2 at any height.

## The arithmetic target is unbounded
S(T) is **unbounded**:
- Selberg: Var S(T) = (1/2π²) log log T · (1+o(1)) → ∞.
- Montgomery (1977): S(T) = Ω_± ( (log T)^{1/3} (log log T)^{-1/3} ) — unconditional Ω result; it exceeds
  any constant infinitely often.

## The consequence (the Hilbert–Pólya gap, as a near-theorem)
Because every F_θ differs from the bare F_0 by ≤2 uniformly, and S(T) → ∞:
```
   F_θ(T) ~ S(T)  as T→∞   ⟺   F_0(T) ~ S(T)   (they differ by a bounded ≤2 term against an unbounded S).
```
So the boundary condition is powerless to *create* the arithmetic fluctuation: **either the operator's
intrinsic (bare) spectrum already encodes S(T) up to a bounded ±2 error, or no self-adjoint extension
does.** The extension freedom is categorically too weak — a bounded ≤2 additive shift — to supply an
unbounded, prime-driven S(T).

The measurement settles which horn: §426 (δ_n = O(1), sub-mean-spacing, bounded) + §427/§466/Task 1
(bare F_0 uncorrelated with S, r ≈ 0.10, X-stable) show the bare F_0 is **O(1)-bounded and not S(T)-shaped**.
A bounded F_0 cannot be within ±2 of an unbounded S(T) for large T. Therefore **no self-adjoint extension of
the archimedean CCM prolate operator reproduces S(T)** — the empirical U(2) NULL (§466) is the finite-height
shadow of this bounded-vs-unbounded obstruction.

## Where S(T) must come from instead
S(T)'s content is the primes (arg ζ = Σ over the Euler product). The extension (a boundary phenomenon)
cannot carry it; it must enter through the **operator's bulk** — the CCM **semilocal** construction (§465,
CCM 2023): replace the scaling operator S by the Jacobi matrix of the measure dm_S(s)=|∏_{v∈S}L_v(½−is)|²ds,
so each prime p adds an |L_p|² Euler factor to the measure. The arithmetic fluctuation is injected by the
measure (the primes), not by a self-adjoint extension of the pure-archimedean operator. This is precisely
the delicate self-adjoint extension of the *semilocal* operator that CCM's own program targets, and the
§451–457 irreducible off-diagonal coupling from the analytic side.

## Status / honesty
- **Rigorous, unconditional:** the ‖F_θ − F_θ'‖_∞ ≤ 2 bound (deficiency (2,2) ⇒ rank-≤2 resolvent
  difference ⇒ counting functions differ by ≤2) and the unboundedness of S(T) (Selberg/Montgomery).
- **Report (empirical):** that the bare F_0 is the O(1)/uncorrelated horn (§426/§427/§466/Task 1, scipy
  floats). Arb-certification of the bare NULL (prov/arb_ST.py) would harden this horn.
- **No Handle for RH.** This is a sharp *negative* structural result: it locates the obstruction (arithmetic
  content is not a boundary/extension degree of freedom of the archimedean operator) and points to the
  semilocal measure as the only route — it does not prove or disprove RH.
