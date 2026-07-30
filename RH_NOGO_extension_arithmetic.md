# A structural no-go: self-adjoint-extension freedom is arithmetic-free

**Session 2026-07-30, generalizing §466/§467.** The CCM-specific NULL (no self-adjoint extension of the
archimedean prolate operator reproduces S(T)) is a special case of a clean, general obstruction to one
flavor of the Hilbert–Pólya program. This note states the general principle and its scope honestly.

## Proposition (extension freedom is a bounded, rank-≤n counting perturbation)
Let D be a densely-defined symmetric operator with **finite** deficiency indices (n,n), n<∞, so its
self-adjoint extensions form a U(n) manifold. Let W(λ) be the common leading Weyl term (extension-invariant:
all extensions share the same essential/absolutely-continuous or asymptotic-density structure; only the
discrete count near the boundary of the spectrum moves). Define each extension's counting fluctuation
F_θ(λ) = N_{A_θ}(λ) − W(λ). Then:

1. **(Krein / finite-rank resolvent.)** For any two self-adjoint extensions A_θ, A_θ', the resolvent
   difference (A_θ − z)^{-1} − (A_θ' − z)^{-1} has rank ≤ n.
2. **(Counting bound.)** Hence |N_{A_θ}(λ) − N_{A_θ'}(λ)| ≤ n for all λ, so
   ‖F_θ − F_θ'‖_∞ ≤ n and diam_{L∞}{ F_θ : θ ∈ U(n) } ≤ 2n.

**Corollary (no arithmetic from the boundary).** Suppose there is a *target* fluctuation Φ(λ) that is
**unbounded** (e.g. Φ = S(T), with Var S = (1/2π²)loglog T → ∞, S = Ω_±((log T/loglog T)^{1/3})). Then:
```
   some extension has  F_θ*(λ) = Φ(λ) + o(Φ)     ⟺     the bare/reference F_0 already does,
```
because all F_θ agree to within the bounded constant n, while Φ → ∞. **The self-adjoint-extension degree of
freedom cannot CREATE an unbounded arithmetic fluctuation that is absent from the bare operator.** Arithmetic
content, if present at all, is an invariant of the *bulk* operator (its symbol / potential / spectral
measure), not of the boundary condition. Boundary tuning moves at most n levels; it cannot manufacture the
prime-driven oscillation S(T).

## Why this bites the "xp / archimedean prolate" Hilbert–Pólya dream
The Berry–Keating / Connes–Consani–Moscovici archimedean prolate operator is the natural "geometric" candidate
whose *mean* density already equals the Riemann–von Mangoldt term ⟨N⟩ (Arb-certified exactly, §425). The hope
is that the *right* self-adjoint extension pins the eigenvalues onto the actual zeros — i.e. that the
extension supplies S(T). The Proposition says it categorically cannot:
- The archimedean operator is **1-D Sturm–Liouville**, so its bare fluctuation F_0 = O(1) (Prüfer
  rotation-number remainder is bounded; §467b). It carries **no** prime content — it is arithmetic-free.
- Its deficiency is **(2,2)** (§424), so extension freedom is a bounded ≤2 counting shift (§467).
- S(T) is unbounded. A bounded F_0 plus a ≤2 shift stays bounded; it cannot track an unbounded S(T).
Hence **no self-adjoint extension** of the archimedean operator reproduces S(T) — confirmed empirically to
r≈0.1 across the whole U(2) family and 3 independent instruments (§466/§467/§467d), and now explained: the
arithmetic was never a boundary degree of freedom to begin with.

## The constructive consequence (where the primes must live)
S(T) = −(1/π) Σ_p Σ_{k≥1} p^{−k/2} sin(kT log p)/k  + smooth  — its content IS the primes (the Euler product
of ζ). By the Proposition this cannot enter through a boundary condition on an arithmetic-free operator; it
must enter through the **bulk**. This is exactly the CCM **semilocal** construction (§465, CCM 2023): replace
the scaling operator by the Jacobi matrix of the measure dm_S(s) = |∏_{v∈S} L_v(½−is)|² ds, so each prime p
adds an |L_p|² Euler factor **to the measure** (the bulk), not the boundary. The delicate self-adjoint
extension that matters is the one of the *semilocal* operator — where the primes are already inside — and
its self-adjointness ⟺ semilocal Weil positivity ⟺ (in the limit over all places) RH. The §451–457
irreducible off-diagonal coupling is the analytic-side image of the same bulk-arithmetic requirement.

## Scope and honesty
- **Rigorous, general:** parts 1–2 (rank-≤n ⇒ ‖F_θ−F_θ'‖_∞ ≤ n) hold for any finite-deficiency symmetric
  operator. The Corollary's "can't create unbounded Φ from the boundary" is an immediate consequence.
- **Rigorous, CCM-specific:** deficiency (2,2) (§424, Arb) and Weyl term = ⟨N⟩ (§425, Arb); F_0 = O(1) via
  1-D Sturm–Liouville (§467b, modulo a certified uniform remainder bound for the singular endpoints).
- **Empirical:** that the bare F_0 is the O(1)/uncorrelated horn (r≈0.1, three instruments, X- and
  height-stable; §466/§467c/§467d).
- **This is a NEGATIVE structural result, not a proof of RH.** It rules out a *class* of approaches (supplying
  the arithmetic via a self-adjoint extension of an arithmetic-free archimedean operator) and redirects to
  the semilocal/bulk route. It neither proves nor disproves RH, and claims no Handle for it.

## Next (constructive test, loaded for the next build)
Directly test bulk-injection: does adding the |L_p|² factor to the measure (S = {p}, one prime) perturb the
prolate level positions with a **log-p-periodic** modulation — the prime-p contribution S_p(T) = −(1/π)Σ_k
p^{−k/2} sin(kT log p)/k that the archimedean operator entirely lacks? First-order measure perturbation of
the base prolate eigenfunctions is the tractable entry (before the full delicate semilocal extension).
Prerequisite (search-before-build): re-read CCM 2023 (arXiv:2310.18423) for the exact L²(ds) semilocal
operator definition and what its finite-S spectrum represents, so the perturbation is set up correctly
(avoid a §462-style false start).
