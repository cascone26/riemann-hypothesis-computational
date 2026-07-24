# RH ↔ Physics/Universe — deep-dive, adversarially verified (2026-07-18)

Follow-up to `RH_PHYSICS_2026-07-16.md` (primon gas / Bost-Connes / GUE / Berry-Keating /
Connes-Consani). This round was a 104-agent deep-research fan-out (5 angles, 22 sources,
71 claims extracted, 25 verified by 3-vote adversarial refutation). Only 5 of 25 claims
survived verification — recording BOTH what survived and what got killed, because the kills
are the value (they stop us chasing overclaimed physics-proves-RH papers).

## CONFIRMED (survived adversarial verification)

1. **Derivative-moment RMT conjectures (Hughes & Pearce-Crump, arXiv:2509.07788, Sep 2025)** —
   NEW, beyond pair correlation. Predicts moments of ζ′ AT THE ZEROS: Σ_{0<γ≤T} ζ′(½+iγ)^k
   ~ (T/2π)(log T/2π)^{k+1}/Γ(k+2) for Re(k)>−3. Proven unconditionally only for k=1,−1;
   general/complex k conjectural. **Directly testable on our verified zeros** (needs ζ′(½+iγ_n),
   computable). This is the freshest, most concrete, genuinely-new test in the whole survey.

2. **Freezing phase transition (Curran, Q. J. Math. 2024)** — proven theorem: the 2nd
   moment-of-moments MoM_T(2,β) has a sharp transition at **β = 1/√2**, separating a
   variance-dominated regime (β<1/√2) from an extreme-value-dominated regime (β>1/√2). This
   is a real statistical-mechanics-style phase transition in ζ-moment data. Computable on a
   finite zero dataset. Physical interpretation OPEN (dynamical? or purely arithmetic?).

3. **Hyperuniform structure factor (Zhang, Martelli & Torquato, arXiv:1801.01541, 2018)** —
   primes AND zeta zeros are "hyperuniform": structure factor S(k) has Bragg-like peaks with
   a diffuse background that → 0 as system size → ∞. Already computed to ~20M zeros. This is
   the concrete realization of Dyson's quasicrystal idea. **CORRECTS our 07-18 assumption that
   the diffraction test was never run at scale — it was (Zhang-Torquato to 20M).**

4. **Dyson quasicrystal — partially done, specific gaps remain (medium confidence).** Structure
   factor consistency confirmed (#3), but NO paper (2018-2026) has: mapped full S(k) at
   ~2π/10⁷ resolution, fit peaks to a specific known quasicrystal lattice class, OR tied the
   diffuse-background DECAY RATE to anything about the zeros being on the line. That last one is
   the only genuinely-unexplored, possibly-non-circular sub-question.

5. **Selberg vs arithmetic gap (high confidence).** Selberg zeta on a hyperbolic surface has a
   genuine self-adjoint Laplacian → RH-analog PROVEN. Arithmetic ζ has no such operator on a
   real geometric space; Connes-Consani adele-class-space is the closest but stays abstract
   (the adele class space is not a manifold with a natural measure + self-adjoint operator in
   Selberg's sense). This IS the project's own Track 2 — no new escape hatch.

## KILLED (overclaims refuted 2-3 votes — do NOT chase these)

- **"Quasicrystal Fourier self-duality FORCES zeros onto the critical line" (Shaughnessy,
  arXiv:2410.03673, 2024)** — REFUTED 0-3. The paper's computational test (scattering spectrum
  of atoms placed at prime positions shows peaks at zero heights) is real, but its claim to
  PROVE RH via Fourier self-duality does not hold up. Same failure mode as every physics-proves-RH
  paper: a consistency signature dressed as a forcing mechanism.
- **Zhang-Torquato explicit peak-location k=mπ/n and peak-height formulas** — refuted/mixed. The
  hyperuniformity observation is solid; the specific closed-form peak formulas are not reliable.
- **Dynamical system with spectral density matching zeros (arXiv:2406.12852)** — has a real
  Lyapunov-instability problem (exponential growth, λ=ln(π²/6)), flagged in its own revision.
  Not viable as stated.
- Several Keating-Snaith moment-formula restatements got 1-2 votes — the underlying CUE/ζ moment
  conjecture (Keating-Snaith 2000) is standard and fine; the extractor's specific formula
  transcriptions were imprecise. Use the primary source, not the extracted claim.

## THE DEEP CONCEPTUAL FINDING — Lee-Yang is the wall, restated in physics

The Lee-Yang circle theorem is a case where "RH for a physical partition function" is genuinely
PROVEN: zeros of a ferromagnetic Ising grand partition function lie exactly on |z|=1. Verified
mechanism: (a) a reflection/multiplicative symmetry Z_N(z)=z_1…z_N·Z_N(1/z), PLUS (b) the
coefficients are real and positive (ferromagnetic coupling). Symmetry + positivity → zeros on
the circle.

Map to ξ(s): the functional equation **ξ(s)=ξ(1−s)** IS ξ's reflection symmetry (axis Re=½
instead of |z|=1). So ξ HAS Lee-Yang's ingredient (a). What it lacks is (b) — the clean
coefficient-positivity that ferromagnetic coupling hands Lee-Yang for free. On the ξ side, the
"coefficients" are the prime/arithmetic data, whose positivity is exactly the Weil-positivity
A+B≥P that `OBSTRUCTION_MAP.md` route 5 identifies as THE wall.

**Conclusion (honest):** rigorously surveyed, the physics does NOT hand us a new forcing
mechanism. The ONE physics idea that IS a genuine forcing mechanism (Lee-Yang) needs a positivity
the arithmetic case doesn't have — which is precisely the project's existing wall. Physics points
back to the same wall from a new direction; it does not go around it. Everything else (GUE,
derivative moments, freezing transition, structure factor, quasicrystal) is a CONSISTENCY
signature — necessary if RH holds, but consistent-with-RH is not proof, same ceiling as the 20M
GUE check.

## What's genuinely RUNNABLE and NEW (vs. another confirmation)

- **[NEW, cheap, could surprise] Derivative-moment test (Hughes-Pearce-Crump)** on our verified
  zeros: compute ζ′(½+iγ_n) for the first ~10⁴–10⁵ zeros, form the moment ratios, check
  convergence to the 1/Γ(k+2) RMT prediction for several k. Nobody has published this on a large
  verified set. If it MATCHES → strongest fine-grained RMT confirmation to date. If it DEVIATES →
  genuinely major. Either way it's new data, not a re-run of pair correlation.
- **[confirmation, lower value] Full S(k) structure factor at high resolution + diffuse-background
  decay** on 20M zeros. Mostly re-does Zhang-Torquato; the only fresh angle is quantifying the
  diffuse-background decay rate and asking if it's a diagnostic (likely still circular).
- **[not runnable, conceptual] Freezing-transition physical interpretation** — open question, not
  a computation we can execute to a proof.

## Open questions worth carrying forward
- Does the β=1/√2 freezing transition hold universally across all L-functions (→ hints a
  universal physical system) or is it purely arithmetic?
- Can the diffuse-background decay of S(k) be tied non-circularly to off-line zero contributions?
- Is there a structural reason NO dynamical system with exactly prime-number periodic orbits has
  ever been constructed (research consensus: likely a fundamental barrier, not an engineering gap)?

**Sources (primary, verified):** arXiv:2509.07788, academic.oup.com/qjmath/article/75/4/1481,
arXiv:1801.01541, arXiv:2410.03673 (computational part only; proof claim refuted),
arXiv:2509.07792, people.maths.bris.ac.uk/~mancs/papers/RMTzeta.pdf (Keating-Snaith).
