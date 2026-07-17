# RH and Real Physics — what's proven vs. what's open (2026-07-16)

Follow-up to the three PRISM rounds (2026-07-15). Those rounds treated physics mostly as
*analogy* (borrow a technique from spin glasses, quasicrystals, ecology). This note is
different: it separates what is **already rigorous, established physics** connected to
ζ(s) from what is **genuinely open experimental territory** worth building. Nothing here
duplicates `connes_spectral_triples_2025.md` (which covers Connes' 1999 trace-formula /
adele-class-space spectral realization and Track 2 of `OBSTRUCTION_MAP.md`) — this note
covers the *statistical-mechanical* and *quantum-chaos* physics upstream of that program,
which the project has not yet written up.

## 1. Rigorous — not analogy, actual theorems

### 1.1 The primon gas (Julia, 1990)
A free Bose gas where each prime p is a single-particle state with energy log(p). By
unique factorization, the n-th multi-particle state has energy log(n) with degeneracy 1.
The grand canonical partition function at inverse temperature β is, exactly:

  Z(β) = Σ_{n≥1} n^{-β} = ζ(β)

This is elementary (it falls straight out of the Euler product), but it means: ζ(s) *is*
literally a partition function of a real statistical-mechanical model, not a formal
trick. The system has a Hagedorn-like divergence at β=1 — exactly the pole of ζ. Also
called the "free Riemann gas."

### 1.2 The Bost–Connes system (1995) — the rigorous version
Bost & Connes ("Hecke algebras, type III factors and phase transitions with spontaneous
symmetry breaking in number theory," *Selecta Math.* 1995) built a genuine C*-dynamical
system (A, σ_t) from the Hecke algebra of the ax+b group over Q and the adeles, whose
partition function is exactly ζ(β). Unlike the primon gas this is a fully rigorous
quantum statistical mechanical system with a real phase transition:
- β ≤ 1: unique KMS state (equilibrium is "boring").
- β > 1: a whole simplex of KMS states — **spontaneous symmetry breaking**, the
  mathematical-physics term for exactly what happens in a real ferromagnet below T_c.
- The symmetry being broken is the idele class group / Gal(Q^ab/Q) ≅ Ẑ*, acting on the
  extremal low-temperature states. Evaluating "ground states" (β→∞) recovers explicit
  class field theory for Q (roots of unity).

This is the direct historical ancestor of Connes' 1999 spectral-realization program
already logged in `connes_spectral_triples_2025.md` — worth reading as the "physics
prequel" to that file, since it's where the phase-transition-as-arithmetic-symmetry-
breaking idea originates.

**What this does NOT give you:** neither system's spectrum is the zeta zeros. Both give
you ζ(s) as a *partition function* (a sum over states), which is a different object from
an *operator whose eigenvalues are the zeros* (Hilbert–Pólya). Real, rigorous, but not
by itself a route to RH — it's a genuine physical realization of the function, not of
the zero set.

### 1.3 GUE statistics — the Montgomery–Odlyzko–Dyson correspondence
Montgomery's pair correlation conjecture (1973, proved for restricted test-function
support) says the zeta zeros' pair correlation matches that of Gaussian Unitary Ensemble
eigenvalues. Famous history: Dyson recognized the formula over tea with Montgomery at
the IAS in 1972. Odlyzko's numerics (1987 onward, now to zero #10^23-scale) confirm
nearest-neighbor spacing matches GUE to high precision.

This statistic is not a math curiosity — it is *the same statistic measured in real
physical systems*: heavy-nuclei slow-neutron resonance spacings (Wigner's original
surmise, 1950s) and later, any quantum system whose classical limit is chaotic with
broken time-reversal symmetry (Bohigas–Giannoni–Schmit conjecture, 1984, now checked in
many real experiments: chaotic microwave cavities, chaotic quantum billiards, Rydberg
atoms in strong fields). The zeros behave statistically like the energy levels of a
real chaotic quantum system — this is measured physics, not metaphor.

**What this does NOT give you:** GUE statistics is a *necessary* signature if a
Hilbert–Pólya operator exists, but matching statistics is consistent with RH, not a
proof of it — a huge random Hermitian matrix also has GUE statistics without being
"the" zeta operator. Statistics alone under-determines the operator.

## 2. Genuinely open — buildable, not yet done

These aren't restatements of 1.1–1.3, they're actual next actions if the goal is new
*physical* evidence rather than another numerical sweep on Olympus.

1. **Physically build a PT-symmetric candidate.** The Bender–Brody–Müller non-Hermitian
   Hamiltonian (already flagged PRISM round 1, item 3) doesn't have to stay a simulation.
   PT-symmetric systems have been physically realized for over a decade as coupled
   optical waveguide arrays and microwave networks (Christodoulides group and others,
   "Observation of parity-time symmetry in optics," *Nat. Phys.* 2010, and follow-ups).
   If that specific Hamiltonian can be cast as a photonic lattice, its spectrum becomes
   something you *measure optically*, not compute — an independent check with a
   completely different error model than floating point.

2. **Quantum graphs as microwave networks.** A network of wires with lengths log(p) for
   each prime is physically buildable (quantum-graph microwave-network experiments are
   an established technique, Kottos–Smilansky theory + Hul/Sirko-style realizations).
   Its measured spectrum is real physics, not a matrix diagonalization — a different
   kind of test of "primes as geometry."

3. **Measure Lee–Yang zeros of a primon-gas-like system directly**, using the real
   experimental technique of Loschmidt-echo / interferometric measurement of a probe
   qubit's decoherence (demonstrated on NMR and trapped-ion platforms for small spin
   models). Applying that measurement method to any physical stand-in for the primon
   gas/Bost–Connes partition function would be *measuring* zeros of a partition function
   with a quantum device, rather than computing them.

4. **Run the existing operator on real quantum hardware, not just Olympus GPUs.**
   Everything in Track 2 so far (`connes_decisive_v3.py`) is classical linear algebra.
   Whether a physical quantum processor's measured spectrum agrees with τ_N is a
   qualitatively different kind of evidence than another floating-point sweep — it
   tests whether the PSD property survives being realized as an actual physical system
   rather than a simulated one.

## Converge — ranked (buildable-now first)
1. Bost–Connes / primon-gas literature deep-read — cheap, clarifies exactly what's
   rigorous vs. hoped-for in the physics-of-ζ program; directly upstream of
   `connes_spectral_triples_2025.md`.
2. GUE-statistics cross-check on τ_N's own spectrum (not just the zeros) — buildable now
   on Olympus, reuses existing eigenvalue data, asks whether the *operator's* spectrum
   (not just known zeros) shows the same universality class.
3. PT-symmetric photonic-lattice feasibility check — literature/design study first,
   before any hardware.
4. Quantum-graph microwave-network design study — same, feasibility before hardware.
5. Quantum-hardware run of a small τ_N block — moonshot, highest cost, most novel
   evidence type if it lands.

**Not re-litigated:** anything from PRISM I–III's converge lists or LEDGER's dead ends
(Hankel/moment, Track 2 PSWF, WKB search, Platt–Trudgian ceiling, re-coordinatized
Li/NB/Jensen/de Branges, Lean formalization, naive perturbation).
