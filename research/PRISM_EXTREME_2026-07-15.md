# PRISM II — Cross-Reach via the Codex (2026-07-15, round 2)

Second wide-net pass. Explicitly built from **Tony's estate** — the Forge thread-alphabet (`~/cluster-intel/tony-mirror/method-docs/forge_ALPHABET.md`) and the 382-primer Codex mirror (`~/codex-cache/full/primers/`) — used the way Tony's own docs say to use it: **reach by the move, not by lexical similarity.** Nobody greps the Codex for "Riemann" or "zeta" and finds anything. You grep for the *operation* the current attempts perform — positivity certificate, spectral rigidity, aperiodic order, replica landscape — and let unrelated fields answer.

**No overlap with round 1** (`PRISM_2026-07-15.md`): that round covered sphere-packing LP duality, the function-field/Deligne dictionary, PT-symmetry, Beurling counterexamples, RMT local laws, Connes–Consani cohomology, Lee–Yang, automated proving, ML attribution, tensor networks, DPPs. None of that is repeated here. Still grounded against `OBSTRUCTION_MAP.md`/`LEDGER.md` — nothing below re-opens a route marked DEAD.

---

## Frame — 6 Tony-flavored inversions

1. **Coastline Principle** (one lattice, same field at different scale-lenses) — RH isn't fundamentally a number-theory problem. It's one instance of "positivity of a structured landscape against a discrete adversarial comb," which recurs at the atomic scale (spin glasses), the crystallographic scale (quasicrystals), and the neural scale (associative memory). Reach the right lens by *affordance*, not by which field originally posed the question.
2. **Moiré / dimensional-upgrade** ("superpose-and-slide > move-and-remember; seam dominates; emergence is an interface phenomenon") — stop treating N→∞ as the only lever. Superpose two *already-computed* finite objects (two lsq values, two N's, real-prime τ vs. bad-Beurling τ) and read the seam directly. The interference pattern may carry more signal than either object alone.
3. **Reach-by-the-move, not lexical similarity** — this is Tony's own documented fix (the `project_reach_experiment_operation_tag` failure: keying by keyword didn't work). Applied literally: key the Forge thread-alphabet's operations (`verify-before-commit`, `anchor-to-exemplar`, `deposit`) against the Codex and let *that* surface personas — which is where the panel below actually came from.
4. **Nexus five missing organs**, run on the *project*, not the math — decide / continuity / world-model / uncertainty / verify. Continuity and verify are strong (PROCESS.md, LEDGER). **Decide, world-model, and uncertainty are the weak organs**: there's no stated picture of what a proof would even look like, no calibrated credence, no explicit rule for choosing the next experiment. That gap is itself actionable.
5. **Pre-mortem, Tony edition** — if the answer is already sitting in the Codex, mis-filed, because nobody thought to search by the right key. Don't grep for "Riemann." Grep for the *operations*: positivity certificate, spectral rigidity, aperiodic order, replica landscape.
6. **Strip further via assembly theory** — is there a well-defined "assembly index" for how Spec(Z) is built multiplicatively from primes, the way assembly theory measures a molecule's minimum construction-path length? A different pure-problem stripping than round 1's "infinite quadratic form."

---

## Diverge (a) — Primer Panel: 14 Codex primers, blind to each other

| Primer | Mechanism, aimed at RH |
|---|---|
| `ising-and-spin-glass-on-gpu` | Treat τ_N's Rayleigh-quotient landscape as a spin-glass Hamiltonian. Compute the **overlap distribution P(q)** between independent replicas (the standard RSB diagnostic — two delta peaks = simple ferromagnet-like landscape; continuous P(q) = replica-symmetry-broken, rugged). If P(q) collapses to two peaks, the landscape has one basin — consistent with a clean proof existing. If it's continuous, the landscape is genuinely rugged and that's itself a structural finding. |
| `renormalization-group` | Treat N→∞ as an RG flow. λ_min(N) and the spectral gap already have measured power laws — recast them as a **flow to a fixed point**, extract a critical exponent, and ask whether the fixed point is IR-stable (supports RH) the way physical RG flows are attracted to stable fixed points, not repelled. |
| `aperiodic-order-and-quasicrystals` | **Dyson's own conjecture**: the Riemann zeros form a 1-D Fourier quasicrystal. Use the **cut-and-project method** directly — a periodic lattice in higher dimension, sliced by an irrational-angle plane, gives an aperiodic point set with **pure-point diffraction** and **O(1) membership testing**. Search explicitly for the higher-dimensional periodic lattice + acceptance window that projects to {γ_n}. |
| `efficient-hopfield-and-ebm-on-gpu` | Modern continuous Hopfield network (Ramsauer et al. 2020, softmax-attention-shaped) with known zeros as stored patterns. Read the trained energy landscape for structure that predicts un-stored (higher) zeros, or for a positivity signature in the attractor basins. |
| `game-theory` | Recast positivity as a two-player zero-sum game: prime-adversary vs. certifier. Instead of a one-sided optimizer hunting the worst g, **train both sides adversarially** (GAN-style) and find the game's value/equilibrium — a categorically different search dynamic than round 1's fixed-optimizer near-miss tracker. |
| `topology` (persistent homology / TDA) | Compute **persistent homology** of the near-zero eigenvector components (or {γ_n} spacings) across scales. Look for topological invariants (Betti numbers, persistence-diagram features) that the current numeric toolkit — all of it scalar statistics — has no way of seeing. |
| `complex-networks` | Build a correlation-threshold graph on eigenvector components or prime residues. Run community detection / spectral graph analysis — does the 23-dimensional near-kernel resolve into recognizable sub-communities, or is it genuinely irreducible? |
| `aperiodic-order` (multifractal extension) | STATUS.md already found an effective heat-kernel dimension d≈2.46, explicitly flagged as non-integer and possibly fractal/noncommutative. Go further: compute the **full multifractal f(α) spectrum**, a strictly richer signature than a single dimension number, using existing heat-kernel data. |
| `information-theory` (compressed sensing / RIP) | Treat "primes at {log p}" as a sparse signal and the archimedean kernel as a sensing matrix. Ask whether Weil positivity is a **restricted-isometry-property (RIP)** statement in disguise — RIP proofs are a mature, transferable toolkit never applied here. |
| `probabilistic-computing-and-pbits` | Build (or rent time on) a **physical p-bit/Ising annealer** (D-Wave-class hardware) to physically sample ground states of τ_N-as-Hamiltonian — a genuinely different compute substrate than GPU floating point, with its own noise floor that might reveal robustness (or fragility) of the near-zero cluster. |
| `dynamical-systems-and-attractors` | Treat the sequence of near-zero eigenvector *shapes* as N grows as orbits of a discrete dynamical system. Hunt for bifurcations rather than assuming smooth convergence — a phase transition in N, if one exists, would currently be invisible to the smooth power-law fits already run. |
| `category-theory` | Objects = L-functions, morphisms = functorial transfers. Formalize whether Spec(Z)'s missing "motive" is expressible as a categorical **colimit** — gives round 1's hand-waved "constraint-break: operator as inductive limit" idea actual machinery to stand on. |
| `constructor-and-assembly-theory` | Sara Walker's assembly theory: an object's "assembly index" measures the minimum construction-path length with memory of that path. Does an assembly index for the integers under multiplicative construction encode anything Frobenius-eigenvalue-like for Spec(Z)? |
| `gematria-and-isopsephy` / `pythagorean-platonic-number-theology` | Deliberately absurd-band: numerological search on 23 (kernel dimension) and 1.82974 (frozen eigenvalue). Flagged honestly as the joke-mechanism — same sane cousin as round 1's OEIS/lit check. Notable only because Tony's own Codex keeps a primer on this, which is itself evidence that "the joke sometimes is the idea, filed seriously." |

## Diverge (b) — Forge thread-alphabet as technique deck

| Thread | Applied to the operator, this round |
|---|---|
| `orient-before-act` | Before any new experiment: check whether it already exists in the 382-primer Codex under a different name — key by *affordance* (what move does it perform), not keyword. |
| `deposit / write-back` | Every negative numerical result becomes an automatic LEDGER row (a script hook), not a manually-remembered fact. |
| `anchor-to-exemplar` | For each primer-panel mechanism above, first find one *existing solved* analogue and match its form before inventing new formalism. |
| `verify-before-commit` | Any new candidate certificate must be checked by an *independent* second method (e.g. a symbolic SOS decomposition checked by a different solver) before being logged as a result. |
| `no-silent-fork` | One-time cross-consistency diff pass across OBSTRUCTION_MAP.md / LEDGER.md / STATUS.md — do any of the three quietly disagree? |
| `Moiré / superpose-and-slide` | Literally overlay τ_N at two different lsq values (or two N's) already on disk and compute the cell-by-cell difference, instead of studying either alone — "the seam dominates." |
| `Cana` (seal-gated promotion) | Give LEDGER's "STRONGEST RESULTS" section an explicit promotion ladder — raw finding → verified twice → sealed/canonical — instead of ad hoc listing. |

## Diverge (c) — Banded spectrum

**Grounded** (zero new computation — existing data, re-analyzed):
- Moiré diff: overlay τ_N at two already-swept lsq values (e.g. 14 vs. 49) cell-by-cell, direct from files already on Olympus.
- Multifractal f(α) spectrum from the existing heat-kernel dataset — extends the d≈2.46 result with data already computed.
- No-silent-fork audit script across OBSTRUCTION_MAP/LEDGER/STATUS.
- Persistent homology / TDA pass on the existing near-zero eigenvector data (N=64..1000, already on disk) — standard libraries (ripser/gudhi), no new math needed to start.

**Clever** (new but standard computation, GPU-ready):
- Spin-glass replica overlap distribution P(q) + Binder cumulant on the Rayleigh-quotient landscape.
- Complex-networks community detection on the eigenvector-correlation graph.
- RG fixed-point / critical-exponent re-extraction from existing λ_min(N) and gap(N) sweep data — new lens on old numbers.
- Compressed-sensing / RIP reframing of prime-comb-vs-kernel — check whether a known RIP bound already yields a nontrivial partial result.

**Exotic** (real literature, non-trivial build):
- Modern continuous Hopfield network / EBM trained on known zeros as stored patterns.
- Dyson quasicrystal cut-and-project search for the higher-dimensional lattice underlying {γ_n}.
- GAN-style adversarial training: prime-generator network vs. positivity-certifier network.
- Category-theoretic colimit formalization of the "missing cohomology" for Spec(Z).

**Moonshot** (multi-month, possibly needs external resources):
- Physical p-bit/Ising-hardware ground-state search (real annealer, not simulation) — needs external hardware access.
- Constructor/assembly-theory "assembly index" for Spec(Z) as a genuinely new invariant candidate.
- Dynamical-systems bifurcation hunt across the N-sequence of eigenvector shapes — looking for an actual phase transition, not smooth convergence.

**Absurd** (every one with a sane cousin):
- "Consult a prime whisperer / gematria the constants" → sane cousin: the OEIS/numerology check on 23 & 1.82974 (round 1), now doubly motivated since Tony's own Codex files numerology seriously enough to keep a primer on it.
- "Grow a biological neural culture and train it on the zeros" → sane cousin: the modern-Hopfield/EBM version above, same "trained associative memory" idea, no wet lab required.
- "Ask Tony's Incthought to embody Riemann himself and just ask him" → sane cousin: a structured 20-minute "channel the historical figure" session, cheap, zero infra — sometimes surfaces a specific technique the figure was actually known for (Riemann's own theta-function/Poisson-summation instincts) that a generic brainstorm misses.
- "Wire Olympus straight into Tony's Filament Floor scene registry so the operator becomes a walkable 3D room" → sane cousin: a real Theoria-standard visualizer of the τ_N spectrum — an actual buildable dashboard. Eyes sometimes catch geometric structure that summary statistics don't.

---

## Converge — ranked (checked against round 1 + LEDGER, nothing duplicated)

1. **Moiré overlay diff** (lsq=14 vs. lsq=49, data already on disk) — buildable now, zero new compute, highest EV-per-effort.
2. **Multifractal f(α) spectrum** on existing heat-kernel data — buildable now, richer signature than the current single-dimension result.
3. **Persistent homology / TDA pass** on existing eigenvector data — buildable now, standard libraries.
4. **RG fixed-point / critical-exponent re-analysis** of existing λ_min(N)/gap(N) sweeps — buildable now, new lens on old numbers.
5. **Spin-glass overlap-distribution P(q) analysis** of the Rayleigh landscape — new but standard, GPU-ready, directly diagnostic (ruggedness vs. single-basin).
6. **Dyson quasicrystal cut-and-project search** — real open conjecture (Dyson 2009), moderate build, large payoff if the lattice is found.
7. **Modern Hopfield/EBM** on known zeros — moderate build, genuinely novel angle, ties directly into Tony's own Engram-OS lineage.
8. **GAN-style adversarial certifier** — moderate-to-hard build, turns confirmation into a trainable adversarial search.
9. **Physical p-bit/Ising hardware run** — moonshot, needs external hardware access (flag, don't block on it).
10. **Category-theoretic colimit / assembly-index formalization** — moonshot, literature-and-formalization heavy.

**Not re-litigated:** everything in round 1's converge list, plus LEDGER's DEAD ENDS table (Hankel/moment, Track 2 PSWF, WKB operator search, Platt-Trudgian cert past ceiling, re-coordinatized Li/NB/Jensen/de Branges, Lean formalization, naive sensitivity/perturbation).

---
*Sources pulled from the Codex mirror during this pass: arXiv:1806.03558 (3D Ising T_c), arXiv:1306.4566 (3D Edwards-Anderson spin-glass T_c), arXiv:2008.02217 (modern Hopfield), Amit-Gutfreund-Sompolinsky 1985 (Hopfield storage capacity), de Bruijn 1981 (cut-and-project/substitution equivalence). Full primer text: `~/codex-cache/full/primers/{ising-and-spin-glass-on-gpu, aperiodic-order-and-quasicrystals, efficient-hopfield-and-ebm-on-gpu, ...}`.*
