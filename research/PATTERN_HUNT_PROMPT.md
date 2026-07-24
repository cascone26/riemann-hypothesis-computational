# SESSION PROMPT — Hunt for new systems that share the Riemann-zero spacing pattern

*(Copy everything below the line into a fresh session. It is self-contained — assume the new session knows nothing about this conversation.)*

---

## Your mission

The nontrivial zeros of the Riemann zeta function have local spacing statistics that match the **GUE** (Gaussian Unitary Ensemble) of random matrix theory — the same "level-repulsion" fingerprint seen in heavy-nucleus energy levels, chaotic quantum billiards, microwave/acoustic cavity resonances, etc. This is the Montgomery–Odlyzko phenomenon.

**Your job:** go find *new* lists of numbers — spectra, eigenvalues, gaps, sequences from math/physics/nature/data — and test whether they carry the SAME GUE fingerprint. Prioritize candidates that, to your knowledge, have **not** already been checked against GUE in the literature. This is exploratory: most candidates will fail or be already-known; the goal is to surface anything genuinely new, and to be brutally honest about which bucket each result lands in.

## What "the same pattern" means — the exact tests to run

For any candidate list of ordered real numbers λ₁ < λ₂ < … < λ_N:

1. **Unfold** the spectrum first (critical — do not skip). Remove the smooth trend so the *average* spacing is 1 everywhere. Fit the smooth counting function N̄(λ) (polynomial or local average) and replace each λ_i with N̄(λ_i). All tests run on the unfolded values.
2. **Nearest-neighbor spacing distribution P(s):** histogram the gaps s_i = (unfolded) λ_{i+1} − λ_i. Compare against:
   - GUE Wigner surmise: P(s) ≈ (32/π²) s² exp(−4s²/π)  → level repulsion, P(0)=0, s² rise.
   - Poisson (uncorrelated): P(s) = e^{−s} → NO repulsion, peak at 0.
   - GOE: (π/2) s exp(−π s²/4) → linear repulsion (weaker).
   A GUE match means: P(0)≈0, quadratic rise, peak near s≈0.9.
3. **Nearest-neighbor spacing ratio** r_i = min(s_i,s_{i-1})/max(s_i,s_{i-1}) — the powerful part is this needs NO unfolding. Mean ⟨r⟩: GUE ≈ 0.6027, GOE ≈ 0.5307, Poisson ≈ 0.3863. This is the cheapest, most robust discriminator — compute it first on every candidate.
4. **Pair correlation** R₂(s): should match the GUE form 1 − (sin(πs)/(πs))².
5. **Spectral rigidity** Δ₃(L) and **number variance** Σ²(L): GUE/chaotic → grows like (1/π²)·log L; integrable/Poisson → grows like L. The log vs linear growth is a strong tell.

**Decision rule:** a candidate "shares the pattern" only if MULTIPLE tests agree on GUE (esp. ⟨r⟩≈0.603 AND P(0)≈0 AND log-growth rigidity). One test alone is not enough. Report the ⟨r⟩ value for every candidate even if you run nothing else.

## Reuse the existing tooling — DO NOT rebuild from scratch

This project already implements every test above. Read these first and reuse/adapt them:
- `~/projects/riemann/pair_correlation_fixed.py` — pair correlation vs Montgomery/GUE
- `~/projects/riemann/check_gue_variance.py` — number variance
- `~/projects/riemann/nn_ratio_exact.py` — nearest-neighbor spacing ratio ⟨r⟩
- `~/projects/riemann/exact_gue.py` — exact GUE spacing CDF via Fredholm determinant of the sine kernel (the gold-standard reference curve)
- `~/projects/riemann/analyze.py` / `analyze_2M.py` — full pipeline on the zeta zeros (your working template)
- `~/projects/riemann/data/` — the verified zeta-zero datasets, use as the known-GUE positive control.

**Always run the zeta zeros through your pipeline as a positive control** and a **Poisson/random list as a negative control**, so you know your unfolding + tests are calibrated before trusting any new result.

## Candidate lists to hunt (pick freely; add your own)

Group A — *number theory / math sequences* (cheap, computable locally, and the most likely to be genuinely un-checked):
- Zeros of other L-functions (Dirichlet L-functions, Dedekind zeta) — known GUE, good extra control.
- Eigenvalues of the Laplacian on various domains / graphs (random regular graphs, expander graphs, Cayley graphs).
- Spacings of: sums of two squares, ideal norms in number fields, class numbers, zeros of Hermite/orthogonal polynomials, roots of random/structured polynomials, Salem/Pisot number conjugates.
- Gaps in sequences: prime gaps themselves (known NOT GUE — good negative example), Ulam sequence, Beatty sequences, continued-fraction denominators.
- Eigenvalue spacings of combinatorial/adjacency matrices of real-world graphs.

Group B — *physics / measured spectra* (need to fetch data):
- Molecular/atomic emission line spectra from public databases (NIST Atomic Spectra Database).
- Nuclear energy-level tables (known GUE — control).
- Quasicrystal diffraction peak positions.
- Any open spectroscopy dataset you can pull.

Group C — *"nobody would think to check" wildcards* (highest novelty, run the cheap ⟨r⟩ test first):
- Eigenvalues of large real-world networks (road networks, protein interaction, neural connectome data if public).
- Time-series "spectra": FFT peak spacings of natural signals (seismic, audio, EEG, financial — with heavy skepticism).
- Spacings in astronomical catalogs (exoplanet orbital resonances, asteroid Kirkwood gaps, galaxy redshift combs, pulsar timing).
- Zeros of special functions no one has spacing-tested.

**Be honest about novelty:** before claiming a candidate is un-checked, do a quick literature/web check. If it's already known GUE (L-functions, nuclei, billiards), label it CONTROL, not discovery.

## Honesty gate (mandatory — this project runs on it)

Follow the Charisma Report-vs-Handle discipline:
- **Report** = "we computed the statistics and they match GUE." That is real and valuable — but it is a numerical/statistical observation, NOT a proof of any connection to RH.
- A statistical GUE match ≠ a proof that the system's spectrum *equals* the zeta zeros, and ≠ progress on RH by itself. Say so plainly for every hit.
- Distinguish sharply: (a) already-known-GUE control, (b) new system that statistically matches GUE, (c) the jackpot — a system provably *isomorphic* to the zeta spectrum (extraordinary claim, extraordinary evidence, almost certainly not happening — but flag if you ever think you see it).
- Report finite-N error bars. Small N fakes both directions. State N for every candidate.
- Log negative results too (things that DON'T match) — a clean list of what fails is as useful as what passes.

## Deliverables

1. A results table: candidate | N | ⟨r⟩ | best-fit ensemble (GUE/GOE/Poisson) | which tests run | novelty (control / known / possibly-new) | one-line honest verdict.
2. Plots for any GUE-matching candidate: P(s) vs GUE/GOE/Poisson curves, and Δ₃(L) log-vs-linear.
3. Save everything under `~/projects/riemann/results/pattern_hunt/` with a dated summary `PATTERN_HUNT_<date>.md`.
4. Append a short entry to `~/projects/riemann/LEDGER.md`: what you tested, what matched, what didn't, and the single most interesting genuinely-new candidate (if any) with its proof pointer (the results file).
5. If nothing new turns up, say so honestly — do not manufacture a discovery.

## Start here

1. Read the reuse-tooling files + one data file to calibrate.
2. Run zeta zeros (positive control) + a random list (negative control) through your pipeline; confirm ⟨r⟩≈0.603 and ≈0.386 respectively. Only proceed once controls pass.
3. Batch through Group A first (cheap, local, computable). Compute ⟨r⟩ for everything; deep-test only the ⟨r⟩≈0.60 hits.
4. Then Group B/C as time allows.
5. Write up honestly.
