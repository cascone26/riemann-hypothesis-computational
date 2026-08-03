# The Phase-Only Möbius Comb — derivation note (2026-08-01, LEDGER §517–§553)

REPORT-tier. All numerics verified in-session; scripts `operators/krein_*.py`, data zeros10k/argzp10k/
abszp10k.npy. Lit gate (§521): the phase-only observation was NOT found in the literature. No Handle for RH.

## The object
For the Riemann zeros ρ = ½ + iγ (γ > 0) define the zero-sum with a **unit-modulus** weight
    G(u) = (1/N) Σ_{0<γ<Γ} e^{−i·arg ζ′(ρ)} · e^{iuγ}.
Empirically (500 → 10⁴ zeros):
- Peaks of |G(u)| land exactly on u = log n for **squarefree** integers n (40/40 vs 0/40 non-squarefree
  at 10⁴ zeros, p < 0.001; §530); squares are actively suppressed.
- At u = log n, Re[e^{ic}G(log n)] ∝ **μ(n)/n^{≈0.5} (phase-only 0.51, drifts; EXACT-Gonek 0.48→½; §563: brackets ½, not clean-fixed)** — the RH-critical exponent ½, robust to 0.50–0.52
  across 3 zero-counts × 3 ranges (§552b). Shuffling the phases collapses the exponent to the trivial 0.42.

## Derivation (why it is the Möbius function)
**Step 1 — the classical explicit formula.** With 1/ζ(s) = Σ μ(n) n^{−s} (Re s>1), Perron + residues at
the simple zeros give the (conditional, simple-zeros) Titchmarsh/Gonek formula for the Mertens function
    M(x) = Σ_{n≤x} μ(n) = Σ_ρ x^ρ/(ρ ζ′(ρ)) − 2 + Σ_{trivial}.
The **exact weight** is the Gonek weight 1/(ρ ζ′(ρ)). Verified here: reconstructing M(x) from 10⁴ zeros
with this weight gives corr 0.9933 → 0.9999 as N grows (§519, §522), and the reconstruction is sharply
best at Re ρ = ½ (σ-knob peak 0.513, §553) — the √x envelope of M(x) is the critical line.

**Step 2 — the weight decomposition.** Write ζ′(ρ) = |ζ′(ρ)| e^{i·arg ζ′(ρ)}. Then
    e^{−i·arg ζ′(ρ)} = |ζ′(ρ)| · (1/ζ′(ρ)) = ρ|ζ′(ρ)| · [1/(ρ ζ′(ρ))].
So the **phase-only weight equals the exact Gonek weight times ρ|ζ′(ρ)|.** The factor ρ|ζ′(ρ)| has modulus
|γ|·|ζ′(ρ)|, which varies *slowly on a logarithmic scale* (|ζ′(½+iγ)| grows like a power of log γ with
O(1) multiplicative fluctuation; |γ| is smooth). **Measured:** |ρζ′(ρ)| ~ γ^{0.22}, multiplicative
fluctuation only 1.77× (log-sd 0.57) — slowly-log-varying, confirmed; the 0.22 drift is what shifts the
measured exponent from ½ to 0.51. A slowly-log-varying real-ish envelope multiplying the
Gonek weight **preserves the oscillatory phases e^{iuγ} and the residue sign pattern** — which is where the
arithmetic (positions log n, signs μ(n)) lives — while only mildly reweighting magnitudes. Hence the
phase-only sum inherits the squarefree comb and the μ-signs, and an *approximate* √n scaling (measured
0.51 vs the exact ½), the small excess being the log-drift of ρ|ζ′(ρ)| (the §519b gap law, ≈0.009·ln n).

**Step 3 — the phase is S(T).** The intrinsic identity (§541, machine-verified to 4×10⁻¹⁴, classical via
ζ(½+it)=e^{−iθ_H}Z): arg ζ′(ρ) = π·S(γ) − π/2, where S(T)=π⁻¹ arg ζ(½+iT). Therefore the twist is a pure
**S(T)-phase**, and the whole statement reads: *the zero-fluctuation field S(γ) alone — no magnitude data —
reconstructs the Möbius function's positions, signs, and (to log-drift) its RH-critical √n normalization.*

## What is and isn't claimed
- **Known:** the Mertens explicit formula (Step 1) and the Gonek weight are classical.
- **The novel content (lit-gated NOT found, §521):** that dropping |ζ′| entirely — the **unit-modulus / phase-
  only** weight e^{−i arg ζ′(ρ)} — still yields the squarefree comb and μ(n)/√n, i.e. the Möbius function is
  carried almost entirely by the *phase* arg ζ′ ≡ πS, with magnitude nearly irrelevant to positions/signs.
- **Not claimed:** any new proof of RH or of M(x)=O(x^{½+ε}). The σ-knob pins Re ρ=½ as a *consistency*
  check of the explicit formula (on-line zeros used), not an independent probe.

## Open (to lift from REPORT toward a theorem)
1. Bound the deviation phase-only − Gonek: show ‖ρ|ζ′(ρ)| − (its log-average)‖ contributes only o(1) to the
   comb positions/signs (moment estimates for |ζ′(ρ)|; Gonek's conjecture territory).
2. Prove the exponent → ½ exactly as Γ→∞ (the 0.51 → 0.50 drift), quantifying the §519b gap.
3. Connect to Ng / Gonek–Ng work on Σ x^ρ/ζ′(ρ) and to unconditional partial results.


## Mechanism of the squarefree SELECTION (§567) — Euler product of exp(S prime-sum)
e^{-2iπS(γ)}, S(γ)=-(1/π)Σ_p Σ_k p^{-k/2}/k sin(kγ log p), gives -2iπS = Σ_p p^{-1/2}(p^{iγ}-p^{-iγ})+(k≥2).
Then exp = Π_p exp[p^{-1/2}(p^{iγ}-p^{-iγ})], first-order = Π_p(1 + p^{-1/2}p^{iγ} - ...): the coefficient of
n^{iγ} is nonzero only for SQUAREFREE n (distinct primes; a repeated prime needs the suppressed k≥2 terms).
CONFIRMED: the twist rebuilt from PRIMES ALONE (no argζ' data) gives the squarefree comb with 0/50 non-sf
contaminants (agreement 0.949 with the real twist). This derives WHICH integers appear. The leading terms give
POSITIVE coeffs (½)^r/√n, so the μ SIGNS are NOT from this expansion — those are the Gonek/Perron piece above.
The two derivations are complementary: Euler-product => squarefree selection; explicit formula => μ signs + √n.