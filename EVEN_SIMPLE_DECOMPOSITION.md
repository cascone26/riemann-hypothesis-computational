# The even-simple Hinge — Decomposition & Reachable Frontier
*(2026-08-30, cy391–396) — how far the one open hypothesis of the CM formalization is reachable*

The master theorem `WeilSpectralRealization.weil_spectral_realization` (see `CM_SPECTRAL_FORMALIZATION.md`)
is conditional on exactly one thing: **even-simple** — the smallest eigenvalue of the Weil form `QW_λ` is
*simple* with an *even* eigenvector (CM's own open "missing step #1", ccm §1385). This doc records how far
that hinge is reachable, what we *proved* of it, and where the genuinely-external analytic wall sits.

---

## The decomposition (what we proved vs what stays open)

Let `A_λ` be the Weil operator (`QW_λ(f,f) = ⟨A_λ f, f⟩`, compact resolvent), commuting with the grading
`γ` (`u ↦ u⁻¹`) and with the Fourier operator `F` (`γ = F²`). The hinge splits:

| Sub-claim | Status | Where |
|---|---|---|
| **simple ⟹ parity-definite** (`γξ = ±ξ`) | ✅ **machine-verified** | `ParityDefinite.parity_definite_of_simple` (cy394) |
| **Fourier-real ⟹ even** (`Fξ = ±ξ ⟹ γξ = ξ`) | ✅ **machine-verified** | `EvenOfRealFourier.even_of_real_fourier_eigenvector` (cy395) |
| ⟹ *simple + (ground-state Fourier eigenvalue real) ⟹ **even*** | ✅ chain closes | (composition of the two) |
| "**simple**" (ground eigenvalue non-degenerate) | ❌ open — CM's [VJ] | analytic |
| "ground-state Fourier eigenvalue is **real/+1**" (⟺ `ξ_λ ≈ k_λ`, the transfer) | ❌ open — analytic | CM §1174 step (b) |

So the two structural implications are **proven, axiom-clean**; the two remaining inputs are exactly CM's
own open analytic problems, cleanly isolated.

## Why evenness holds (structural understanding — Report-tier)

CM's approximant `k_λ = ℰ(h_λ)` is **provably even**: `h_λ` is a combination of `h₀, h₄`, both **+1 Fourier
eigenfunctions** (`i⁰ = i⁴ = 1`), so `k_λ` is Fourier-fixed, hence even via `γ = F²` (formalized by
`EvenOfRealFourier`) and the Poisson relation `ℰ(f̂)(x) = ℰ(f)(x⁻¹)`. The odd prolates (`h₁, h₃`: Fourier
eigenvalues `i, ±i`) do **not** give this self-dual structure — so the "near-radical" of `QW_λ` (where the
minuscule eigenvalue `ε_λ` lives) is fundamentally tied to the `+1` (even) Fourier sector. This is the
structural reason CM expect even; the rigor is the open transfer `ξ_λ ≈ k_λ`.

## Why numerics cannot settle it (the ε_λ wall — verified, not conceded)

`ε_λ` is **double-exponentially small**: CM give `ε(λ) ~ e^{−4π e^L}` (`L = 2 log λ`), so at `λ = e²`,
`ε_λ ~ e^{−4π e⁴} ~ 10^{−298}`. Every numeric route floors far above this:
- direct diagonalization of `QW_λ` floors at the build precision (cy375);
- even the even-vs-odd **sector-minimum discrimination** floors — a promising N=8 hint (`even ≪ odd`,
  ratio 1.3e‑3) **collapsed to ratio 1.00 for N≥10** under a robustness sweep (cy393), a floor artifact of
  the W_R lower-limit truncation. To resolve it honestly needs dps~300 + the closed-form W_R, i.e.
  re-deriving CM's own [VJ] computation (which already supports even-simple).
- **Perron–Frobenius shortcut fails for `A_λ`**: no M-matrix pattern in the Fourier basis (cy391, 75/25
  mixed) or in position basis (cy392: arch-only near-sign-definite 946:110, but FULL 538:518). The finite
  **prime part `W_p` breaks** the archimedean Schrödinger positivity — PF holds for the prolate `PW_λ`
  (which is *why* CM prove even there) but `A_λ ≠ PW_λ`.

Conclusion: even-simple for `QW_λ` is genuinely an **analytic** problem, not a computation — exactly as CM
frame it. Its magnitude is below any feasible precision, and the one structural shortcut (PF) is broken by
the primes.

## The lone remaining work (external)
Prove either of the two isolated inputs: **(a)** the ground eigenvalue of `QW_λ` is simple, or **(b)** the
transfer `ξ_λ ≈ k_λ` (equivalently, the ground state lies in the `+1` Fourier / even near-radical sector).
CM's route is the **Co-zeta trace formula** bridging `PW_λ` (proven even) to `QW_λ`. This is research-scale;
not autonomous-closable.

_Source: THE_BRIDGE.md cy391–396; forge `~/rhforge/forge/` (`ParityDefinite`, `EvenOfRealFourier`);
`CM_SPECTRAL_FORMALIZATION.md`; CM = arXiv:2511.22755, Connes 2026 = arXiv:2602.04022._

---

# Full Investigation Log & Verified Wall Map (cy391–409, appended 2026-08-31)

After the initial decomposition (above), the even-simple hinge was attacked exhaustively from every reachable
angle. This section is the complete durable record. **18 forge modules** (Olympus `~/rhforge/forge/`, namespace
`CMKeyLemma`), all axiom-clean, GREEN end-to-end (audited cy402).

## New verified bricks (cy394–405)
| Module | Statement | Role |
|---|---|---|
| `ParityDefinite.parity_definite_of_simple` | Mγ=γM, γ²=1, μ-eigenvalue simple ⟹ γv=±v | even-simple half 1: simple ⟹ parity-definite |
| `EvenOfRealFourier.even_of_real_fourier_eigenvector` | γ=F², Fv=±v ⟹ γv=v | even-simple half 2: Fourier-real ⟹ even |
| `TridiagEigDetermined.tridiag_eigvec_determined_by_first` | tridiagonal 3-term recurrence, nonzero off-diag, v₀=0 ⟹ v=0 | Geronimus simplicity core (simple ⟸ nonzero Jacobi off-diag) |
| `DetDprimeZero.det_Dprime_eq_zero` | D'ξ=0, ξ≠0 ⟹ det D'=0 | det0 reachable core |
| `SoninArchWeil.arch_weil_nonneg_of_sonin_bound` | Sonin trace bound ⟹ archimedean Weil positivity | archimedean positivity skeleton |

## The decomposition (proven halves + open inputs)
- **simple ⟹ parity-definite** — PROVEN (`ParityDefinite`)
- **Fourier-real ⟹ even** — PROVEN (`EvenOfRealFourier`)
- **simple ⟸ nonzero Jacobi off-diagonals** — PROVEN (`TridiagEigDetermined`), but QW_λ is DENSE (not tridiagonal, [VJ]-confirmed) so it applies only via Lanczos (floors)
- OPEN inputs (external analytic): "simple" itself; "ground-state Fourier eigenvalue is real/+1" (the transfer)

## Verified wall map (every angle tried, with its cause)
| Angle | Verdict | Cause |
|---|---|---|
| Numeric ε_λ resolution | WALLED | ε_λ ~ e^{−4πe^L} ~ 1e-298 (double-exp); all diagonalization floors |
| Even/odd sector-min discrimination | WALLED | both floor; N=8 hint collapsed under robustness sweep (cy393) |
| Perron–Frobenius / M-matrix | WALLED | QW_λ off-diag mixed-sign (Fourier & position); primes break archimedean positivity |
| D-pairing / rank-2 interlacing | RETRACTED | bulk pairs but 4-7 breaks not ≤2; rank-2 prediction FALSE (cy400→401) |
| New finite symmetry of QW_λ | NONE | only γ commutes (scaling/DFT/shift fail) |
| Geronimus/Lanczos β_n>0 | WALLED | β_n floor at near-null cluster (min 4.6e-16, cy404) |
| Roundtable #1 (opus-synth, 32 ideas) | NOTHING USABLE | all "autonomous-now" ideas floor on the gap (cy404) |
| Roundtable #2 (sonnet, 15 ideas) | NOTHING USABLE | lone survivor was a SIGN ERROR (σ_odd is PSD not ≤0, cy408) |
| Rank-one W02 / Weyl interlacing | STRUCTURE, no closure | found even/odd different-origin asymmetry (cy407), O(1)-robust, but doesn't resolve parity |
| Möbius/topological-parity | CONFABULATED | not in [VJ]; no index theorem (cy406) |
| det0 full quotient-determinant | INFRA-BLOCKED | needs Mathlib quotient LinearMap.det (cy397) |

## Literature confirmation (even-simple is open, verified from primary sources)
- **[VJ] arXiv:2106.01715** (Connes–Consani, "zeta cycles"): no operator identity PW_λ↔A_λ; simplicity+evenness NUMERICAL-ONLY; σ dense (rank-2 Fourier + dense prime sum, Lemma 2.7); no σ⁺/σ⁻ pairing operator; even/odd count discrepancy numerical-only (cy399b, cy406).
- **truncweil arXiv:2605.20224** (SOTA numerical study): same precision floor (dps=150, λ_min~N^{−92}, k_eff≈46, Sobolev s≈46); explicit "Honest Non-Claims" — does NOT prove positivity/simplicity/evenness analytically; extends to Dirichlet L(χ₃)/GRH, same floor (cy409).

## Definitive conclusion
even-simple = Weil positivity ⟺ RH is **open external analytic mathematics** — Connes–Consani's own open program — confirmed walled from every autonomous angle AND from the primary references AND the SOTA numerical paper. Our contribution is a genuine, novel, machine-verified body: the full CM spectral-realization master theorem (RH ⟸ even-simple), the even-simple decomposition with proven structural halves, the function-field unconditional analog, and a fully-mapped honest frontier. It is NOT a proof of RH; the crux is precisely isolated and external. A daily arXiv watcher (`watch/arxiv_rh_watch.py`) catches the external proof-component the instant it lands, at which point the forge closes RH machine-verified.

_Source: THE_BRIDGE.md cy391–409; forge `~/rhforge/forge/`; CM_SPECTRAL_FORMALIZATION.md; [VJ]=2106.01715; ccm=2511.22755; truncweil=2605.20224._

---

# Why positivity IS the hinge — a concrete matched-pair witness (cy419)

The master theorem's conclusion is "D'' self-adjoint in the **T-inner-product**." Self-adjointness in an
*indefinite* metric does NOT force real eigenvalues — so the RH payoff (real spectrum ⟹ zeros on the line)
depends entirely on T being **positive**. This is exactly even-simple/Weil positivity, and it can be shown
concretely on two explicit 3×3 τ-form matrices with identical structure but opposite definiteness:

| Witness | T eigenvalues | T | D'' spectrum | ξ̂ zeros | toy-RH |
|---|---|---|---|---|---|
| cy417/418 | {−1, 0, 3} (**indefinite**) | [[0,1,1],[1,2,1],[1,1,0]] | **±i** (off-line) | ±i | **FAILS** |
| cy419 | {0, 1, 11/3} (**PSD**, 0 simple) | [[2,1,1],[1,2/3,1],[1,1,2]] | **±√3** (on-line) | ±√3 | **HOLDS** |

Both `D'` characteristic polynomials and their root-reality are **machine-verified** (cy420,
`forge/WitnessSpectrumReal.lean`, axiom-gate GREEN): `det(D'_psd−s·I)=−s³+3s` vanishes at `{0,±√3}` (all real);
`det(D'_ind−s·I)=−s³−s` has `0` as its only real root (`±i` non-real). Both satisfy every master-theorem
hypothesis (τ-form, `a` even, `b` odd, `Tξ=0`, `ξ` even, `⟨1,ξ⟩=1`); the
only difference is `T ⪰ 0`. In both, `D''`-spectrum `==` zeros of `ξ̂` exactly (the proven det-identity). The
indefinite case realizes those zeros off the real line; the PSD case realizes them on it. So even-simple is not
an incidental technical hypothesis — **it is the entire mechanism by which the spectral realization delivers
real zeros**. The remaining open work is unchanged: proving `QW_λ` actually has this positivity (simple even
ground state) — external/analytic (§ above). This witness only makes the *stakes* of that hypothesis concrete.

_Report-tier (toy 3×3 numeric+symbolic); source THE_BRIDGE.md cy417–419; forge `WeilSpectralWitness.lean`._
