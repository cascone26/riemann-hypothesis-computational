# RH_MAP.md — the campaign's durable map (riemann4, banked 2026-07-21, LEDGER through §129)

**Charisma labeling:** Report = researched/described/computed-by-us; Handle = a theorem actually proven (by us or,
where noted, in the literature and verified-against-source by us). Doc offline 2026-07-21 → full current-schema
Charisma stocking deferred; stable Report/Handle core applied. This map is the deliverable, NOT a proof.

## 1. THE ONE WALL, SIX FACES (operator-positivity basin — exhausted)
RH ⟺ Weil positivity, and it has (at least) six proven-equivalent faces. Each is UNCONDITIONAL in a regime, STUCK at
the same limit (support a→∞ / parameter ω→0 / λ→∞ / all-primes). No face is a shortcut past the limit (verified §121).

| Face | Handle (proven) | Report (our finding) | Stuck at |
|---|---|---|---|
| k_λ ≈ θ_x (CCM 2511.22755 §7) | even-simple⟹on-line (Connes) | our §94-103 reconstruction | a→∞ approximation |
| self-adjoint extension selection | Thm 1.5 reality+simplicity (Suzuki) | §113-116 (Herglotz claim RETRACTED §117) | limit-circle b.c. at ∞ |
| E_ξ=Ξ+iΞ′ Hermite–Biehler / Θ_ξ inner | classical | §117 (m_a is INNER not Herglotz) | RH itself |
| Jacobi / orthogonal-poly positivity (CM) | — | §119 semilocal Jacobi (CM's deferred computation) | all-primes |
| canonical-system Hamiltonian H(t)≥0 (Suzuki 1204) | ω>1 construction | — | ω→0 |
| Parseval / Gram factorization G=U*U | — | §121 (S_x carries zeros via Θ_ξ) | the zeros in S_x |

## 2. CROSS-BASIN VERDICT (every other approach, verified §123-125)
- ZERO-DENSITY estimates (Guth-Maynard 2405.20552, 2024): NOT an RH-path (N(σ,T)>0 always). Handle: it's a real advance,
  wrong target.
- POSITIVE-PROPORTION on line (Levinson-Conrey): the RH-branch, stuck ~41.7% (Pratt 2020); bottleneck moved off
  mollifier-length to linear-combo optimization (Conrey 2508.11108, 2025).
- LeClair spectral-flow (2406.01828): = KNOWN RH-equivalent (Lagarias −Re(ζ'/ζ)≲½log(t/2π)), physics-dressed, unitarity-
  at-σ=½ IS RH. Cross-connection: Lagarias pointwise ⟷ our §81 band-limited |v̂|² are Fourier duals.
- Hilbert-Pólya/Berry-Keating: density-mismatch, no candidate operator. dBN Λ≤0: dormant. Function-field/Connes-site,
  moments/RMT, Beurling-Nyman: SAME WALL relabeled.
- Meta-state (verified, not assumed): every currently-known RH approach is stuck at a genuine obstruction.

## 3. NOVEL DEPOSITS worth keeping (Report-tier unless noted)
- §81 Q_prime = (1/π)∫Re[ζ'/ζ(½+iξ)]|v̂|²dξ multiplier formula (Handle: verified, ratio 1.0000).
- §112 prime-bend: primes inject zeta-like gap-irregularity into the localized spectrum (CV(full)=0.393≈zeta 0.389).
- §119 semilocal Jacobi matrix (CM 2310.18423's DEFERRED computation, done + validated): primes suppress low-index
  recurrence coefficients; all-primes limit DEGENERATES (§120, measure concentrates at s=0).
- §122 arch↔prime cancellation: arch O(−e^a) ↔ prime +O(e^a), λ_a>0 clean to a=5 (Handle: computed).
- §129 off-line-zero DAMAGE LAW (Lens G, rule-13-confirmed novel intersection): damage ~ δ²·f(a)·[γ₀-FLAT]. ⟹
  density⟹effective-Weil-bound is USELESS (obstructions compound); detector redundant with zero-verification.

## 4. TOOLS (operators/dbn_dynamics/ + operators/)
compact G_a (g_a_compact.py — breaks the large-a conditioning wall, §110); mpmath hi-prec Weil (mpmath_weil_hiprec.py,
§109); suzuki_screw_lambda_a.py (core λ_a); suzuki_W_zeros.py; semilocal_jacobi.py; lens_g_execute.py; prime-bend;
Weyl m-function (weyl_mfunction.py). Sources: membrane_pulls/ (suzuki 2606/2301, connes 2602, RT 1801, CM 2310, GM 2405,
LeClair 2406). Data: 100k zeros, high-precision.

## 5. SELF-CORRECTIONS (the discipline that held)
- §92: "E = Lyapunov functional RT lacked" RETRACTED — E is exactly Rodgers-Tao's Hamiltonian.
- §117: "m_a Herglotz, reality-for-free" RETRACTED — m_a is unimodular/inner (|v̂_+|=|v̂_-| exact). Lesson: check a
  property's DEFINITIONAL one-liner (|m|≡1) before the interpretive overlay.
- §127: Lens-G float64 attempt cancellation-corrupted → §129 zero-side reformulation fixed it.

## 6. TIER-3 FORWARD (the real hunt — see NEXT_SESSION_RH.md / the strategic note)
The function-field analogy (Weil-Deligne template — the ONE place RH-type statements are PROVEN) + the missing Λ≤0
monotone functional (Perelman/Rodgers-Tao template). Worked at the object level, not just labeled. Not started.
