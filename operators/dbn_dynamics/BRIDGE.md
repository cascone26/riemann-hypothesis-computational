# BRIDGE.md — The Suzuki-screw ↔ Weil ↔ dBN-flow ↔ log-gas ↔ prime-Ξ bridge

Consolidated theorem-shaped writeup of the dynamical arc (LEDGER §85–92), grounded in a
term-by-term comparison with the actual Rodgers–Tao source (arXiv:1801.05914). Written
2026-07-20 (riemann4, Opus) after the §92 verify-before-ship read.

**Honesty banner (Charisma Report/Handle):** the log-gas half of this bridge is RODGERS–TAO's,
not ours — see §1. Our genuine, RT-absent contribution is the three-item delta in §3. Every
"verified" below points at a tool in `operators/dbn_dynamics/`. The RH-closing step is OPEN and
its honest shape is §4.

---

## 0. The objects (one classical picture)

- Pólya/dBN: `Φ(u) = Σ_n (2π²n⁴e^{9u} − 3πn²e^{5u})e^{−πn²e^{4u}}`, even, super-exp decay.
  `H_t(z) = ∫₀^∞ e^{tu²}Φ(u)cos(zu)du`, backward heat `∂_t H_t = −∂_{zz}H_t`. `H_0(z)=⅛ξ(½+iz/2)`.
  RH ⟺ all zeros of `H_0` real ⟺ **Λ ≤ 0** (de Bruijn–Newman constant). [RT eqs 78–103; our §65.]
- Zeros `x_j(t)` of `H_t` obey the repulsion ODE `∂_t x_k = 2 Σ_{j≠k} 1/(x_k − x_j)`. [RT eq (ode) 711; CSV.]
- Suzuki screw function `g(t)` = our Weil `A+B−P`; localized form `λ_a = inf Q(v)/‖v‖²` on `(−a,a)`,
  RH ⟺ `λ_a ≥ 0 ∀a`. [Suzuki 2606.09096; our `operators/suzuki_screw_lambda_a.py`.]

These are the SAME `ξ`. The bridge is that the dBN flow of the zeros and the screw/Weil positivity
of the primes are two faces of one object.

---

## 1. The log-gas half — THIS IS RODGERS–TAO (full credit, NOT ours)

The `§92` read established that the following are all explicit in RT Section 1/4. We rediscovered
them independently (§86–91); they are **not** novel. Recorded here so the bridge is honest.

| Our label | Claim | RT location | Status |
|-----------|-------|-------------|--------|
| §86 | flow = gradient flow of log-gas potential | RT 723–726, `H=Σ_{j≠k}log 1/|x_k−x_j|`, `∂_t x_k=−∂_{x_k}H` | RT's. Our `E=Σlog|·|` = `−H`. |
| §87 | log-energy concave (`Hess=−L`, `L` PSD) | RT 727/154 "the Hamiltonian is convex" | RT's. `H` convex ⟺ our `E` concave. |
| §90 | `d/ds Σx_j² = 4N(2N−1)`, `Σx_j` conserved | RT Lemma Identities (iv), Virial identity, line 763 | RT's virial identity. |
| §90/§91 | Weil/energy derivative = definite local part + indefinite cross part | RT line 771 + Lemma (iii)(v) | RT's; controlling cross terms = RT §6–9. |
| §91 | flow energy = 2-point(+)/3-point(−) near-cancellation | RT line 771 "monotone if cross terms controlled" | RT's known crux. |

**Verified-in-repo (our numerics reproducing RT's structure):** `E_concavity.py` (Hess E = −L to
2e-9, L PSD); `weil_flow_deriv.py` (virial constant 4N(2N−1)=319200 to the digit); `flow_energy_
structure.py` (2-pt +1681 / 3-pt −1179 split). These CONFIRM RT, they don't extend RT.

**RT's theorem and its DIRECTION.** RT prove **Λ ≥ 0** (Theorem 1). Mechanism (RT 132–135, endgame
1571): *assume Λ<0 ⇒ backward heat flow relaxes zeros to a local arithmetic progression ⇒ contradicts
CGGGH/Montgomery pair correlation (a positive proportion of gaps are below mean spacing).* This is a
**lower bound only** — the method structurally cannot produce an upper bound on Λ.

---

## 2. Why RT does not touch RH

RH = **Λ ≤ 0**, an UPPER bound. RT's engine bounds Λ from below by contradiction with known zero
statistics. There is no version of "relax-and-contradict-pair-correlation" that yields Λ≤0: the
forward flow always disperses toward equilibrium, and equilibrium is *consistent* with Λ≤0, not
contradicted by it. So RT's dynamical machinery, as wielded, is a Λ≥0 instrument. **RH needs an
ingredient RT's framework structurally lacks.**

---

## 3. The delta — what is genuinely ours and RT-absent

1. **A POSITIVITY functional on the flow: Suzuki's `λ_a ≥ 0` (the screw/Weil form).**
   RT use NO positive-definite object anywhere — only heat relaxation + external statistics. The screw
   kernel `g` is a *positive screw increment* (`g'' ≥ 0` as a measure on the atomic part) exactly when
   RH holds, and the modular/theta structure of `Φ` that could force this is present **only at λ=0**
   (heat flow for λ≠0 kills the Euler product; §47d-e). This is a genuinely different handle, pointed
   at the correct (upper-bound) direction.
   Tools: `suzuki_screw_lambda_a.py`, `screw_flow.py`.
2. **Prime-explicit velocity `v_j = Ξ″(γ_j)/Ξ′(γ_j)` (§88).** The flow velocity of each zero is a
   prime-determined quantity via the Hadamard/Euler duality. RT connect the flow to arithmetic ONLY
   through pair-correlation inputs, never through the explicit formula. Honest status: correct but
   likely folklore (direct from Hadamard) — Report-tier until it does load-bearing work.
   Tool: `velocity_bridge.py` (verified ~2–6% truncation error).
3. **Relaxation spectrum split (§89).** The log-gas Laplacian `L` at the zeta configuration splits
   into a rigid-long sector (slowest eigenmode 0.04868 = equal-spacing exactly) and a GUE-short sector
   (fast close-pair modes). Λ localizes to the fast/close-pair modes. This is a diagnostic RT don't
   compute; it says WHERE in the spectrum the Λ-obstruction lives.
   Tool: `relaxation_spectrum.py`.

---

## 4. The sharpened frontier (OPEN — honest shape)

The live question is **not** "find a Lyapunov functional" — RT have it. It is:

> **Does the screw/Weil positivity `λ_a ≥ 0` (which RT lack) constrain the flow in the Λ≤0
> direction — i.e. does positivity of the screw kernel at `t=0` obstruct the zeros from having
> flowed up from a configuration with an off-line zero (`Λ>0`)?**

Concretely: RT flow FORWARD to force equilibrium-and-contradiction for `Λ<0`. The RH side must run
the argument the other way — positivity of `g` at `t=0` (the only time the modular structure is
present) as the obstruction to `Λ` being positive. The pair/triple-correlation cross terms that RT
control statistically (§6–9) are the same cross terms that appear in our screw-form derivative (§91);
the bet is that **positivity gives a DETERMINISTIC handle on those cross terms where RT use
statistics** — the individual-isolating ingredient (§59/§62) the averaging routes lacked.

### 4a. Probe results (§93 — EXECUTED, hypothesis ELIMINATED)
Both probes ran (`screw_vs_crossterm.py`, `offline_crossterm.py`). Outcome = a clean NULL that narrows the map:
- **Probe #1 (backward-flow)** is not numerically observable: §84/§89 already showed Λ is a high-zero/close-pair
  phenomenon; low zeros stay real both flow directions, so the positivity-margin transition is invisible at
  computable height. (Honest, not a failure — it's why the tail is hard.)
- **Probe #2 (cross-term control):** on real zeros `|CROSS_K| < |LOCAL_K|` in ALL windows (0.12–0.93) — RT's
  indefinite cross-term is already dominated by the definite sum-of-squares term, so energy monotonicity looks
  deterministic. **But #2b (decisive):** pinching a pair off-line does NOT push the ratio past 1 — it plunges
  (0.12→0.0002), because a close complex pair makes the local term blow up like 1/δ⁴. So **`|CROSS|<|LOCAL|` is
  RH-BLIND**: RT's energy monotonicity holds for off-line configs too, carrying no reality information.

**⇒ Conclusion:** the screw POSITIVITY (Λ≤0/RH content) does NOT couple to the dBN flow through the energy
monotonicity — that structure is RH-blind. The dynamical bridge is DIAGNOSTIC, not the crack. The RH content sits
in individual-reality (screw `λ_a≥0`), orthogonal to the flow's lower-bound (Λ≥0) machinery. **Do not re-enter the
flow-energy route for Λ≤0.** The relocated spearhead (LEDGER §93) is the individual-isolation mechanism: Connes'
finite-place→global Euler-product convergence (§59/§73), the one framework that isolates the individual curve/ζ.

---

## Provenance
LEDGER §85–92. RT source: `membrane_pulls/rt_1801/fmp-template.tex`. Tools: this directory (8 files).
The §91 headline ("E = the Lyapunov functional RT lacked") is **retracted** by §92; corrected here.
Companion: memory `rh-suzuki-screw-function-frontier`, `rh-positivity-global-metatheorem`.
