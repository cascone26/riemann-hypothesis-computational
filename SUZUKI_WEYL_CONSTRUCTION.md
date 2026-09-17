# The Suzuki–Weyl m-function construction (LEDGER §982–985, 2026-08-12)

**One-line:** the pole-count = κ lever that §980–981 named as the single non-circular route to RH — and
called *genuinely unbuilt* — now composes correctly on a real object: the correct Weyl m-function built from
Suzuki's canonical system is real-symmetric, reads κ=0 in the built regime, and faithfully detects off-line
zeros. The apparent obstructions ("κ is fragile", "naive m counts zeros not the pole") were **measurement
artifacts of a wrong construction**, now diagnosed and fixed. This does **not** prove RH — the a→∞ limit
(Suzuki Cor 1.6) remains the open analytic step — but it removes false walls and de-risks the lever.

## The frontier before this arc (§980–981)
The Weil form Q = arch + pole − prime is Π₁ iff RH (κ = #negative squares ≤ 1, absorbed by ζ's s=1 pole).
The only non-circular lever identified: reformulate as a generalized-Nevanlinna m-function and count its UHP
poles (Krein–Langer: κ = pole-count). Blockers on record:
- **"naive Weyl ratio m = v̂₊/v̂₋ counts in-window zeros (~10–12), not the pole"** (§974/§978).
- **"κ is fragile — jumps 1→5→28 under ±0.1% Λ perturbation"** (kappa_fragility_analysis.py).
Both read as reasons the lever might be hopeless. Both turned out to be artifacts.

## §982 — the instrument was broken (function-field validation)
Tested the pole-count=κ counter in the **proven** function-field case (RH_ff is a theorem; zeros controlled).
The project's own §957 harness gave a meaningless control (κ=8 where κ=0 required). Two bugs:
1. **Sign.** m(z)=Σ 1/(z−uᵢ) with residues +1 is *anti*-Herglotz (Im m<0 in UHP) → kernel negative-definite →
   counter returns full rank. Fix: proper Herglotz N(z)=Σ res/(uᵢ−z).
2. **Structure.** A single non-real pole 1/(β−z) is **not** a generalized-Nevanlinna function (fails
   N(z̄)=N(z)‾). Off-circle Frobenius zeros must enter as **real-symmetric conjugate pairs** {u, ū} — exactly
   what the functional equation forces (α at r>√q ⟺ partner q/α at r<√q). Each pair = one negative square.

Corrected file `operators/ff_nevanlinna_CORRECTED.py`: control κ=0; sweep κ = off-circle count **exactly**
(0,1,2,3,4,5); κ **invariant under ±1% perturbation** (40 trials). → "κ fragility" retired as a broken-counter
artifact; §981's "robust integer, not a fragile margin" **confirmed**. Resolution floor (honest): near-critical-
line poles undercount unless test-points sampled near the axis (y↓0.05) — a numerical limit mapping to the real
difficulty of resolving near-line zeros, not a mechanism defect.

## §983 — ζ-side companion (real ζ zeros)
`operators/zeta_kappa_mechanism.py`. ζ's s=1 pole lies **off** the critical line, so in height-coordinates
(ζ(½+iz)) it is a conjugate pair at ±i/2 = one negative square; on-line zeros are real heights = zero negative
squares. Result: on-line + pole → κ=1 (robust ±0.5% jitter); each off-line pair → +1 (κ=2,3,4). So **RH ⟺ κ=1**
confirmed on ζ's actual data with the corrected instrument.

## §984 — the correct construction from Suzuki's operator (the "unbuilt" step)
`operators/suzuki_weyl_kappa.py`. From Suzuki's W_family (v̂₊, v̂₋; characteristic function W(a,θ;z), Thm 1.5
zeros all real), form the **Weyl inner function** Θ(z) = −(z−i)v̂₊/((z+i)v̂₋) and the Weyl m-function as its
**Cayley transform** m(z) = i(1+Θ)/(1−Θ). Controls:
- **|Θ|→1 as Im z→0** (0.42→0.81→0.98) — genuine inner function. ✓
- **Real-symmetry** (the crux): naive v̂₊/v̂₋ defect **3.5–5.2** (not generalized-Nevanlinna — the §982 category
  error, measured) vs Weyl–Cayley m defect **4e-16** (real-symmetric to machine precision). ✓
- **κ (corrected counter):** naive → 20–23 (garbage, "counts zeros" — §974/§978, now *explained*); Weyl–Cayley
  m → **0** (kernel min-eig ≈ −5e-15 = numerical zero) — correct: Suzuki's ω>1 operator is self-adjoint.
- **Faithfulness:** injecting 0/1/2/3 off-line conjugate pairs → κ = 0/1/2/3 exactly; base κ=0 stable across
  lam_shift. The construction genuinely detects off-line structure.

## §985 — convergence toward ζ
`operators/suzuki_a_scaling.py`. Across support a=1.2→3.0: correct-m κ stays **0** (kernel PSD) while the
W-zeros lock onto the true ζ-zeros with monotonically shrinking error (0.043→0.0095). The correct construction
behaves exactly as the theory demands as it approaches ζ.

## §986 — the density-wall test (confirms the field's wall, does not cross it)
`operators/suzuki_density_wall.py`. Does the correct construction escape the density-mismatch wall
(Kim 2607.24830, "fundamental")? Best-θ W-zero match vs ζ-zeros, window scaled with a:

| a | zmax | #ζ-zeros in window | #W matched (<0.15) | RvM N(zmax) |
|---|------|--------------------|--------------------|-------------|
|1.5|42.5| 7 | 7 | 7.0 |
|2.0|50.0|10 |10 | 9.4 |
|2.5|57.5|12 |12 |12.0 |
|3.0|65.0|14 |14 |14.7 |

The construction matches **every** ζ-zero up to the height its support a can resolve — matched count tracks
RvM N(T) exactly when a is scaled to the window. This is the density relation (W-zero density a/π ↔ ζ density
(1/2π)log(T/2π)) behaving precisely as the field understands it: finite a resolves finitely many zeros; the full
set needs a→∞. So §986 **confirms** the construction is faithful up to its density-limited height and respects
the a/π↔log-density law — it does not cross the wall. Not a breakthrough; a correct characterization of it.

## §987 — the winding-number (argument-principle) κ handle
`operators/suzuki_winding.py`. The κ of the correct construction now has a second, **topological** realization:
the winding number of (Θ−1) around a UHP contour. Results: winding[Θ−1] = κ = 0 at every a (matches the kernel
count), computed from Θ **alone** (no assumption on zero locations); cleanly separates winding[Θ] = zero-count/
density (~a/π·2R) from winding[Θ−1] = pole-count = κ (the §974/§978 distinction, resolved topologically); and
faithful — injecting k UHP zeros → winding = k. **Conceptual gain:** κ is a winding number = a robust integer
stable under continuous deformation, so the a→∞ wall reframes as a **continuity condition** ("no zero of (Θ−1)
crosses into the UHP as a→∞") — a topological invariant rather than the raw density mismatch. Not a proof; a
sharper handle on the same open limit.

## Honest ceiling & standing frontier
Suzuki Thm 1.5 makes W's zeros real for **all finite a** unconditionally → at finite a the correct m
necessarily reads κ=0, carrying **no RH content**. RH lives only in the **a→∞ (ω→1) limit** (Cor 1.6), where the
truncation vanishes and ζ's true structure (possible off-line zeros) would appear. **The one open analytic step
is now sharply posed on a correct, verified object:** control κ of the correct Weyl m in the a→∞ limit. Nothing
here proves RH; what changed is that the lever is now built, faithful, and free of the false obstructions that
made it look hopeless.

## Files
`operators/ff_nevanlinna_CORRECTED.py` · `operators/zeta_kappa_mechanism.py` ·
`operators/suzuki_weyl_kappa.py` · `operators/suzuki_a_scaling.py`. LEDGER §982–985.
