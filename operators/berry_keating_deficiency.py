#!/usr/bin/env python3
"""
berry_keating_deficiency.py — RESEARCH_QUEUE #6 (2026-07-20). The convergence point of §37/§38/§39.

§39 (graph-zeta) showed: in every SOLVED-RH analog, positivity is FREE because a bounded self-adjoint
operator supplies a real spectrum. For ζ that operator is the missing Hilbert–Pólya operator, whose
classical limit is Berry–Keating H = xp. The precise unfinished question (Bender–Brody–Müller 2017) is
SELF-ADJOINTNESS, which is decided by DEFICIENCY INDICES n_± = dim ker(H* ∓ i):
  • n_+ = n_- = 0  ⟹ essentially self-adjoint: UNIQUE extension (but xp on (0,∞) then has CONTINUOUS
    spectrum — no discrete zeros).
  • n_+ = n_- = k > 0 ⟹ a k²-parameter FAMILY of self-adjoint extensions, each with (generically)
    DISCRETE spectrum — ROOM for a boundary condition to encode the primes.
  • n_+ ≠ n_-  ⟹ NO self-adjoint extension exists.

The existing repo script `auto_berry_keating_xp.py` FORCES Hermiticity ((H+H†)/2), sidestepping this
entirely; `operator_candidates.md` only CITES the indices. This script actually COMPUTES them (analytic
solution of H*ψ=±iψ + L²-normalizability at each endpoint) for the physically-motivated xp variants,
and — for the case with a family of extensions — computes the extension spectrum vs the boundary phase.

Payoff: locates PRECISELY where "the primes must enter" in the Hilbert–Pólya program — as the boundary
condition selecting one self-adjoint extension — making §39's mechanism diff concrete and computable.
"""

import numpy as np
import json

# --------------------------------------------------------------------------------------
# The symmetrized dilation operator D = (xp+px)/2 = -i(x d/dx + 1/2),  p = -i d/dx.
# H*ψ = ±i ψ  ⟺  -i(x ψ' + ψ/2) = ±i ψ  ⟺  x ψ' = -(1/2 ± 1) ψ  ⟺  ψ = x^{-(1/2±1)}.
#   +i:  ψ_+ = x^{-3/2}
#   -i:  ψ_- = x^{+1/2}
# Deficiency index n_± = 1 if ψ_± ∈ L²(domain), else 0.
# --------------------------------------------------------------------------------------

def is_L2_powerlaw(power, a, b):
    """Is ψ = x^power in L²(a,b)?  ∫|x^power|² dx = ∫ x^{2·power} dx.  ANALYTIC endpoint test
    (numerical quad with a finite cutoff cannot detect endpoint divergence — that was the bug):
      • near 0 (a=0): converges ⟺ 2·power > -1 ⟺ power > -1/2
      • near ∞ (b=∞): converges ⟺ 2·power < -1 ⟺ power < -1/2
      • a bounded closed subinterval (0<a<b<∞) of a continuous power: always L².
    """
    ok_at_a = (power > -0.5) if a == 0 else True         # only 0 is a singular endpoint here
    ok_at_b = (power < -0.5) if np.isinf(b) else True     # only ∞ is a singular endpoint here
    return ok_at_a and ok_at_b


def deficiency_dilation(a, b, label):
    """Deficiency indices of D=(xp+px)/2 on L²(a,b). ψ_+=x^{-3/2}, ψ_-=x^{1/2}."""
    n_plus = int(is_L2_powerlaw(-1.5, a, b))
    n_minus = int(is_L2_powerlaw(0.5, a, b))
    return n_plus, n_minus, label


def analytic_check():
    print("=" * 78)
    print("DEFICIENCY INDICES of the symmetrized dilation D=(xp+px)/2 = -i(x d/dx + 1/2)")
    print("  H*ψ=+iψ ⟹ ψ_+ = x^{-3/2};   H*ψ=-iψ ⟹ ψ_- = x^{+1/2}")
    print("=" * 78)
    cases = [
        (0.0, np.inf, "(0,∞) half-line  [Berry–Keating bare]"),
        (1.0, 60.0, "[1,60] finite interval  [B–K phase-space truncation]"),
        (0.0, 1.0, "(0,1] — singular at 0 only"),
        (1.0, np.inf, "[1,∞) — singular at ∞ only"),
    ]
    results = {}
    for a, b, label in cases:
        n_p, n_m, _ = deficiency_dilation(a, b, label)
        if n_p == n_m == 0:
            verdict = "essentially self-adjoint (UNIQUE ext; but xp ⟹ CONTINUOUS spectrum)"
        elif n_p == n_m:
            verdict = f"({n_p}²={n_p*n_p})-param FAMILY of self-adjoint extensions (DISCRETE spectra — room for primes)"
        else:
            verdict = "NO self-adjoint extension (n_+≠n_-)"
        print(f"\n  {label}")
        print(f"    n_+ (x^-3/2 ∈ L²?) = {n_p} ,  n_- (x^1/2 ∈ L²?) = {n_m}   ⟹ {verdict}")
        results[label] = {"n_plus": n_p, "n_minus": n_m, "verdict": verdict}
    return results


# --------------------------------------------------------------------------------------
# For the finite-interval case (deficiency (1,1)): the self-adjoint extensions of D on [a,b]
# are parametrized by a U(1) boundary condition ψ(b) = e^{iφ} (b/a)^{...} ψ(a) coupling the two
# endpoints. Compute the spectrum as a function of φ. D generates dilations; on log-coordinate
# t=log x, D = -i(d/dt) (translation!) on [log a, log b], i.e. momentum on a CIRCLE of
# circumference Λ=log(b/a) with boundary phase φ. Its spectrum is EXACTLY
#     E_n(φ) = (2π n + φ)/Λ ,  n ∈ ℤ   — evenly spaced, gap 2π/Λ.
# This is the KEY structural fact: the truncated Berry–Keating operator is momentum-on-a-circle;
# every self-adjoint extension has a RIGID arithmetic-progression spectrum (gap set by the
# phase-space size Λ=log(b/a)), NOT the irregular ζ zeros. The primes cannot be encoded by a
# single boundary phase — they would need the boundary condition to be NONLOCAL / energy-dependent.
# --------------------------------------------------------------------------------------

def extension_spectrum(a, b, phi, nmax=8):
    """Spectrum of the self-adjoint extension of D on [a,b] with boundary phase phi:
    in log-coordinate D = -i d/dt on a circle of circumference Λ=log(b/a). E_n=(2πn+phi)/Λ."""
    Lam = np.log(b / a)
    return np.array([(2 * np.pi * n + phi) / Lam for n in range(-nmax, nmax + 1)])


def spectrum_analysis():
    print("\n" + "=" * 78)
    print("EXTENSION SPECTRUM of the truncated B–K operator (deficiency (1,1) case)")
    print("  In t=log x coords, D = -i d/dt = momentum on a CIRCLE (circumference Λ=log(b/a)).")
    print("=" * 78)
    a, b = 1.0, 60.0
    Lam = np.log(b / a)
    print(f"  interval [1,60]: Λ = log(60) = {Lam:.4f}, spectral gap 2π/Λ = {2*np.pi/Lam:.4f} (RIGID)")
    for phi in [0.0, np.pi / 2, np.pi]:
        sp = extension_spectrum(a, b, phi, nmax=4)
        print(f"    φ={phi:.3f}: E_n = {np.round(sp, 3)}")
    # compare to actual ζ zero gaps (irregular)
    zeros = [14.1347, 21.0220, 25.0109, 30.4249, 32.9351, 37.5862, 40.9187]
    zgaps = np.diff(zeros)
    print(f"\n  ζ zero gaps (first few): {np.round(zgaps, 3)}  — mean {zgaps.mean():.3f}, "
          f"std {zgaps.std():.3f} (IRREGULAR)")
    print(f"  B–K extension gaps: CONSTANT {2*np.pi/Lam:.3f} (std 0).")
    print("  ⟹ NO single boundary phase reproduces the irregular ζ spacings. The primes cannot")
    print("     enter as a local U(1) boundary phase; they require a NONLOCAL / energy-dependent")
    print("     boundary condition — precisely the oscillatory prime/periodic-orbit correction that")
    print("     `auto_berry_keating_xp.py` isolated as 'missing' on top of the smooth Weyl term.")
    return {"Lambda": float(Lam), "bk_gap": float(2 * np.pi / Lam),
            "zeta_gap_mean": float(zgaps.mean()), "zeta_gap_std": float(zgaps.std())}


def main():
    r1 = analytic_check()
    r2 = spectrum_analysis()

    print("\n" + "=" * 78)
    print("SYNTHESIS (ties §37/§38/§39 → the missing operator, now LOCALIZED):")
    print("=" * 78)
    print("""
  • Bare xp on (0,∞): essentially self-adjoint (n_±=0) but CONTINUOUS spectrum — no zeros. Dead end,
    exactly as Berry–Keating/Endres–Steiner found.
  • Truncated xp on [a,b]: deficiency (1,1) ⟹ a U(1) family of self-adjoint extensions WITH discrete
    spectra — but each is momentum-on-a-circle: a RIGID arithmetic progression (gap 2π/log(b/a)),
    never the irregular ζ zeros. A local boundary phase is NOT enough.
  • CONCLUSION (the localized missing ingredient): the Hilbert–Pólya operator cannot be a self-adjoint
    EXTENSION of truncated xp via a local boundary condition. The primes must enter as a NONLOCAL,
    energy-dependent boundary condition (equivalently, the periodic-orbit sum). This is the SAME
    "prime/oscillatory correction on top of the smooth self-adjoint spectrum" that appears in every
    formulation — §39's 'self-adjointness carries positivity' is true, but the self-adjoint operator
    that has the RIGHT (irregular, zero-matching) spectrum is exactly the nonlocal object we lack.
  • This SHARPENS the target: not 'find any self-adjoint operator' (truncated xp already is one) but
    'find the self-adjoint operator whose boundary/nonlocal data IS the prime periodic-orbit sum' —
    unifying the Berry–Keating (operator) and explicit-formula (primes) sides into one object.
""")
    out = {"deficiency": r1, "spectrum": r2}
    with open("results/berry_keating_deficiency.json", "w") as f:
        json.dump(out, f, indent=2)
    print("Results: results/berry_keating_deficiency.json")


if __name__ == "__main__":
    main()
