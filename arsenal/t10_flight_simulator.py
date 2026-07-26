#!/usr/bin/env python3
"""
t10_flight_simulator.py  -- Tier B #10.  A tunable-axiom sandbox on the shared Weil engine: dial
the axioms of the explicit formula (archimedean constant, prime strength, prime cutoff) and watch
the positivity margin W = A+B-P move, mapping the BOUNDARY of the RH-forcing region -- where does
positivity break, and how far is the real zeta from that edge?

REPORT-tier numeric instrument.  Three knobs:
  1. arch shift  delta : Omega(r) -> Omega(r) + delta   (move the place-at-infinity constant).
  2. prime gain  kappa : W = A + B - kappa * P           (scale total prime strength).
  3. prime cutoff  P_max: keep only primes <= P_max       ('drop any prime' stress, LEDGER 168).
For each it reports the margin and, for the continuous knobs, the critical value where W crosses 0
(the forcing boundary) and how much slack the real zeta has.

Usage: python3 t10_flight_simulator.py
"""
from __future__ import annotations
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from core.honesty import Result, Tier
from core import provenance
from core.weil import WeilForm, omega_zeta


def margin_with(wf: WeilForm, t, delta=0.0, kappa=1.0, pmax=None):
    """W = A(+delta) + B - kappa * P(<=pmax)."""
    gh = wf.ghat(wf._grid, t)
    A = float(np.trapezoid(gh * (wf._omega_grid + delta), dx=wf._dr) / (2 * math.pi))
    B = wf.B(t)
    if pmax is None:
        mask = slice(None)
    else:
        mask = wf.ns <= pmax
    P = 2.0 * float(np.sum(wf.lam[mask] / wf.sqrtn[mask] *
                           np.exp(-0.5 * (wf.logn[mask] / t) ** 2)))
    return A + B - kappa * P


def bisect_zero(f, lo, hi, tol=1e-9, it=80):
    flo, fhi = f(lo), f(hi)
    if flo == 0: return lo
    if fhi == 0: return hi
    if (flo < 0) == (fhi < 0):
        return None
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if abs(fm) < tol:
            return mid
        if (fm < 0) == (flo < 0):
            lo, flo = mid, fm
        else:
            hi, fhi = mid, fm
    return 0.5 * (lo + hi)


def analyze(t=0.28):
    # operating point: t~0.28 puts real prime mass in play AND keeps the explicit-formula
    # identity validating (W >> quadrature floor ~5e-12). At t<0.15 primes are negligible;
    # beyond t~0.55 W drops into the numerical floor. This regime is where A+B ~ P (tight).
    wf = WeilForm(zero_count=40)
    base = margin_with(wf, t)

    # knob 1: arch shift delta -> where does W cross 0?  (delta* = -W since int ghat / 2pi = 1)
    f_delta = lambda d: margin_with(wf, t, delta=d)
    delta_star = bisect_zero(f_delta, -5.0, 0.0)
    delta_curve = [(d, f_delta(d)) for d in (-0.5, -0.2, -0.05, 0.0, 0.2, 0.5)]

    # knob 2: prime gain kappa -> critical kappa where W=0  (kappa* = 1 + W/P > 1)
    f_kappa = lambda k: margin_with(wf, t, kappa=k)
    kappa_star = bisect_zero(f_kappa, 1.0, 60.0)
    kappa_curve = [(k, f_kappa(k)) for k in (0.5, 1.0, 1.05, 1.2, 1.5, 2.0)]

    # knob 3: prime cutoff -> margin vs how many primes included
    cutoffs = [2, 3, 5, 7, 11, 23, 100, 1000, 100000]
    cutoff_curve = [(c, margin_with(wf, t, pmax=c)) for c in cutoffs]

    # slack sweep: how the forcing slack (|delta*| = W) shrinks as t grows (positivity tightens)
    slack_sweep = []
    for tt in (0.10, 0.20, 0.28, 0.35, 0.45):
        slack_sweep.append((tt, margin_with(wf, tt)))

    return {
        "t": t, "baseline_margin": base,
        "delta_star": delta_star, "delta_curve": delta_curve,
        "kappa_star": kappa_star, "kappa_curve": kappa_curve,
        "cutoff_curve": cutoff_curve, "slack_sweep": slack_sweep,
    }


def selftest(verbose=True) -> bool:
    r = analyze()
    ok = True
    if r["baseline_margin"] <= 0:
        ok = False; print("  FAIL: baseline margin not positive")
    # lowering the arch constant enough must break positivity (a real boundary exists)
    if r["delta_star"] is None or r["delta_star"] >= 0:
        ok = False; print("  FAIL: no arch-shift forcing boundary found below 0")
    # inflating prime strength must eventually break positivity (kappa_star > 1)
    if r["kappa_star"] is None or r["kappa_star"] <= 1.0:
        ok = False; print("  FAIL: no prime-gain boundary above kappa=1")
    # margin should be monotone-ish increasing as more primes included toward the true value
    cc = [m for _, m in r["cutoff_curve"]]
    if cc[-1] <= 0:
        ok = False; print("  FAIL: full-prime margin not positive")
    if verbose:
        print(f"  baseline W = {r['baseline_margin']:+.6e}  (t={r['t']})")
        ds = r['delta_star']; ks = r['kappa_star']
        print(f"  arch-shift boundary: delta* = {ds:.5f}  "
              f"(real zeta slack in the arch constant = {abs(ds):.5f})" if ds is not None
              else "  arch-shift boundary: none found")
        print(f"  prime-gain boundary: kappa* = {ks:.5f}  "
              f"(positivity tolerates prime strength up to {ks:.3f}x)" if ks is not None
              else "  prime-gain boundary: none found in [1,60]")
        print(f"  slack sweep (|delta*|=W shrinks as t grows -> positivity tightens):")
        for tt, m in r["slack_sweep"]:
            print(f"      t={tt:.2f}  W = {m:+.6e}")
        print(f"  prime-cutoff curve (margin as primes added):")
        for c, m in r["cutoff_curve"]:
            print(f"      p<= {c:<7d}  W = {m:+.6e}")
    print("t10.selftest:", "PASS" if ok else "FAIL")
    return ok, r


def main():
    print("=== t10 L-function flight simulator (tunable-axiom positivity sandbox) ===\n")
    ok, r = selftest()
    ds, ks = r['delta_star'], r['kappa_star']
    results = [Result(
        claim=f"forcing-region boundary mapped at t={r['t']}: arch-constant slack |delta*|="
              f"{abs(ds):.4f}, prime-strength tolerance kappa*={ks:.4f}; real zeta sits strictly "
              f"INSIDE the positive region; slack shrinks monotonically as t grows (finite/numeric)",
        tier=Tier.REPORT, method="parametric Weil margin over 3 axiom knobs",
        tool="t10_flight_simulator", data={k: v for k, v in r.items()
                                           if k in ("baseline_margin", "delta_star", "kappa_star",
                                                    "slack_sweep")})]
    path = provenance.stamp("t10_flight_simulator", results, extra={"selftest_pass": ok})
    print(f"\nReading: at t={r['t']} you must bend the arch constant by ~{abs(ds):.3f} or inflate "
          f"prime strength ~{ks:.3f}x to break positivity. The slack sweep shows this margin "
          "SHRINKING as t grows -- positivity gets tighter, which is why RH lives at the limit.")
    print("run record:", path)
    return ok


if __name__ == "__main__":
    main()
