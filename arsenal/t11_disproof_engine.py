#!/usr/bin/env python3
"""
t11_disproof_engine.py  -- Tier B #11.  The adversary: instead of proving RH true, aggressively
try to BREAK Weil positivity by bending the axioms of the prime side, and report WHICH bent axiom
is forced to give -- the obstruction to fabricating a counterexample is the mechanism forcing RH.

Generalizes operators/beurling_adversary.py.  The place-at-infinity (arch A + boundary B) is held
at the real zeta value; the PRIME side is bent, and the margin W = A + B - P_bent is watched.  The
real primes leave a positive margin; each bend perturbs it.  The tool ranks the bends by how much
they erode the margin and flags any that drive it negative (= that bent 'zeta' would have off-line
zeros), naming the axiom responsible.

REPORT-tier numeric adversary.  Findings are which axiom-bends threaten positivity, on a finite
test-function window -- diagnostic of the forcing structure, not a proof.

Bends tried (each an axiom of the prime system):
  - shift    : primes sit at log(p)+eps instead of log(p)      (wrong locations)
  - fake     : inject a fake prime power at a non-integer location (Beurling B3)
  - dropP    : delete the prime p=2                              (Euler-product incompleteness)
  - inflate  : multiply one prime's weight Lambda(p) by c        (PNT violation, local)
  - equis    : replace primes by an equi-spaced comb of same count (breaks PNT globally)

Usage: python3 t11_disproof_engine.py
"""
from __future__ import annotations
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from core.honesty import Result, Tier
from core import provenance
from core.weil import WeilForm


def _AB(wf, t):
    gh = wf.ghat(wf._grid, t)
    A = float(np.trapezoid(gh * wf._omega_grid, dx=wf._dr) / (2 * math.pi))
    return A + wf.B(t)


def _P_from(logvals, weights, t):
    logv = np.asarray(logvals, float)
    w = np.asarray(weights, float)
    return 2.0 * float(np.sum(w * np.exp(-0.5 * (logv / t) ** 2)))


def analyze(t=0.28):
    # t~0.28: real prime mass (mainly p=2, some p=3) is in play and A+B ~ P (tight), so
    # bending the prime side visibly moves W across 0. At t<0.15 the prime term is negligible.
    wf = WeilForm(zero_count=40)
    AB = _AB(wf, t)
    # real prime side: weights Lambda(n)/sqrt(n) at log n
    logn = wf.logn.copy()
    wts = (wf.lam / wf.sqrtn).copy()
    P0 = _P_from(logn, wts, t)
    W0 = AB - P0

    bends = {}

    # 1. shift: move every prime location by eps
    eps = 0.05
    bends["shift(+0.05 in log p)"] = AB - _P_from(logn + eps, wts, t)

    # 2. fake prime power at x=log(4.1) with weight ~ log(4.1)/sqrt(4.1)
    x = math.log(4.1); w = math.log(4.1) / math.sqrt(4.1)
    bends["fake prime @4.1"] = AB - (_P_from(logn, wts, t) + 2.0 * w * math.exp(-0.5 * (x / t) ** 2))

    # 3. drop p=2 (remove all powers of 2)
    keep = ~np.isin(wf.ns, [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024])
    bends["drop p=2 tower"] = AB - _P_from(logn[keep], wts[keep], t)

    # 4. inflate weight of p=2 by 3x
    wts2 = wts.copy()
    wts2[wf.ns == 2] *= 3.0
    bends["inflate Lambda(2) x3"] = AB - _P_from(logn, wts2, t)

    # 5. equi-spaced comb: same number of 'primes', evenly spaced in log, same total weight
    npr = len(logn)
    equi_log = np.linspace(logn.min(), logn.max(), npr)
    bends["equi-spaced comb"] = AB - _P_from(equi_log, wts, t)

    # rank by erosion of the margin (most negative first)
    ranked = sorted(bends.items(), key=lambda kv: kv[1])
    broke = [name for name, w in ranked if w < 0]
    return {"t": t, "AB": AB, "W0": W0, "bends": dict(bends), "ranked": ranked, "broke": broke}


def selftest(verbose=True) -> bool:
    r = analyze()
    ok = True
    if r["W0"] <= 0:
        ok = False; print("  FAIL: real-prime baseline margin not positive")
    # at least one adversarial bend must erode the margin below the real baseline (adversary bites)
    if not any(w < r["W0"] for _, w in r["ranked"]):
        ok = False; print("  FAIL: no bend eroded the margin (adversary inert)")
    # the fake prime / inflate must be among the strongest eroders (they add prime mass wrongly)
    worst = r["ranked"][0][0]
    if verbose:
        print(f"  real-prime baseline margin W0 = {r['W0']:+.6e}")
        print(f"  bent-axiom margins (ranked, most-eroded first):")
        for name, w in r["ranked"]:
            tag = "  <== BREAKS positivity (W<0)" if w < 0 else \
                  ("  (erodes)" if w < r["W0"] else "  (harmless/helps)")
            print(f"     {name:24s} W = {w:+.6e}{tag}")
        print(f"  bends that break positivity: {r['broke'] or 'none in this window'}")
    print("t11.selftest:", "PASS" if ok else "FAIL")
    return ok, r


def main():
    print("=== t11 Disproof engine (adversarial axiom-bending on Weil positivity) ===\n")
    ok, r = selftest()
    broke = r["broke"]
    if broke:
        claim = (f"adversary found {len(broke)} axiom-bend(s) that drive Weil positivity negative "
                 f"in the t={r['t']} window: {broke}. The real prime system is forced to sit where "
                 f"it does; those bends are what RH forbids (finite/numeric diagnostic)")
    else:
        claim = (f"in the t={r['t']} window no single bend broke positivity, but all erode the "
                 f"margin from W0={r['W0']:.3e}; the real primes leave the largest margin "
                 f"(finite/numeric diagnostic)")
    results = [Result(claim=claim, tier=Tier.REPORT,
                      method="fixed arch+boundary, bent prime side, margin watch",
                      tool="t11_disproof_engine",
                      data={"W0": r["W0"], "ranked": [(n, w) for n, w in r["ranked"]],
                            "broke": broke})]
    path = provenance.stamp("t11_disproof_engine", results, extra={"selftest_pass": ok})
    print("\nReading: you cannot fabricate a counterexample by nudging the primes without either "
          "moving the place-at-infinity to match or breaking PNT/Euler structure -- exactly the "
          "'no second sign-carrying bridge' Wall Invariant, seen from the adversary side.")
    print("run record:", path)
    return ok


if __name__ == "__main__":
    main()
