#!/usr/bin/env python3
"""Arb certification for the log3-rung findings (§476/§476b -> Handle).
CERT A (load-bearing indefiniteness): at 4W=1.32 (past log3), the p2-only Weil form (n=3^k terms removed)
  is INDEFINITE: exhibit a float direction v and rigorously enclose Q_no3(v) < 0. One quadratic-form value,
  robust to entry radii. Promotes §476b(1) to a rigorous statement.
CERT B (gate 2, p=2 rung PD): at 4W=0.95 in [log2, log3), interval-LDL^T certify the full Weil form PD on
  the probe span (CC's proven rung, certified in our independent instrument).
Reuses operators/arb_fejer_cert2.py machinery with W/TW/c_tri/RC re-set per run.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import arb_fejer_cert2 as ac
import l2_fejer as fj

KILL3 = set()

def set_W(wval, rc):
    ac.W = arb(wval); ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W; ac.RC = arb(rc)

# patch prime_entry to honor KILL3
def prime_entry(cj, ck):
    s = arb(0); fourw = 4*ac.W
    for n in range(2, ac.NM+1):
        if ac.LAM[n] is None or n in KILL3: continue
        u = arb(n).log()
        if u > fourw: continue
        s = s + ac.LAM[n]/arb(n).sqrt()*ac.gpair_closed_arb(cj, ck, u)
    return 2*s

def entry(cj, ck):
    return ac.pole_entry(cj, ck) + ac.arch_entry(cj, ck) - prime_entry(cj, ck)

def certA():
    global KILL3
    print("=== CERT A: p2-only form INDEFINITE at 4W=1.32 (load-bearing p=3) ===", flush=True)
    fourW = 1.32; s = 4.0
    # float pass to get the witness direction on the SAME probe set
    fj.W = fourW/4.0
    cs_f = fj.probes((14.0, 44.0), s)
    scale = np.ones_like(fj.LAMW); scale[np.isin(fj.NS, [3, 9, 27, 81])] = 0.0
    f = fj.build(cs_f, lam_scale=scale)
    FU = 0.5*(f['full']+f['full'].T)
    ev, V = np.linalg.eigh(FU)
    print(f"float: n={len(cs_f)} probes (s={s}), min eig of no-3 form = {ev[0]:.4e}", flush=True)
    if ev[0] >= 0:
        print("float min-eig >= 0 at this spacing -- indefiniteness not visible on this span; widen span/spacing.")
        return
    v = V[:, 0]
    # rigorous Q(v) with arb machinery
    set_W(fourW/4.0, 300)
    KILL3 = {3, 9, 27, 81}
    cs = ac.probes(arb(s))
    n = len(cs); assert n == len(cs_f)
    va = [arb(float(x)) for x in v]
    Q = arb(0)
    for j in range(n):
        for k in range(j, n):
            e = entry(cs[j], cs[k])
            Q = Q + (arb(1) if j == k else arb(2))*va[j]*va[k]*e
        print(f"  row {j} done", flush=True)
    KILL3 = set()
    print(f"RIGOROUS enclosure Q_no3(v) = {Q}", flush=True)
    if Q < 0:
        print("CERT A PASSED: Q_no3(v) < 0 rigorously => p2-only Weil form INDEFINITE past log3 (HANDLE).")
    else:
        print("CERT A INCONCLUSIVE: enclosure does not exclude 0 (raise prec/RC).")

def certB():
    print("\n=== CERT B: full Weil form PD at 4W=0.95 (p=2 rung, CC-proven, our independent cert) ===", flush=True)
    set_W(0.95/4.0, 400)
    ok = ac.certify("4.0")
    print("CERT B", "PASSED (HANDLE)" if ok else "INCONCLUSIVE (radii vs margin — raise RC/prec)")

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "AB"
    if "A" in which: certA()
    if "B" in which: certB()
