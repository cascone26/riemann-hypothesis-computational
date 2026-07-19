#!/usr/bin/env python3
"""CRUX DIAGNOSTIC: min-eigenvalue scaling of the REAL Weil-positivity operator tau=A+B-P.

RH <=> tau stays PSD as BOTH N (test-basis half-size, matrix is 2N+1) and c (prime cutoff, L=log c)
grow without bound. tau is known PD at finite (N,c). The honest open question in THIS framework is
whether lambda_min(tau) is bounded below by a positive constant (RH-favorable / uniform positivity)
or is DECAYING toward 0 as we refine (the danger sign — positivity could fail in the limit).

This does not prove or disprove anything; it measures the actual trend on the actual object, which is
exactly the kind of real signal (not a re-verification) worth having. Prints a grid + a fitted trend.
"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
import numpy as np
from critical_eigvec_id import build_tau_np

def mineig(N, c):
    tau, ns = build_tau_np(N, c)
    tau = (tau + tau.T) / 2.0            # symmetrize against float noise
    w = np.linalg.eigvalsh(tau)
    return float(w[0]), float(w[-1]), tau.shape[0]

def main():
    print("=== tau min-eigenvalue scaling (real operator) ===", flush=True)
    # 1) refine N at fixed c: does more test-space kill positivity?
    print("\n[A] fixed c=14, growing N (test-basis refinement):", flush=True)
    c = 14
    rowsN = []
    for N in [16, 24, 32, 48, 64, 96, 128, 192, 256]:
        t0 = time.time()
        lo, hi, D = mineig(N, c)
        rowsN.append((N, lo))
        print(f"  N={N:4d} dim={D:4d}  lambda_min={lo:+.6e}  lambda_max={hi:+.4e}  cond~{hi/abs(lo):.2e}  ({time.time()-t0:.1f}s)", flush=True)

    # 2) refine c at fixed N: does adding more primes kill positivity?
    print("\n[B] fixed N=96, growing c (prime-cutoff refinement):", flush=True)
    N = 96
    rowsC = []
    for c in [8, 14, 25, 50, 100, 200, 400]:
        t0 = time.time()
        lo, hi, D = mineig(N, c)
        rowsC.append((c, lo))
        print(f"  c={c:4d}  lambda_min={lo:+.6e}  lambda_max={hi:+.4e}  ({time.time()-t0:.1f}s)", flush=True)

    # trend read: is lambda_min heading to 0 (danger) or leveling off (favorable)?
    def trend(rows, label):
        xs = np.array([r[0] for r in rows], float)
        ys = np.array([r[1] for r in rows], float)
        allpos = bool(np.all(ys > 0))
        # log-log slope of |lambda_min| vs the refinement parameter (only meaningful if monotone-ish)
        with np.errstate(all="ignore"):
            slope = np.polyfit(np.log(xs), np.log(np.abs(ys)), 1)[0]
        last3 = ys[-3:]
        leveling = bool(np.std(last3) < 0.25 * abs(np.mean(last3)))  # last 3 within 25% => plateau
        print(f"\n  [{label}] all lambda_min>0: {allpos} | log-log slope: {slope:+.3f} | "
              f"last-3 plateau: {leveling} | last value: {ys[-1]:+.3e}", flush=True)
        return allpos, slope, leveling, float(ys[-1])

    aN = trend(rowsN, "grow N")
    aC = trend(rowsC, "grow c")

    print("\n=== READ ===", flush=True)
    if aN[0] and aC[0]:
        if aN[2] or aC[2] or (aN[1] > -0.15 and aC[1] > -0.15):
            print("FAVORABLE: lambda_min stays positive and shows plateau / slow decay — consistent with a "
                  "UNIFORM positive lower bound (the RH-favorable scenario). Worth a rigorous bound attempt.", flush=True)
        else:
            print("WATCH: positive so far but lambda_min is DECAYING with refinement (steep log-log slope) — "
                  "the limit is not obviously bounded away from 0. This is the honest danger direction to probe.", flush=True)
    else:
        print("ANOMALY: a NEGATIVE lambda_min appeared at finite (N,c) — that would be a genuine surprise "
              "(numerical or real). Flag for Arb interval re-verification before believing it.", flush=True)

if __name__ == "__main__":
    main()
