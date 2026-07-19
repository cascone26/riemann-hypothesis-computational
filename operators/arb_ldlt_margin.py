#!/usr/bin/env python3
"""CERTIFIED POSITIVITY-MARGIN SCALING of the real Weil operator tau.

Additive extension of external/groskin_2607.02828/code/arb_ldlt_certify.py (imports its exact,
rigorous Arb assembly build_arb_tau — does NOT modify it). The canonical certifier only reports the
PD boolean (n_neg==0). Here we ALSO capture the SMALLEST PIVOT of the interval LDL^T: its rigorous
lower bound is a certified positivity margin. Tracking how that margin scales with N (test-basis
half-size) and c (prime cutoff) is the honest crux question of this framework:

  margin plateaus / decays slowly  -> consistent with a UNIFORM positive lower bound (RH-favorable)
  margin decays fast toward 0       -> positivity may fail in the N,c -> infinity limit (danger)
  a pivot straddles 0 at high prec  -> genuine near-singularity, escalate

Rigorous: every pivot is an Arb ball; we report min over i of the ball's LOWER bound. If any pivot
ball straddles 0 the point (N,c,prec) is UNDETERMINED (needs more precision), reported honestly.
"""
import sys, os, time, json, argparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "external", "groskin_2607.02828", "code"))
from flint import arb, arb_mat, ctx
from arb_ldlt_certify import build_arb_tau


def inertia_with_margin(A, DIM):
    """Interval LDL^T that also returns (n_pos, n_neg, undet_idx, min_pivot_lowerbound, min_pivot_mid)."""
    d = [None] * DIM
    Lf = [[arb(0)] * DIM for _ in range(DIM)]
    n_pos = n_neg = 0
    min_low = None; min_mid = None
    for i in range(DIM):
        s = A[i, i]
        for k in range(i):
            s = s - Lf[i][k] * Lf[i][k] * d[k]
        d[i] = s
        if s > 0:
            n_pos += 1
            low = float(s.mid()) - float(s.rad())   # rigorous lower bound on this pivot's true value
            mid = float(s.mid())
            if min_low is None or low < min_low:
                min_low = low; min_mid = mid
        elif s < 0:
            n_neg += 1
        else:
            return n_pos, n_neg, i, None, None
        for j in range(i + 1, DIM):
            t = A[j, i]
            for k in range(i):
                t = t - Lf[j][k] * Lf[i][k] * d[k]
            Lf[j][i] = t / d[i]
    return n_pos, n_neg, None, min_low, min_mid


def certify_margin(c, N, prec):
    ctx.prec = prec
    A, DIM = build_arb_tau(c, N, prec)
    npos, nneg, undet, mlow, mmid = inertia_with_margin(A, DIM)
    return dict(c=c, N=N, dim=DIM, prec=prec, n_pos=npos, n_neg=nneg,
                undetermined_pivot=undet,
                certified_pd=bool(undet is None and nneg == 0),
                min_pivot_lower=mlow, min_pivot_mid=mmid)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--prec", type=int, default=1000)
    p.add_argument("--out", type=str, default="results/margin_scaling.json")
    args = p.parse_args()

    results = {"grow_N": [], "grow_c": []}
    print("=== CERTIFIED positivity-margin scaling (Arb interval LDL^T) ===", flush=True)

    print("\n[A] fixed c=14, growing N:", flush=True)
    for N in [8, 16, 24, 32, 48, 64, 96, 128]:
        t0 = time.time()
        r = certify_margin(14, N, args.prec)
        results["grow_N"].append(r)
        tag = "PD" if r["certified_pd"] else ("UNDET@%s" % r["undetermined_pivot"] if r["undetermined_pivot"] is not None else "INDEF n_neg=%d" % r["n_neg"])
        print("  N=%4d dim=%4d  min_pivot_lower=%s  [%s]  (%.0fs)"
              % (N, r["dim"], ("%.3e" % r["min_pivot_lower"]) if r["min_pivot_lower"] is not None else "n/a", tag, time.time() - t0), flush=True)

    print("\n[B] fixed N=64, growing c:", flush=True)
    for c in [8, 14, 25, 50, 100]:
        t0 = time.time()
        r = certify_margin(c, 64, args.prec)
        results["grow_c"].append(r)
        tag = "PD" if r["certified_pd"] else ("UNDET@%s" % r["undetermined_pivot"] if r["undetermined_pivot"] is not None else "INDEF n_neg=%d" % r["n_neg"])
        print("  c=%4d dim=%4d  min_pivot_lower=%s  [%s]  (%.0fs)"
              % (c, r["dim"], ("%.3e" % r["min_pivot_lower"]) if r["min_pivot_lower"] is not None else "n/a", tag, time.time() - t0), flush=True)

    os.makedirs("results", exist_ok=True)
    json.dump(results, open(args.out, "w"), indent=2)
    print("\nwrote %s" % args.out, flush=True)

    # honest trend read
    import math
    def read(rows, key, label):
        good = [(r[key], r["min_pivot_lower"]) for r in rows if r["certified_pd"] and r["min_pivot_lower"] and r["min_pivot_lower"] > 0]
        if len(good) < 3:
            print("  [%s] insufficient certified points for a trend (%d)" % (label, len(good)), flush=True); return
        xs = [math.log(x) for x, _ in good]; ys = [math.log(y) for _, y in good]
        n = len(xs); sx = sum(xs); sy = sum(ys); sxx = sum(x*x for x in xs); sxy = sum(x*y for x, y in zip(xs, ys))
        slope = (n*sxy - sx*sy) / (n*sxx - sx*sx)
        print("  [%s] certified margin log-log slope = %+.3f (%s)" % (label, slope,
              "near-flat: uniform-bound-favorable" if slope > -0.3 else
              "mild decay" if slope > -1.0 else "steep decay: probe the limit"), flush=True)
    print("\n=== TREND ===", flush=True)
    read(results["grow_N"], "N", "grow N")
    read(results["grow_c"], "c", "grow c")


if __name__ == "__main__":
    main()
