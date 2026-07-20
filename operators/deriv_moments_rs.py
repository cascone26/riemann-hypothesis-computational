#!/usr/bin/env python3
"""
deriv_moments_rs.py -- HKO derivative-moment test at LARGE T via fast Riemann-Siegel.

Key fact: at a zero rho=1/2+i*gamma, zeta(rho)=0, so
    |zeta'(1/2 + i gamma)| = |Z'(gamma)|
where Z(t) = e^{i theta(t)} zeta(1/2+it) is the (real) Riemann-Siegel Z function.
So we evaluate Z'(gamma_n) at the KNOWN zeros (no root-finding) via vectorized
Riemann-Siegel (main sum + C0 remainder), central-difference for Z'.

This replaces mpmath (54 ms/zero -> ~30 h for 2M) with a numpy evaluator (minutes),
extending the log-lever from T~10^4 (L~7.4) to T~10^6 (L~12) so the exponent k(k+2)
(k=1->3, k=2->8, k=3->15) can actually be resolved.

Validated against mpmath zeta' on the first 300 zeros before trusting the sweep.
"""
import numpy as np
import json, time, sys

TWO_PI = 2 * np.pi

def theta(t):
    return 0.5*t*np.log(t/TWO_PI) - 0.5*t - np.pi/8 + 1.0/(48*t) + 7.0/(5760*t**3)

def rs_Z_chunk(t):
    """Riemann-Siegel Z(t) for array t (single chunk), main sum + C0 remainder."""
    a = np.sqrt(t/TWO_PI)
    N = np.floor(a).astype(np.int64)
    p = a - N
    th = theta(t)
    Nmax = int(N.max())
    ns = np.arange(1, Nmax+1)
    sqn = np.sqrt(ns)
    logn = np.log(ns)
    ang = th[:, None] - t[:, None]*logn[None, :]
    terms = np.cos(ang) / sqn[None, :]
    mask = ns[None, :] <= N[:, None]
    main = 2.0 * np.sum(terms*mask, axis=1)
    C0 = np.cos(TWO_PI*(p*p - p - 1.0/16)) / np.cos(TWO_PI*p)
    rem = np.where((N-1) % 2 == 0, 1.0, -1.0) * a**(-0.5) * C0
    return main + rem

def Zprime(gam, h=1e-3, chunk=4000):
    """|Z'(gamma)| via central difference, chunked to bound memory."""
    out = np.empty(len(gam))
    for i in range(0, len(gam), chunk):
        g = gam[i:i+chunk]
        zp = (rs_Z_chunk(g+h) - rs_Z_chunk(g-h)) / (2*h)
        out[i:i+chunk] = zp
    return out

# ---------- validation vs mpmath on first 300 zeros ----------
def validate():
    import mpmath as mp
    mp.mp.dps = 30
    print("[val] loading first 300 zeros (31 dp) + mpmath zeta' ...", flush=True)
    gam = []
    with open("/Users/scones/projects/riemann/data/zeros_100k.txt") as f:
        for i, line in enumerate(f):
            if i >= 300: break
            gam.append(float(line.strip()))
    gam = np.array(gam)
    zp_rs = np.abs(Zprime(gam))
    zp_mp = np.array([float(mp.fabs(mp.zeta(mp.mpf('0.5')+1j*mp.mpf(str(g)), derivative=1))) for g in gam])
    rel = np.abs(zp_rs - zp_mp) / np.maximum(zp_mp, 1e-12)
    print(f"[val] |Z'| vs mpmath |zeta'|: median rel err={np.median(rel):.2e}  "
          f"max={np.max(rel):.2e}  mean={np.mean(rel):.2e}", flush=True)
    return np.median(rel) < 5e-3

if __name__ == "__main__":
    NZ = int(sys.argv[1]) if len(sys.argv) > 1 else 2000000
    if not validate():
        print("[val] FAILED — Riemann-Siegel |Z'| disagrees with mpmath >0.5%. Aborting.", flush=True)
        sys.exit(1)
    print(f"[val] PASS — Riemann-Siegel evaluator trusted.\n", flush=True)

    print(f"[rs] loading {NZ} zeros (9 dp Odlyzko 2M) ...", flush=True)
    g = np.loadtxt("/Users/scones/projects/riemann/data/zeros_odlyzko_2M", max_rows=NZ)
    print(f"[rs] {len(g)} zeros, range {g[0]:.3f}..{g[-1]:.1f}", flush=True)

    t0 = time.time()
    absZp = np.abs(Zprime(g))
    print(f"[rs] |Z'| for {len(g)} zeros in {time.time()-t0:.0f}s", flush=True)
    absZp2 = absZp**2   # |zeta'|^2

    KS = [1, 2, 3]
    # ladder of T (cumulative over first m zeros), wide range for clean exponent
    ladder = [2000, 5000, 10000, 25000, 50000, 100000, 200000, 400000,
              700000, 1000000, 1500000, len(g)]
    ladder = sorted(set(m for m in ladder if m <= len(g)))
    out = {"n_zeros": len(g), "prediction_power": {k: k*(k+2) for k in KS}, "ladder_m": ladder, "ks": {}}
    print("\n" + "="*70)
    print("HKO DERIVATIVE MOMENTS at LARGE T (power k(k+2): 3, 8, 15)")
    print("="*70)
    for k in KS:
        w = absZp2**k
        Ls, logMs = [], []
        for m in ladder:
            Mk = float(np.mean(w[:m]))
            L = np.log(g[m-1]/TWO_PI)
            Ls.append(np.log(L)); logMs.append(np.log(Mk))
        Ls = np.array(Ls); logMs = np.array(logMs)
        A = np.vstack([Ls, np.ones_like(Ls)]).T
        slope_full = np.linalg.lstsq(A, logMs, rcond=None)[0][0]
        half = len(Ls)//2
        slope_top = np.linalg.lstsq(A[half:], logMs[half:], rcond=None)[0][0]
        pred = k*(k+2)
        print(f"\n k={k}: predicted exponent {pred}")
        print(f"   full-range slope = {slope_full:.3f}   top-half slope = {slope_top:.3f}   "
              f"(Lmax={np.exp(Ls[-1]):.2f})")
        out["ks"][k] = {"predicted": pred, "slope_full": float(slope_full),
                        "slope_top": float(slope_top),
                        "logL": Ls.tolist(), "logM": logMs.tolist()}
    p = "/Users/scones/projects/riemann/results/deriv_moments_rs.json"
    json.dump(out, open(p, "w"), indent=2)
    print(f"\nsaved -> {p}")
