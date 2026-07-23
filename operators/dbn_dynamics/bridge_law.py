#!/usr/bin/env python3
"""§178d: the dist(k_λ,θ_x) DECAY LAW vs c — the Connes §6.6 correction rate.

§178b/c established: Connes' analytical k_λ = E(h_λ) IS the CvS near-null eigenvector ψ (overlap 1−2e-14 at c=13).
The residual L2 dist(k_λ,ψ) DROPS with N (ψ under-resolution, §178c: 3e-3→9e-5→2e-7 at N=40,60,80 for c=13) and the
M(quadrature) axis is FLAT (ruled out as the floor). So the only knob left is N (CvS resolution), gated by dps
(precision): too-low dps STARVES the near-null and 1-overlap gets WORSE (§178c: N=100 dps90 regressed to 2.5e-12).

This driver: for each c, sweep N upward with dps auto-scaled to the near-null depth, at fixed M. Per c it reports the
CLEANEST row (minimum 1−overlap = best-resolved ψ) as the tightest dist(k_λ,θ_x) bound, and flags whether that value
has PLATEAUED (dist stable across the last two clean N ⇒ the TRUE correction θ_x−k_λ is visible) or is still
N-LIMITED (dist still dropping ⇒ the value is an UPPER BOUND, true correction smaller). Then fits log10(dist) vs c
over the plateaued points → the decay rate = the §6.6 k_λ→θ_x correction rate (the analytic target).

Run from ~/projects/riemann:  cvs_venv/bin/python3 operators/dbn_dynamics/bridge_law.py <c-list> <N-list> [M] [Tint]
  e.g.  ... bridge_law.py 5,8,13 40,60,80,100 12000 400
dps is auto-scaled per (c,N); pass a matching dps-list as arg 5 to override."""
import sys, time, math
import numpy as np
sys.path.insert(0, ".venv_cvs/lib/python3.14/site-packages")
from bridge_overlap import cvs_ground, k_lambda_on, f_psi, overlap

LN10 = math.log(10.0)

def depth_est(c):
    """Fuchs envelope magnitude (Connes §6.4 / §171): |6.37 − 4πc/ln10 + (9/2)ln c/ln10| = expected |log10 λmin|."""
    return abs(6.37 - 4 * math.pi * c / LN10 + 4.5 * math.log(c) / LN10)

def dps_for(c, N):
    """Auto-scale precision to resolve the near-null eigenVECTOR: depth + a resolution margin growing with N.
    Anchored to §178c empirics (c=13: N=80 clean at dps90; N=100 needed ~130) — deliberately conservative (safe-slow)."""
    return max(60, int(depth_est(c) + 0.7 * N + 15))

def run_c(c, Ns, M, T, dps_override=None):
    lam = math.sqrt(c)
    print(f"\n=== c={c}  λ={lam:.4f}  depth_est(|log10 λmin|)≈{depth_est(c):.1f}  M={M} T={T} ===")
    print(f"{'N':>5} {'dps':>5} {'log10λmin':>10} {'dist':>13} {'1-overlap':>13}  {'note':>16}  {'t':>6}")
    rows = []
    prev_1mov = None
    for i, N in enumerate(Ns):
        dps = dps_override[i] if dps_override else dps_for(c, N)
        t0 = time.time()
        lmin_log, ev, L = cvs_ground(c, N, T, dps)
        x = np.linspace(0.0, L, M + 1)
        fp = f_psi(x, ev, L)
        kA = k_lambda_on(x - L / 2.0, lam)
        ov, d = overlap(fp, kA, x)
        one_m_ov = 1.0 - ov
        # dps-starvation flag: 1-overlap got WORSE than the previous (lower) N ⇒ near-null under-resolved at this dps
        starved = prev_1mov is not None and one_m_ov > prev_1mov
        note = "STARVED(raise dps)" if starved else "clean"
        print(f"{N:>5} {dps:>5} {lmin_log:>10.2f} {d:>13.3e} {one_m_ov:>13.3e}  {note:>16}  {time.time()-t0:>5.0f}s",
              flush=True)
        rows.append(dict(N=N, dps=dps, lmin_log=lmin_log, dist=d, one_m_ov=one_m_ov, starved=starved))
        prev_1mov = one_m_ov
    # cleanest row = min 1-overlap among non-starved
    clean = [r for r in rows if not r["starved"]] or rows
    best = min(clean, key=lambda r: r["one_m_ov"])
    # plateau test: is dist stable across the two deepest clean N (within 30%)?
    clean_sorted = sorted(clean, key=lambda r: r["N"])
    plateau = False
    if len(clean_sorted) >= 2:
        a, b = clean_sorted[-2]["dist"], clean_sorted[-1]["dist"]
        ratio = a / b if b > 0 else float("inf")
        plateau = 0.7 <= ratio <= 1.43
    kind = "PLATEAU (true correction visible)" if plateau else "N-LIMITED (dist=upper bound, true corr smaller)"
    print(f"  → c={c}: best-resolved dist={best['dist']:.3e} at N={best['N']} (1-ov={best['one_m_ov']:.2e}) — {kind}")
    return dict(c=c, dist=best["dist"], N=best["N"], one_m_ov=best["one_m_ov"], plateau=plateau)

if __name__ == "__main__":
    clist = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [5, 8, 13]
    Ns    = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [40, 60, 80, 100]
    M     = int(sys.argv[3]) if len(sys.argv) > 3 else 12000
    T     = int(sys.argv[4]) if len(sys.argv) > 4 else 400
    dpsov = [int(x) for x in sys.argv[5].split(",")] if len(sys.argv) > 5 else None
    print(f"§178d dist(k_λ,θ_x) DECAY LAW  c={clist} N={Ns} M={M} T={T}  dps={'auto' if not dpsov else dpsov}")
    summary = [run_c(c, Ns, M, T, dpsov) for c in clist]
    print("\n================ §178d SUMMARY ================")
    print(f"{'c':>4} {'best dist':>12} {'at N':>5} {'1-overlap':>12}  {'status':>10}")
    for s in summary:
        print(f"{s['c']:>4} {s['dist']:>12.3e} {s['N']:>5} {s['one_m_ov']:>12.2e}  {'plateau' if s['plateau'] else 'N-limit'}")
    plat = [s for s in summary if s["plateau"]]
    fitset = plat if len(plat) >= 2 else summary  # prefer plateaued points; else fit the upper-bound trend (note it)
    if len(fitset) >= 2:
        cs = np.array([s["c"] for s in fitset], float)
        ld = np.array([math.log10(s["dist"]) for s in fitset])
        A = np.vstack([cs, np.ones_like(cs)]).T
        (slope, intc), *_ = np.linalg.lstsq(A, ld, rcond=None)
        tag = "PLATEAUED points (true §6.6 correction rate)" if len(plat) >= 2 else \
              "UPPER-BOUND points (N-limited — rate is a LOWER bound on the true decay)"
        print(f"\nFIT log10 dist(c) = {slope:.4f}·c + {intc:.3f}   ⇒  dist ~ 10^({slope:.3f}c) = e^({slope*LN10:.3f}c)")
        print(f"     over {tag}")
        print(f"     decay rate α = {-slope*LN10:.4f} per unit c  (dist ~ e^(−α c));  compare to Fuchs 4π/ln10·ln10=4π≈{4*math.pi:.3f}? — for the analytic read")
    else:
        print("\nNeed ≥2 comparable points to fit; extend the c-list.")
    print("\nNEXT: if PLATEAUED at small c but N-LIMITED at large c, push N (with scaled dps) at the large-c points to")
    print("convert upper bounds → plateaus; the clean law across c is the §6.6 k_λ→θ_x correction rate (analytic target).")
