#!/usr/bin/env python3
"""§178c: pin the TRUE dist(k_λ, θ_x). §178b saw overlap=1.000000, L2 dist≤7e-4 decreasing with c. Is that floor the
REAL correction θ_x−k_λ (the object Connes §6.6 must bound), or an artifact of (a) CvS eigenvector ψ under-resolution at
finite N, or (b) trapezoid quadrature at finite M? float64 prolates are ~1e-15 accurate — 11 orders below 1e-4 — so the
floor is NOT prolate precision. Sweep N (CvS resolution) × M (quadrature) at fixed c: if dist stabilizes ⇒ that stable
value = true dist(k_λ,θ_x); if it keeps dropping with N ⇒ ψ was under-resolved and the true residual is smaller."""
import sys, time
import numpy as np
sys.path.insert(0, ".venv_cvs/lib/python3.14/site-packages")
from bridge_overlap import cvs_ground, k_lambda_on, f_psi, overlap

if __name__ == "__main__":
    c   = int(sys.argv[1]) if len(sys.argv) > 1 else 13
    dps = int(sys.argv[2]) if len(sys.argv) > 2 else 90
    T   = int(sys.argv[3]) if len(sys.argv) > 3 else 400
    Ns  = [int(x) for x in sys.argv[4].split(",")] if len(sys.argv) > 4 else [40, 60, 80, 100]
    Ms  = [int(x) for x in sys.argv[5].split(",")] if len(sys.argv) > 5 else [6000, 12000, 24000]
    lam = np.sqrt(c)
    print(f"§178c dist-floor sweep  c={c} λ={lam:.3f} dps={dps} T={T}\n")
    print(f"{'N':>5} {'log10λmin':>10} " + " ".join(f"{'dist@M='+str(M):>14}" for M in Ms) + f"  {'1-overlap@Mmax':>15}")
    for N in Ns:
        t0 = time.time()
        lmin_log, ev, L = cvs_ground(c, N, T, dps)
        row = []
        ov_last = None
        for M in Ms:
            x = np.linspace(0.0, L, M + 1)
            fp = f_psi(x, ev, L)
            kA = k_lambda_on(x - L / 2.0, lam)
            ov, d = overlap(fp, kA, x)
            row.append(d); ov_last = ov
        dstr = " ".join(f"{d:>14.3e}" for d in row)
        print(f"{N:>5} {lmin_log:>10.2f} {dstr}  {1-ov_last:>15.3e}  ({time.time()-t0:.0f}s)", flush=True)
    print("\nREAD: across a ROW (fixed N, rising M) → quadrature convergence; if dist flat in M, quadrature is NOT the floor.")
    print("Down a COLUMN (fixed M, rising N) → if dist keeps dropping, ψ was under-resolved (true residual smaller); if it")
    print("stabilizes, that value = the TRUE dist(k_λ,θ_x) = the §6.6 correction magnitude at this c. Compare across c for the law.")
