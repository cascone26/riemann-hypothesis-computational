#!/usr/bin/env python3
"""
extended_powerlaw.py — Extend sc power-law fit to lsq=100-350.

Current fit: sc ≈ 3.073e-3 / L^1.9151 (14 points, lsq=5-100)
Goal: Extend to lsq=100-350 to:
  (a) Validate L^{-2} convergence at larger scale
  (b) Identify if resonance minima (lsq=200,250,312) break the power law
  (c) Get N-saturation behavior at each scale

Uses N=min(5*lsq, 2000) on Olympus or N=min(2*lsq, 600) on Mac.
"""
import os, json, time
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
THRESH = 0.02
ON_OLYMPUS = (os.uname().nodename == 'olympus' or os.path.exists('/home/cobo'))
G = 60000 if ON_OLYMPUS else 40000
N_MULT = 5 if ON_OLYMPUS else 2
N_CAP = 2000 if ON_OLYMPUS else 600

OUTFILE = os.path.expanduser("~/projects/riemann/results/extended_powerlaw.out")

def vm_table(kmax):
    tab = []
    for k in range(2, int(kmax) + 1):
        n, p, f, ok = k, 0, 2, True
        while f * f <= n:
            if n % f == 0:
                if p == 0: p = f
                elif p != f: ok = False; break
                while n % f == 0: n //= f
            f += 1
        if n > 1:
            if p == 0: p = n
            elif p != n: ok = False
        if ok and p > 0:
            tab.append((k, float(mp.log(p))))
    return tab

def measure(lsq, N):
    L  = float(mp.log(lsq))
    sh = float(np.sinh(L / 4) ** 2)
    VM = vm_table(lsq)
    idx = np.arange(-N, N + 1, dtype=float)
    ns  = idx
    R   = TWO_PI * N / L + 80.0
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]; w = np.full(G, h); w[0] = w[-1] = 0.5 * h
    Om = np.real(_digamma(0.25 + 0.5j * grid)) - LOGPI
    d = (TWO_PI * ns / L)[:, None] - grid[None, :]
    tiny = np.abs(d) < 1e-12
    E = np.where(tiny, np.sqrt(L),
                 (1 / np.sqrt(L)) * np.exp(0.5j * d * L)
                 * 2 * np.sin(0.5 * d * L) / (d + 1e-300 * tiny))
    arch = ((E * (w * Om)) @ E.conj().T).real / (2.0 * np.pi)
    B = (32.0 * L * sh * (L ** 2 - 16.0 * np.pi ** 2 * np.outer(ns, ns)) /
         np.outer(L ** 2 + 16.0 * np.pi ** 2 * ns ** 2,
                  L ** 2 + 16.0 * np.pi ** 2 * ns ** 2))
    if VM:
        lkpp = np.array([np.log(k) for k, _ in VM])
        wts  = np.array([lk * k ** (-0.5) for k, lk in VM])
        ph   = TWO_PI * np.outer(ns, lkpp) / L
        S    = np.sin(ph) @ wts
        dn   = ns[:, None] - ns[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            P = (S[None, :] - S[:, None]) / (np.pi * dn)
        P[dn == 0] = 0.0
        np.fill_diagonal(P, np.cos(ph) @ (2 * wts * (1 - lkpp / L)))
    else:
        P = np.zeros((len(idx), len(idx)))
    AB = arch + B
    tau = (AB - P); tau = (tau + tau.T) / 2.0
    ev, V = scipy.linalg.eigh(tau)
    lam0 = float(ev[0])
    null_mask = np.abs(ev) < THRESH
    dim = int(null_mask.sum())
    if dim > 0:
        Vn = V[:, null_mask]
        ABp = Vn.T @ AB @ Vn
        Pp  = Vn.T @ P  @ Vn
        sc = float(np.linalg.norm(ABp - Pp, 'fro') / (np.linalg.norm(ABp, 'fro') + 1e-100))
    else:
        sc = float('nan')
    return sc, lam0, dim

def main():
    # Extend power law fit: cover 100-350 with focus on non-resonance lsq
    # Avoid lsq=200,250,312 for the MAIN fit (they are resonance points)
    # Include them labeled separately
    lsq_list = [
        110, 120, 130, 140, 150,   # 100-150 range
        165, 175, 190,              # below γ₁ resonance (should be spikes if γ₁ n<0 lattice)
        200,                        # γ₁ main resonance (special)
        215, 225, 235,              # between γ₁ n=0 and n=1
        250,                        # γ₁ n=1 (special)
        270, 290, 310,              # between γ₁ n=1 and n=2
        312,                        # γ₁ n=2 (special)
        330, 350,                   # above γ₁ n=2
    ]

    def N_for(lsq):
        return min(N_MULT * lsq, N_CAP)

    results = []
    t0 = time.time()
    mode = "OLYMPUS" if ON_OLYMPUS else "MAC"
    print(f"Extended power-law scan [{mode}]: lsq=100-350, N=min({N_MULT}*lsq,{N_CAP})")
    print(f"{'lsq':>5} {'N':>5} {'sc':>12} {'pred':>12} {'ratio':>7} {'flag':>10} {'dt':>5}")
    print("-" * 75)
    with open(OUTFILE, 'w') as f:
        f.write(f"extended_powerlaw [{mode}] {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    resonance_pts = {200: 'γ₁n=0', 250: 'γ₁n=1', 312: 'γ₁n=2'}

    for lsq in lsq_list:
        N = N_for(lsq)
        L = float(mp.log(lsq))
        t1 = time.time()
        try:
            sc, lam0, dim = measure(lsq, N)
        except Exception as e:
            print(f"{lsq:5d}: ERROR {e}")
            continue
        dt = time.time() - t1
        sc_pred = 3.073e-3 / L ** 1.9151
        ratio = sc / sc_pred if sc_pred > 0 and not np.isnan(sc) else float('nan')
        flag = resonance_pts.get(lsq, '')
        # γ₁ n<0 lattice: predict spikes at 160(n=-1), 128(n=-2) — label nearby
        if 150 <= lsq <= 190:
            flag = 'γ₁precursor'
        line = f"{lsq:5d} {N:5d}: sc={sc:.4e}  pred={sc_pred:.4e}  ratio={ratio:7.3f}  {flag:12}  ({dt:.0f}s)"
        print(line, flush=True)
        with open(OUTFILE, 'a') as f:
            f.write(line + "\n")
        results.append({'lsq': lsq, 'N': N, 'sc': sc, 'ratio': round(ratio,5) if not np.isnan(ratio) else None,
                        'lam0': lam0, 'flag': flag, 'dt': round(dt,1)})

    # Fit power law EXCLUDING resonance points
    non_res = [(r['lsq'], r['sc']) for r in results if r['flag'] == '' and r['sc'] is not None and r['sc'] > 0]
    if len(non_res) >= 5:
        import numpy.linalg as nla
        lsqs = np.array([x[0] for x in non_res])
        scs  = np.array([x[1] for x in non_res])
        Ls   = np.log(lsqs)
        log_sc = np.log(scs)
        A    = np.column_stack([np.ones_like(Ls), Ls])
        sol  = nla.lstsq(A, log_sc, rcond=None)[0]
        alpha = -sol[1]
        C    = np.exp(sol[0])
        residuals = log_sc - (sol[0] + sol[1] * Ls)
        se   = np.std(residuals) / np.sqrt(len(Ls))
        print(f"\nPower-law fit ({len(non_res)} non-resonance points, lsq={lsqs[0]:.0f}-{lsqs[-1]:.0f}):")
        print(f"  sc ≈ {C:.4e} / L^{alpha:.4f} ± {se:.4f}")
        print(f"  vs prior: sc ≈ 3.073e-3 / L^1.9151")
        with open(OUTFILE, 'a') as f:
            f.write(f"\nPower-law fit: sc ≈ {C:.4e} / L^{alpha:.4f} ± {se:.4f}  ({len(non_res)} pts)\n")

    total = time.time() - t0
    print(f"\nTotal: {total:.1f}s  →  {OUTFILE}")
    jfile = OUTFILE.replace('.out', '.json')
    with open(jfile, 'w') as f:
        json.dump({'results': results, 'total': round(total,1), 'mode': mode}, f, indent=2)

if __name__ == "__main__":
    main()
