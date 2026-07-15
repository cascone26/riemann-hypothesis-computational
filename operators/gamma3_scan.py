#!/usr/bin/env python3
"""
gamma3_scan.py — Test γ₃ resonance at lsq≈625.5 (γ₃=25.011, γ₃²=625.5).

Geometric lattice predictions for γ₃:
  n=-2: lsq = 625.5 × e^{-2π/25.011} = 625.5 × e^{-0.2513} = 625.5 × 0.7778 = 486.4
  n=-1: lsq = 625.5 × e^{-π/25.011}  = 625.5 × e^{-0.1257} = 625.5 × 0.8819 = 551.5
  n=0:  lsq = 625.5 (main resonance)
  n=+1: lsq = 625.5 × e^{+π/25.011}  = 625.5 × 1.1339   = 709.4

γ₁ lattice at these lsq:
  At lsq=486: γ₁ n ≈ log(486/200)/(π/14.135) = 3.94 → not a lattice point
  At lsq=552: γ₁ n ≈ log(552/200)/(π/14.135) = 4.47 → between n=4 and n=5
  At lsq=625: γ₁ n ≈ log(625/200)/(π/14.135) = 5.01 → EXACTLY n=5 lattice point!

γ₂ lattice at these lsq:
  At lsq=625: γ₂ n ≈ log(625/442)/(π/21.022) = 2.43 → not a lattice point
  At lsq=552: γ₂ n ≈ log(552/442)/(π/21.022) = 1.51 → not a lattice point

On Olympus: use N=min(5*lsq, 3200). On Mac: use N=min(3*lsq, 1000).
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
N_CAP = 3200 if ON_OLYMPUS else 1000

OUTFILE = os.path.expanduser("~/projects/riemann/results/gamma3_scan.out")

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
    gamma3 = 25.01085758
    gamma3_sq = gamma3**2  # 625.54

    lsq_list = [
        400, 440, 460,             # below γ₂ main resonance (control)
        486,                        # γ₃ n=-2 prediction
        520, 540,                   # interpolation
        552,                        # γ₃ n=-1 prediction
        575, 600,                   # interpolation
        625, 626,                   # γ₃ n=0 main resonance
        660, 700,                   # γ₃ n>0 predictions
        710,                        # γ₃ n=+1 prediction
        750, 800,                   # recovery
    ]

    def N_for(lsq):
        return min(5 * lsq if ON_OLYMPUS else 3 * lsq, N_CAP)

    results = []
    t0 = time.time()
    mode = "OLYMPUS" if ON_OLYMPUS else "MAC"
    print(f"γ₃ resonance scan [{mode}]: lsq=400-800 for deep minimum at lsq≈625.5")
    print(f"γ₁=14.135 γ₂=21.022 γ₃=25.011 → lsq_resonance={gamma3_sq:.1f}")
    print(f"{'lsq':>5} {'N':>5} {'sc':>12} {'pred':>12} {'ratio':>7} {'lam0':>12} {'dim':>6} {'dt':>5}")
    hdr = "-" * 80
    print(hdr)
    with open(OUTFILE, 'w') as f:
        f.write(f"gamma3_scan [{mode}] {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(hdr + "\n")

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
        line = f"{lsq:5d} {N:5d}: sc={sc:.4e}  pred={sc_pred:.4e}  ratio={ratio:.3f}  lam0={lam0:.2e}  dim={dim}  ({dt:.0f}s)"
        print(line, flush=True)
        with open(OUTFILE, 'a') as f:
            f.write(line + "\n")
        results.append({'lsq': lsq, 'N': N, 'sc': sc, 'ratio': round(ratio,5) if not np.isnan(ratio) else None,
                        'lam0': lam0, 'dim': dim, 'dt': round(dt,1)})

    total = time.time() - t0
    print(f"\nTotal: {total:.1f}s")
    jfile = OUTFILE.replace('.out', '.json')
    with open(jfile, 'w') as f:
        json.dump({'lsq_list': lsq_list, 'results': results, 'total': round(total,1),
                   'gamma3': gamma3, 'gamma3_sq': gamma3_sq, 'mode': mode}, f, indent=2)
    print(f"JSON: {jfile}")

if __name__ == "__main__":
    main()
