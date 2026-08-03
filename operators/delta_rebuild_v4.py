"""Delta rebuild v4: PURE exact-Lambda grid scan (no RS anywhere).
Grid 0.08 over [8,260] (min known gap 0.24), sign changes -> exact bisection, all
parallel with top-level workers (spawn-safe). Gates: S-band, LMFDB leading zeros, readout."""
import numpy as np
from scipy.special import loggamma as slg
from multiprocessing import Pool
TWO_PI = 2*np.pi
TMAX = 260.0

_poly = [1] + [0]*60
for _n in range(1, 61):
    for _ in range(24):
        _new = _poly[:]
        for _i in range(60, _n-1, -1):
            _new[_i] -= _poly[_i-_n]
        _poly = _new
TAU = [int(_poly[n-1]) for n in range(1, 61)]

def theta_v(t):
    return np.imag(slg(6.0 + 1j*np.asarray(t, dtype=float))) - np.asarray(t, dtype=float)*np.log(TWO_PI)

def Lam(t):
    import mpmath as mp
    mp.mp.dps = 60 + int(0.75*t)
    s = mp.mpf(6) + 1j*mp.mpf(t)
    tot = mp.mpf(0)
    for n, tn in enumerate(TAU, start=1):
        x = mp.mpf(TWO_PI*n)
        tot += tn*mp.re(mp.gammainc(s, x, mp.inf)*mp.power(x, -s)
                        + mp.gammainc(12-s, x, mp.inf)*mp.power(x, s-12))
    return tot

def lam_sign(t):
    return 1 if Lam(t) > 0 else -1

def refine(bracket):
    a, b, fa_sign = bracket
    for _ in range(22):
        m = 0.5*(a+b)
        fm = 1 if Lam(m) > 0 else -1
        if fm == fa_sign: a = m
        else: b = m
        if b - a < 5e-6: break
    return 0.5*(a+b)

def sgn_slope(gz):
    return 1.0 if Lam(gz+0.002) > Lam(gz-0.002) else -1.0

if __name__ == "__main__":
    grid = list(np.arange(8.0, TMAX, 0.08))
    with Pool(8) as p:
        signs = p.map(lam_sign, grid, chunksize=16)
    signs = np.array(signs)
    idx = np.where(signs[:-1] != signs[1:])[0]
    brackets = [(grid[i], grid[i+1], signs[i]) for i in idx]
    print("exact-scan sign changes: %d  (theta-count predicts ~%.0f)" %
          (len(brackets), (theta_v(TMAX)-theta_v(8.0))/np.pi), flush=True)
    with Pool(8) as p:
        zeros = np.array(sorted(p.map(refine, brackets, chunksize=4)))
    print("zeros: %d | first five: %s" % (len(zeros), np.round(zeros[:5], 5)))
    print("LMFDB truth:            9.22238 13.90755 17.44277 19.65651 22.33610")

    k = np.arange(1, len(zeros)+1) - 0.5
    S = k - theta_v(zeros)/np.pi
    print("S band: [%.3f, %.3f] drift %.3f  (PASS if width<2, drift<0.7)" %
          (S.min(), S.max(), abs(S[-10:].mean() - S[:10].mean())), flush=True)

    with Pool(8) as p:
        sg = np.array(p.map(sgn_slope, list(zeros), chunksize=4))
    alt = bool(np.all(sg[1:] == -sg[:-1]))
    print("sgn Lambda' strict alternation:", alt)
    ph = -np.pi/2 - theta_v(zeros) + np.pi*(sg < 0)
    ph = np.mod(ph + np.pi, 2*np.pi) - np.pi
    np.save("Delta_zeros_complete.npy", zeros)
    np.save("Delta_argLp_complete.npy", ph)
    print("saved Delta_zeros_complete.npy / Delta_argLp_complete.npy")
    print("DONE")
