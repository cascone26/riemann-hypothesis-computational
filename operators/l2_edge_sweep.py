#!/usr/bin/env python3
"""§335 instrument: is the dense-regime certificate collapse (§286-287) an L2/INTEGRABILITY-threshold
phenomenon (sampling/frame/prolate plunge) rather than a SIGN phenomenon?

Roundtable-#4 lens (LEDGER §334): the complementary series is unitary => positivity frames are blind
at the tempered edge; the true threshold is an integrability/decay-rate condition with no slack.
Test on OUR instrument: the Weil form on symmetrized Gaussian probes (band [14,44], w=3 baseline,
exactly the §286-287 setup).

Closed forms (LEDGER §284): probe F_j(r) = e^{-(r-c_j)^2/2w^2} + e^{-(r+c_j)^2/2w^2} (even).
  B_full = pole + arch - prime   (analytic side)
  B_zero[j,k] = 2 Sum_{gamma>0} F_j(gamma) F_k(gamma)   (zero side; explicit formula => equal)
Self-validation: ||B_full - B_zero|| is the correctness certificate (core/weil.py discipline).

Discriminators:
  G_overlap[j,k] = int F_j F_k dr            -- pure L2 Gram, no arithmetic (sign-blind control 1)
  B_archpole     = pole + arch               -- no primes            (sign-blind control 2)
  B_fake         = pole + arch - prime_bent  -- Lambda(2) inflated x1.05 (t11-style sign probe)

Modes: validate | sweepA | sweepC | fake | spectrum
"""
import sys, os
import numpy as np
from scipy.special import digamma

# ---------------- shared machinery ----------------
W_DEFAULT = 3.0
BAND = (14.0, 44.0)

# first 29 nontrivial zeros (gamma < 101), standard values
ZEROS = np.array([
    14.134725142, 21.022039639, 25.010857580, 30.424876126, 32.935061588,
    37.586178159, 40.918719012, 43.327073281, 48.005150881, 49.773832478,
    52.970321478, 56.446247697, 59.347044003, 60.831778525, 65.112544048,
    67.079810529, 69.546401711, 72.067157674, 75.704690699, 77.144840069,
    79.337375020, 82.910380854, 84.735492981, 87.425274613, 88.809111208,
    92.491899271, 94.651344041, 95.870634228, 98.831194218])

def vonmangoldt(N):
    lam = np.zeros(N + 1)
    is_comp = np.zeros(N + 1, dtype=bool)
    for p in range(2, N + 1):
        if not is_comp[p]:                    # p is prime
            is_comp[2 * p::p] = True
            lp = np.log(p); pk = p
            while pk <= N:
                lam[pk] = lp; pk *= p
    return lam

NMAX = 200
LAM = vonmangoldt(NMAX)
NS = np.arange(2, NMAX + 1)
LOGN = np.log(NS)
LAMW = LAM[2:] / np.sqrt(NS)          # Lambda(n)/sqrt(n)

RGRID = np.arange(-90.0, 90.0, 0.002)
OMEGA = np.real(digamma(0.25 + 0.5j * RGRID)) - np.log(np.pi)

def probes(band, s):
    c0, c1 = band
    return np.arange(c0, c1 + 1e-9, s)

def Fmat(cs, w, r):
    """probe matrix P[j, i] = F_j(r_i), symmetrized Gaussians."""
    d1 = (r[None, :] - cs[:, None]) / w
    d2 = (r[None, :] + cs[:, None]) / w
    return np.exp(-0.5 * d1**2) + np.exp(-0.5 * d2**2)

def g_pair(cj, ck, w, u):
    """(1/2pi) int F_j F_k cos(ur) dr, closed form (LEDGER §284 generalized off-diagonal).
    u may be an array."""
    pref = w * np.sqrt(np.pi) / (2 * np.pi)
    e_near = np.exp(-(cj - ck)**2 / (4 * w**2))
    e_far  = np.exp(-(cj + ck)**2 / (4 * w**2))
    gauss = np.exp(-u**2 * w**2 / 4.0)
    return pref * gauss * 2.0 * (e_near * np.cos(u * (cj + ck) / 2.0)
                                 + e_far * np.cos(u * (cj - ck) / 2.0))

def build_forms(cs, w, lam_scale=None):
    """returns dict of matrices. lam_scale: optional array over NS multiplying LAMW (for bent primes)."""
    n = len(cs)
    P = Fmat(cs, w, RGRID)
    dr = RGRID[1] - RGRID[0]
    arch = (P * OMEGA[None, :]) @ P.T * dr / (2 * np.pi)
    G = P @ P.T * dr                                  # pure overlap Gram
    # pole term: h(i/2)+h(-i/2) = 2 Re[F_j(i/2) F_k(i/2)], F complex at imaginary argument
    z = 0.5j
    fp = np.exp(-((z - cs)**2) / (2 * w**2)) + np.exp(-((z + cs)**2) / (2 * w**2))
    pole = 2.0 * np.real(np.outer(fp, fp))
    # prime term
    lw = LAMW if lam_scale is None else LAMW * lam_scale
    prime = np.zeros((n, n))
    for j in range(n):
        for k in range(j, n):
            val = 2.0 * np.sum(lw * g_pair(cs[j], cs[k], w, LOGN))
            prime[j, k] = prime[k, j] = val
    full = pole + arch - prime
    # zero side
    FZ = Fmat(cs, w, ZEROS)
    bzero = 2.0 * FZ @ FZ.T
    return dict(full=full, archpole=pole + arch, overlap=G, prime=prime, zero=bzero)

def mineig(M):
    return float(np.linalg.eigvalsh(M)[0])

def diag_margin(M):
    d = np.diag(M)
    off = np.sum(np.abs(M), axis=1) - np.abs(d)
    return float(np.min(d - off))

def nrank(M, tol=1e-10):
    ev = np.linalg.eigvalsh(M)
    return int(np.sum(ev > tol * max(1.0, ev[-1])))

def zero_count(cs, w):
    """effective Weyl count: zeros within [min c - 2w, max c + 2w]."""
    lo, hi = cs.min() - 2 * w, cs.max() + 2 * w
    return int(np.sum((ZEROS >= lo) & (ZEROS <= hi)))

# ---------------- modes ----------------
def validate():
    print("VALIDATE vs banked §286-287 (band [14,44], w=3):")
    print("expect: s=6 margin +0.267 minEig 0.886 | s=5 minEig 0.241 | s=4 minEig ~8.6e-5 | s=3 ~1e-18")
    for s in (6.0, 5.0, 4.0, 3.0):
        cs = probes(BAND, s)
        f = build_forms(cs, W_DEFAULT)
        resid = np.max(np.abs(f['full'] - f['zero']))
        print(f"  s={s:.1f} n={len(cs):2d}  diagMargin={diag_margin(f['full']):+8.4f}  "
              f"minEig(full)={mineig(f['full']):+.4e}  ||full-zero||_max={resid:.2e}")
    print("residual is the correctness certificate: small => construction VALIDATED (explicit formula closes).")

def sweepA():
    print("SWEEP A: w=3, band [14,44], fine spacing. full vs sign-blind controls vs rank.")
    print(f"{'s':>5} {'n':>3} {'Nzero':>5} {'minEig(full)':>13} {'minEig(A+P)':>12} "
          f"{'minEig(ovl)':>12} {'rank(zero)':>10} {'resid':>9}")
    for s in np.arange(6.5, 2.49, -0.1):
        cs = probes(BAND, s)
        f = build_forms(cs, W_DEFAULT)
        resid = np.max(np.abs(f['full'] - f['zero']))
        print(f"{s:5.1f} {len(cs):3d} {zero_count(cs, W_DEFAULT):5d} "
              f"{mineig(f['full']):+13.4e} {mineig(f['archpole']):+12.4e} "
              f"{mineig(f['overlap']):+12.4e} {nrank(f['zero']):10d} {resid:9.1e}", flush=True)

def sweepC():
    print("SWEEP C: scaling collapse. For each w: minEig(full) vs probe count / zero count.")
    print(f"{'w':>4} {'s':>5} {'n':>3} {'Nz':>3} {'ratio':>6} {'minEig(full)':>13} {'minEig(ovl)':>12}")
    for w in (2.0, 2.5, 3.0, 3.5, 4.0):
        for s in np.arange(6.5, 2.49, -0.25):
            cs = probes(BAND, s)
            f = build_forms(cs, w)
            nz = zero_count(cs, w)
            print(f"{w:4.1f} {s:5.2f} {len(cs):3d} {nz:3d} {len(cs)/nz:6.3f} "
                  f"{mineig(f['full']):+13.4e} {mineig(f['overlap']):+12.4e}", flush=True)

def fake():
    print("FAKE-PRIME PROBE (t11-style): Lambda(2) x 1.05. Where does minEig go NEGATIVE vs spacing?")
    scale = np.ones_like(LAMW)
    # inflate all powers of 2
    for i, n in enumerate(NS):
        m = n
        while m % 2 == 0: m //= 2
        if m == 1: scale[i] = 1.05
    print(f"{'s':>5} {'n':>3} {'minEig(true)':>13} {'minEig(bent)':>13} {'delta':>10}")
    for s in np.arange(6.5, 2.49, -0.1):
        cs = probes(BAND, s)
        ft = build_forms(cs, W_DEFAULT)
        fb = build_forms(cs, W_DEFAULT, lam_scale=scale)
        mt, mb = mineig(ft['full']), mineig(fb['full'])
        print(f"{s:5.1f} {len(cs):3d} {mt:+13.4e} {mb:+13.4e} {mb-mt:+10.2e}", flush=True)

def spectrum():
    print("SPECTRUM: full sorted eigenvalues at several spacings (prolate-plunge shape check).")
    for s in (4.5, 4.0, 3.5, 3.0, 2.5):
        cs = probes(BAND, s)
        f = build_forms(cs, W_DEFAULT)
        ev = np.linalg.eigvalsh(f['full'])[::-1]
        nz = zero_count(cs, W_DEFAULT)
        print(f"s={s:.1f} n={len(cs)} Nzero={nz}: " + " ".join(f"{e:+.3e}" for e in ev), flush=True)

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "validate"
    {"validate": validate, "sweepA": sweepA, "sweepC": sweepC,
     "fake": fake, "spectrum": spectrum}[mode]()
