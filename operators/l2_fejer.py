#!/usr/bin/env python3
"""§376: SECOND INSTRUMENT — Fejer/sinc^2 test class (compact frequency support => FINITE exact
prime side). Tests whether the §335/§373/§375 phenomena are Gaussian-special or intrinsic.

Probe in r (=gamma) domain: F_c(r) = phi(r-c) + phi(r+c) (even), phi(x) = sinc^2(b x/pi) with
sinc(t)=sin(pi t)/(pi t), i.e. phi(x)=(sin(b x/2)/(b x/2))^2 ... use the clean pair:
  phi(x) = (sin(W x)/(W x))^2   [Fejer], FT_u phi = (1/W) tri(u/(2W)) supported on |u|<2W
    where tri(t)=max(0,1-|t|). So phi(x) = (W/ (2pi))? -- fix normalization by the RESIDUAL gate.
Product F_j F_k: its u-transform is a sum of shifted (phi*phi)-transforms; phi(x)^2-type products
of two Fejer kernels have FT = convolution of triangles = compactly supported on |u|<4W. So the
prime sum runs only over n with log n < 4W: FINITE, EXACT, no tail.

Sides (mirror the Gaussian instrument, VALIDATE by ||analytic - zero|| residual):
  zero:  2 sum_{gamma>0} F_j(gamma) F_k(gamma)
  arch:  (1/2pi) int F_j F_k Omega dr   (numerical, sinc^4 decay integrable)
  pole:  2 Re[F_j(i/2) F_k(i/2)]        (sinc^2 at imaginary arg = (sinh/.)^2, real)
  prime: 2 sum_n Lambda(n)/sqrt(n) gpair_fej(cj,ck,log n),  gpair_fej = (1/2pi) int F_jF_k cos(ur) dr
GATE: if ||full - zero||/||full|| is small (~1e-4 or better given numerical arch), instrument VALID.
"""
import sys
import numpy as np
from scipy.special import digamma

W = 1.2            # Fejer half-bandwidth in u (prime cutoff log n < 4W ~ 4.8 => n <= ~120)
BAND = (14.0, 44.0)

ZEROS = np.array([
    14.134725142, 21.022039639, 25.010857580, 30.424876126, 32.935061588,
    37.586178159, 40.918719012, 43.327073281, 48.005150881, 49.773832478,
    52.970321478, 56.446247697, 59.347044003, 60.831778525, 65.112544048,
    67.079810529, 69.546401711, 72.067157674, 75.704690699, 77.144840069,
    79.337375020, 82.910380854, 84.735492981, 87.425274613, 88.809111208,
    92.491899271, 94.651344041, 95.870634228, 98.831194218])

def vonmangoldt(N):
    lam = np.zeros(N + 1); isc = np.zeros(N + 1, dtype=bool)
    for p in range(2, N + 1):
        if not isc[p]:
            isc[2 * p::p] = True
            lp = np.log(p); pk = p
            while pk <= N:
                lam[pk] = lp; pk *= p
    return lam

NMAX = 200
LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX + 1)
LOGN = np.log(NS); LAMW = LAM[2:] / np.sqrt(NS)

RGRID = np.arange(-200.0, 200.0, 0.0025)
OMEGA = np.real(digamma(0.25 + 0.5j * RGRID)) - np.log(np.pi)

def fej(x):
    # (sin(Wx)/(Wx))^2, with the removable singularity at 0
    y = W * x
    out = np.ones_like(y)
    nz = np.abs(y) > 1e-12
    out[nz] = (np.sin(y[nz]) / y[nz]) ** 2
    return out

def probes(band, s):
    return np.arange(band[0], band[1] + 1e-9, s)

def Frow(c, r):
    return fej(r - c) + fej(r + c)

def build(cs, lam_scale=None):
    n = len(cs)
    dr = RGRID[1] - RGRID[0]
    P = np.array([Frow(c, RGRID) for c in cs])      # n x len(r)
    arch = (P * OMEGA[None, :]) @ P.T * dr / (2 * np.pi)
    G = P @ P.T * dr
    # prime side: gpair(cj,ck,u) = (1/2pi) int F_j F_k cos(u r) dr, numerical (exact-ish, fine grid)
    lw = LAMW if lam_scale is None else LAMW * lam_scale
    prime = np.zeros((n, n))
    cosu = np.cos(np.outer(LOGN, RGRID))            # nprime x len(r)
    for j in range(n):
        FjP = P[j][None, :] * P                      # n x len(r): F_j * F_k for all k
        # gpair for all k, all u: (1/2pi) int (F_jF_k)(r) cos(u r) dr
        gp = (FjP @ cosu.T) * dr / (2 * np.pi)       # n x nprime
        prime[j, :] = 2.0 * (gp * lw[None, :]).sum(axis=1)
    prime = 0.5 * (prime + prime.T)
    # pole: F_c(i/2) = fej(i/2 - c)+fej(i/2 + c); fej(x)= (sin(Wx)/(Wx))^2 complex
    def fej_c(x):
        y = W * x
        return (np.sin(y) / y) ** 2
    z = 0.5j
    fp = fej_c(z - cs) + fej_c(z + cs)
    pole = 2.0 * np.real(np.outer(fp, fp))
    FZ = np.array([Frow(c, ZEROS) for c in cs])
    bzero = 2.0 * FZ @ FZ.T
    full = pole + arch - prime
    return dict(full=full, archpole=pole + arch, overlap=G, prime=prime, zero=bzero)

def mineig(M): return float(np.linalg.eigvalsh(M)[0])

if __name__ == "__main__":
    print(f"FEJER INSTRUMENT (W={W}, prime cutoff log n<{4*W:.1f} => n<={int(np.exp(4*W))}). VALIDATION:")
    print(f"{'s':>5} {'n':>3} {'minEig(full)':>13} {'minEig(zero)':>13} {'resid':>10} {'rel':>9}")
    for s in (6.0, 5.0, 4.0, 3.0, 2.5, 2.0):
        cs = probes(BAND, s)
        f = build(cs)
        resid = np.max(np.abs(f['full'] - f['zero']))
        rel = resid / (np.max(np.abs(f['full'])) + 1e-300)
        print(f"{s:5.1f} {len(cs):3d} {mineig(f['full']):+13.4e} {mineig(f['zero']):+13.4e} "
              f"{resid:10.2e} {rel:9.2e}", flush=True)
    print("GATE: rel residual small => Fejer instrument VALID (analytic side = zero side).")
