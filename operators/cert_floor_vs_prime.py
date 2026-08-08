#!/usr/bin/env python3
"""§834: decompose the knife-edge margin into arch-FLOOR vs PRIME norm (the §833 structural heart).
For spans around the collapse, in the G-metric (generalized eigenproblem), measure:
  floor  = λ_min(pole+arch, G)         — the archimedean PD floor (should stay ≳0.1, §816)
  pnorm  = λ_max(prime, G)             — the prime form's top eigenvalue (indefinite, eats the floor)
  margin = λ_min(pole+arch−prime, G)   — the full Weil margin
Question: is margin ≈ floor − pnorm (the assembly bound m ≥ floor − pnorm), and does the COLLAPSE
happen when pnorm rises to meet floor? If so, the certifiable Handle is m(window) ≥ floor − pnorm with
§816-CERT (floor) + §803A (prime) constants. Float, fast."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam
NMAX = 400; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)
RGRID = np.arange(-260.0, 260.0, 0.0025); dr = 0.0025
OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)

def fej(x, W):
    y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
    out[nz] = (np.sin(y[nz])/y[nz])**2; return out
def Frow(c, r, W): return fej(r-c, W) + fej(r+c, W)

def parts(cs, W):
    n = len(cs); fourW = 4*W
    P = np.array([Frow(c, RGRID, W) for c in cs]); G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]
    prime = np.zeros((n, n))
    if ln.size:
        cosu = np.cos(np.outer(ln, RGRID))
        for j in range(n):
            gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
            prime[j, :] = 2.0*(gp*lw[None, :]).sum(axis=1)
        prime = 0.5*(prime+prime.T)
    def fc(c):
        y = W*(0.5j-c); y2 = W*(0.5j+c)
        return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in cs]); pole = 2.0*np.real(np.outer(fp, fp))
    archpole = pole + arch
    floor = float(eigh(archpole, G, eigvals_only=True)[0])          # arch PD floor
    pnorm = float(eigh(prime, G, eigvals_only=True)[-1])            # prime top eigenvalue
    margin = float(eigh(archpole - prime, G, eigvals_only=True)[0]) # full Weil margin
    return floor, pnorm, margin

print("="*96)
print("§834 FLOOR vs PRIME decomposition — is margin ≈ floor − pnorm, and does collapse = pnorm→floor?")
print("floor=λmin(arch+pole,G)  pnorm=λmax(prime,G)  margin=λmin(Weil,G)")
print("="*96)
for fourW in [1.90, 2.20]:
    W = fourW/4
    print("-"*96)
    print(f"4W={fourW}, window [14,44], sweep N:")
    print(f"  {'N':>4} {'floor':>11} {'pnorm':>11} {'floor-pnorm':>12} {'margin':>13} {'margin≈floor-pnorm?':>20}")
    for N in range(5, 10):
        CEN = list(np.linspace(14, 44, N)); fl, pn, m = parts(CEN, W)
        approx = fl - pn
        close = "yes" if abs(approx - m) < 0.3*max(abs(m), 1e-6)+0.02 else "no"
        print(f"  {N:>4} {fl:>11.5f} {pn:>11.5f} {approx:>12.5f} {m:>+13.5e} {close:>20}")
print("="*96)
print("If margin tracks floor−pnorm and collapse coincides with pnorm rising to floor → the certifiable")
print("Handle is m ≥ floor − pnorm, using §816-CERT (floor≥~0.15) + §803A (prime bound). floor is the")
print("arch density = zero density (§833); pnorm is the prime Euler side. Their crossing IS the razor.")
print("="*96)
