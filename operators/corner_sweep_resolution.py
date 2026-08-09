#!/usr/bin/env python3
"""§842: the corner-race workflow left the margin SIGN 'undetermined at 4W≥3' (agents saw ~1e-9 and called
it float noise). But §841 showed the 4W=4 margin is RESOLUTION-STABLE (real, not noise). Settle it across
the whole sweep: compute margin at [14,44] for 4W=2.5..5.5 at TWO integration resolutions (dr=0.004, 0.002).
If they agree to ≥3 sig figs → the margin is genuinely-positive (integration error < the margin) = float-
TRUSTWORTHY. Where they diverge → genuinely at the noise floor. This maps the trustworthy-positive range,
upgrading the agents' conservative 'undetermined' where it's actually resolved."""
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
NMAX = 3000; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)

def margin(cs, W, dr):
    RGRID = np.arange(-300.0, 300.0, dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    n = len(cs); fourW = 4*W
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    P = np.array([fej(RGRID-c)+fej(RGRID+c) for c in cs]); G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]
    prime = np.zeros((n, n)); cosu = np.cos(np.outer(ln, RGRID))
    for j in range(n):
        gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*lw[None, :]).sum(axis=1)
    prime = 0.5*(prime+prime.T)
    def fc(c):
        y = W*(0.5j-c); y2 = W*(0.5j+c); return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in cs]); pole = 2.0*np.real(np.outer(fp, fp))
    return float(eigh(pole+arch-prime, G, eigvals_only=True)[0])

print("="*88)
print("§842 corner margin resolution-stability sweep — [14,44], trustworthy-positive range")
print("="*88)
print(f"{'4W':>5} {'N':>4} {'margin(dr=.004)':>17} {'margin(dr=.002)':>17} {'agree?':>8} {'verdict':>22}")
for fourW in [2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]:
    W = fourW/4; Ncap = 15*fourW/np.pi; N = max(4, round(0.75*Ncap))
    CEN = list(np.linspace(14, 44, N))
    m1 = margin(CEN, W, 0.004); m2 = margin(CEN, W, 0.002)
    rel = abs(m1-m2)/max(abs(m1), 1e-300)
    agree = rel < 1e-3
    if m1 > 0 and m2 > 0 and agree: verdict = "TRUSTWORTHY +positive"
    elif m1 > 0 and m2 > 0: verdict = "positive, near noise"
    else: verdict = "AT NOISE FLOOR"
    print(f"{fourW:>5.1f} {N:>4} {m1:>+17.6e} {m2:>+17.6e} {str(agree):>8} {verdict:>22}", flush=True)
print("="*88)
print("agree (rel<1e-3) + both>0 → integration error < margin → genuinely positive (float-trustworthy).")
print("This settles the workflow's 'undetermined ≥3' where it's actually resolved.")
print("="*88)
