#!/usr/bin/env python3
"""§846: verify the §844 lead before a workflow. The 'lift' = margin − min_window S(r) (how far
band-limiting raises λmin above the always-negative Weil symbol-min). Is the lift a PROLATE/
concentration quantity (tying §844 to CCM's prolate program) or something else? Compute, per (window,4W):
  S_min = min_window[(1/2π)Ω − P]   (symbol min, §844)
  margin = λmin(Q,G)                 (true Weil margin)
  lift = margin − S_min              (band-limiting gain)
  |S_neg|_L1 = window-integral of max(0,−S)   (symbol negative mass)
  conc_top = top eigenvalue of concentration A=2∫_win F F (in G-metric) — the prolate concentration factor
and check which the lift tracks: lift ≈ |S_neg| (extremal), lift ≈ (1−conc)·something (prolate), or neither."""
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

def analyze(a, b, W, dr=0.004):
    fourW = 4*W; L = b-a
    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]
    Psym = np.zeros_like(RGRID)
    for u, w in zip(ln, lw): Psym += w*np.cos(u*RGRID)
    Psym *= (1.0/np.pi)
    S = OMEGA/(2*np.pi) - Psym
    win = (RGRID >= a) & (RGRID <= b)
    S_min = float(S[win].min())
    S_neg_L1 = float(np.maximum(0, -S[win]).sum()*dr) / L    # avg negative depth
    N = max(4, round(0.75*L*2*W/np.pi)); CEN = list(np.linspace(a, b, N))
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    Pm = np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G = (Pm @ Pm.T)*dr
    arch = (Pm*OMEGA[None, :]) @ Pm.T * dr/(2*np.pi)
    cosu = np.cos(np.outer(ln, RGRID)); prime = np.zeros((N, N))
    for j in range(N):
        gp = ((Pm[j][None, :]*Pm) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*lw[None, :]).sum(axis=1)
    prime = 0.5*(prime+prime.T)
    def fc(c):
        y = W*(0.5j-c); y2 = W*(0.5j+c); return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in CEN]); pole = 2.0*np.real(np.outer(fp, fp))
    margin = float(eigh(pole+arch-prime, G, eigvals_only=True)[0])
    # concentration operator A = 2 int_window F_j F_k
    Aconc = (Pm[:, win] @ Pm[:, win].T)*dr*2.0
    conc_top = float(eigh(Aconc, G, eigvals_only=True)[-1])   # top prolate concentration eigenvalue
    lift = margin - S_min
    return S_min, S_neg_L1, margin, lift, conc_top, N

print("="*112)
print("§846 lift structure: is lift=margin−S_min a prolate/extremal quantity?")
print("="*112)
print(f"{'window':>12} {'4W':>4} {'S_min':>8} {'|Sneg|avg':>9} {'margin':>12} {'lift':>8} {'conc_top':>8} {'lift/|Smin|':>11}")
for (a, b, fourW) in [(14,44,2.0),(14,44,3.0),(14,44,4.0),(100,130,3.0),(100,130,4.0),(400,430,4.0)]:
    W = fourW/4
    Smin, Sneg, mg, lift, ct, N = analyze(a, b, W)
    print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {Smin:>8.4f} {Sneg:>9.4f} {mg:>+12.4e} {lift:>8.4f} {ct:>8.4f} {lift/abs(Smin):>11.4f}")
print("="*112)
print("lift/|S_min| ~ 1 → margin is the small residual of near-cancellation (RH = lift≥|S_min|).")
print("conc_top ~ 1 → probes fully concentrated in window (prolate). Track how lift relates to conc_top and |Sneg|.")
print("="*112)
