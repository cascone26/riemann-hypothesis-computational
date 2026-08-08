#!/usr/bin/env python3
"""
Find which N value at γ=69.5 gives the cos_max=0.999984 (the worst misalignment).
"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh

NMAX = 400
def vonmangoldt(N):
    lam = np.zeros(N+1)
    isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True
            lp = np.log(p)
            pk = p
            while pk <= N:
                lam[pk] = lp
                pk *= p
    return lam

LAM = vonmangoldt(NMAX)
NS = np.arange(2, NMAX+1)
LOGN = np.log(NS)
LAMW = LAM[2:]/np.sqrt(NS)
RGRID = np.arange(-320.0, 320.0, 0.0025)
dr = 0.0025
OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)

def fej(x, W):
    y = W*x
    out = np.ones_like(y)
    nz = np.abs(y) > 1e-12
    out[nz] = (np.sin(y[nz])/y[nz])**2
    return out

def Frow(c, r, W):
    return fej(r-c, W) + fej(r+c, W)

def mats(cs, W):
    n = len(cs)
    fourW = 4*W
    P = np.array([Frow(c, RGRID, W) for c in cs])
    G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)

    mask = LOGN < fourW
    ln = LOGN[mask]
    lw = LAMW[mask]
    prime = np.zeros((n, n))
    if ln.size:
        cosu = np.cos(np.outer(ln, RGRID))
        for j in range(n):
            gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
            prime[j, :] = 2.0*(gp*lw[None, :]).sum(axis=1)
        prime = 0.5*(prime+prime.T)

    def fc(c):
        y = W*(0.5j-c)
        y2 = W*(0.5j+c)
        return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in cs])
    pole = 2.0*np.real(np.outer(fp, fp))

    return pole+arch, prime, G

gamma = 69.5
window_length = 30
height_scale = max(1.0, np.log(max(gamma, 30) / 30.0))
fourW_target = 2.8 * (1.0 + height_scale * 0.3)
W = fourW_target / 4.0
Ncap = 30 * 2 * W / np.pi

hw = window_length / 2.0
a = gamma - hw
b = gamma + hw

print(f"Finding worst N at γ={gamma}")
print(f"Window [{a:.1f}, {b:.1f}], Ncap={Ncap:.2f}")
print()
print(f"{'N':>3} {'cos':>12} {'angle°':>8}")
print("-"*30)

worst_N = None
worst_cos = -1.0

for N in range(7, 17):
    CEN = list(np.linspace(a, b, N))
    AP, PR, G = mats(CEN, W)

    wf, Vf = eigh(AP, G)
    v_a = Vf[:, 0]

    wp, Vp = eigh(PR, G)
    v_p = Vp[:, -1]

    cosang = abs(float(v_a @ (G @ v_p)))
    ang = np.degrees(np.arccos(min(1.0, cosang)))

    print(f"{N:3d} {cosang:12.6f} {ang:8.2f}")

    if cosang > worst_cos:
        worst_cos = cosang
        worst_N = N

print("-"*30)
print(f"\nWorst N: {worst_N} with cos={worst_cos:.6f}")
