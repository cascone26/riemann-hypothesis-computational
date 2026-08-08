#!/usr/bin/env python3
"""§835: quantify the MISALIGNMENT-PROTECTED positivity (§834's RH-heart). Past the razor the prime
top-eigenvalue exceeds the arch floor, yet the Weil margin stays 0⁺ — because the prime's top
eigenvector v_p and the arch's bottom eigenvector v_a are not aligned. Measure, in the G-metric,
the angle between v_p and v_a across N, and the "protection" = margin − (floor − pnorm) (how much the
true margin beats the crude scalar bound). If protection rises as alignment falls → the positivity in
the collapse regime is carried by eigenvector misalignment (the non-scalar object = RH-hard). Float."""
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

def mats(cs, W):
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
    return pole+arch, prime, G

def gsolve_ext(M, G):
    """generalized eig (M,G): return (eigenvalues asc, eigenvectors columns G-orthonormal)."""
    w, V = eigh(M, G); return w, V

print("="*100)
print("§835 MISALIGNMENT probe — is the past-razor positivity carried by v_prime ⟂ v_arch?")
print("angle = G-inner-product angle between prime-top eigvec and arch-bottom eigvec")
print("protection = margin − (floor − pnorm)  (how much true margin beats the crude scalar bound)")
print("="*100)
for fourW in [2.20]:
    W = fourW/4
    print(f"4W={fourW}, window [14,44]:")
    print(f"  {'N':>4} {'floor':>9} {'pnorm':>9} {'floor-pnorm':>12} {'margin':>12} {'protection':>11} {'angle(deg)':>10} {'cos':>7}")
    for N in range(5, 10):
        CEN = list(np.linspace(14, 44, N)); AP, PR, G = mats(CEN, W)
        wf, Vf = gsolve_ext(AP, G); floor = float(wf[0]); v_a = Vf[:, 0]        # arch bottom
        wp, Vp = gsolve_ext(PR, G); pnorm = float(wp[-1]); v_p = Vp[:, -1]      # prime top
        wm, _ = eigh(AP-PR, G); margin = float(wm[0])
        # G-inner product angle (vectors are G-orthonormal already: <v,v>_G=1)
        cosang = abs(float(v_a @ (G @ v_p)))
        ang = np.degrees(np.arccos(min(1.0, cosang)))
        protection = margin - (floor - pnorm)
        print(f"  {N:>4} {floor:>9.5f} {pnorm:>9.5f} {floor-pnorm:>12.5f} {margin:>+12.4e} {protection:>+11.4f} {ang:>10.2f} {cosang:>7.3f}")
print("="*100)
print("If, past the crossing (floor−pnorm<0), the angle is LARGE (v_p⟂v_a) and protection≈|floor−pnorm|")
print("(margin stays ~0⁺ while the scalar bound goes negative) → the residual positivity IS the")
print("eigenvector misalignment. That joint arch⊕prime geometry — not any scalar bound — is the RH-hard object.")
print("="*100)
