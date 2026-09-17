#!/usr/bin/env python3
"""§952 IS THE ONE NEGATIVE SQUARE = THE POLE DIRECTION? (makes the Krein-Langer decomposition EXPLICIT)

§889/§947-949: M_σ=arch−prime has exactly ONE negative square; Q=M_σ+pole (rank-≤2 PSD) is PSD ⟺ RH; §948 saw
the negative eigenvector 'rotate' window-to-window. Hypothesis: it rotates BECAUSE the pole vector fp (which
depends on the window centers) rotates with it — the single negative direction of M_σ IS the pole direction.
If overlap(neg-eigvec(M_σ), pole-range) ≈ 1 uniformly, then:
  • the Krein-Langer decomposition is EXPLICIT + UNIFORM: K₋ = projector on the pole direction (known in closed
    form, NOT extracted from Q) — no circularity through 'Q is PSD';
  • RH-corner ⟺ pole absorbs exactly that direction ⟺ the scalar criterion p·(v·w)² ≥ μ (§891);
  • 'κ(M_σ)=1' becomes 'arch−prime is negative ONLY along the pole direction' — a structural, testable claim.
Measure: (a) cos-angle between M_σ's most-negative eigenvector w and the pole's 2-dim range span{Re fp, Im fp};
(b) how much of w lies in that pole range (projection norm). Across sliding windows + bandwidths. HANDLE."""
import numpy as np
from scipy.special import digamma

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam
NMAX = 8000; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)

def mats(a, b, W, N, dr=0.004):
    fourW = 4*W
    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    CEN = list(np.linspace(a, b, N))
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    P = np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]
    cosu = np.cos(np.outer(ln, RGRID)); prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*lw[None, :]).sum(axis=1)
    prime = 0.5*(prime+prime.T)
    def fc(c):
        y = W*(0.5j-c); y2 = W*(0.5j+c); return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in CEN])   # complex pole vector; pole = 2 Re(fp fp^T)
    return arch, prime, G, fp

def neg_eig_and_pole_overlap(a, b, W, N):
    arch, prime, G, fp = mats(a, b, W, N)
    gw, gv = np.linalg.eigh(G); keep = gw > gw.max()*1e-12
    B = gv[:, keep]/np.sqrt(gw[keep])                 # G^{-1/2}
    Msig = B.T @ (arch - prime) @ B                   # in G-orthonormal coords
    ev, evec = np.linalg.eigh(Msig)
    w = evec[:, 0]                                    # most-negative eigenvector (already G-orthonormal coords)
    mineig = ev[0]
    # pole range in the SAME G-orthonormal coords: pole = 2 Re(fp fp^T); its range = span{Re(B^T G fp?), ...}
    # In G-metric, the pole quadratic form on coeff x is 2|fp·x|^2 (since fp is a linear functional on coeffs).
    # Convert fp (a covector in coeff space) to the G-orthonormal coords y where x = B y: fp·x = (fp·B) y.
    fB = fp @ B                                        # complex length-keep vector: pole form = 2|fB·y|^2
    re, im = np.real(fB), np.imag(fB)
    # orthonormal basis of the (<=2)-dim real pole range
    M = np.vstack([re, im]); Q, _ = np.linalg.qr(M.T)  # columns span the pole range
    proj = Q @ (Q.T @ w)
    overlap = np.linalg.norm(proj)/np.linalg.norm(w)   # fraction of w inside the pole range (1 = fully aligned)
    # also the pole's PSD lift value along w: 2|fB·w|^2 vs |mineig|
    lift = 2*abs(np.dot(fB, w))**2
    return mineig, overlap, lift

print("="*96)
print("§952 Does M_σ's single negative direction lie in the POLE range? overlap→1 ⟹ explicit Krein-Langer K₋.")
print("="*96)
print(f"{'window':>12} {'4W':>4} {'N':>3} | {'min eig(M_σ)':>13} {'overlap w/ pole range':>21} {'pole lift 2|fp·w|²':>18} {'lift≥|μ|?':>9}")
cfgs = [(14,44,3.0),(14,44,4.0),(14,44,5.0),(14,44,6.0),(15,45,4.0),(16,46,4.0),(18,48,4.0),
        (14,34,4.0),(60,90,5.0),(120,150,5.0)]
for (a,b,fourW) in cfgs:
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,round(0.8*Ncap))
    mineig, ov, lift = neg_eig_and_pole_overlap(a,b,W,N)
    ok = "YES" if lift >= abs(min(mineig,0))*0.999 else "no"
    tag = "" if mineig < -1e-11 else "  (M_σ⪰0 here, no neg dir)"
    print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {N:>3} | {mineig:>+13.3e} {ov:>21.3f} {lift:>18.3e} {ok:>9}{tag}")
print("="*96)
print("READ: overlap≈1 ⟹ the one negative square IS the pole direction ⟹ Krein-Langer K₋ is explicit/uniform")
print("(closed-form pole covector fp), not extracted from Q ⟹ decomposition is non-circular. And lift≥|μ| along w")
print("⟺ the pole absorbs it ⟺ RH-corner (scalar criterion §891). Low overlap ⟹ neg dir is NOT the pole (the")
print("κ=1 miracle is deeper than 'missing pole') — also decisive, redirects the structural search. Honest either way.")
print("="*96)
