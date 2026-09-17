#!/usr/bin/env python3
"""§951 IS κ (negative squares) A MORSE/TOPOLOGICAL COUNT OF THE SYMBOL? (mechanism of the bounded-inertia)

§949: κ=#neg(M_σ) is arithmetic — REAL primes pin κ=1, perturbations grow κ. WHY does the real symbol give
exactly 1? Hypothesis: κ is a TOPOLOGICAL count of the symbol σ(r)=Ω/2π−P(r)'s negative structure within the
window — e.g. the number of DEEP negative wells the band-limited space can resolve. If κ ≈ (#negative wells of
σ weighted by width vs band resolution), then the bounded-inertia is an index/Morse phenomenon, and 'κ=1 for
real primes' becomes 'the real symbol has ONE resolvable negative well' — a provable structural statement about
σ = archimedean − von Mangoldt, independent of zeros. Test:
  - compute σ(r) on the window; find connected components of {σ<0} and their depths/widths;
  - band resolution ~ π/(2W) (Fejér cell). #wells wider than the cell = 'resolvable' count.
  - compare that predicted count to actual κ, for REAL vs perturbed/scrambled symbols.
Also report #neg(arch) and #neg(prime) separately (is either PSD? — grounds the Krein-Langer decomposition).
HANDLE (own code)."""
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

def symbol(RGRID, freqs, weights):
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    Pr = np.zeros_like(RGRID)
    for f, w in zip(freqs, weights):
        Pr += w*np.cos(RGRID*f)
    Pr /= np.pi
    return OMEGA/(2*np.pi) - Pr

def wells(sig, r, cell):
    """connected components of {sig<0}; return (n_components, n_resolvable) where resolvable = width>=cell."""
    neg = sig < 0
    comps = 0; resolvable = 0; i = 0; n = len(neg); dr = r[1]-r[0]
    while i < n:
        if neg[i]:
            j = i
            while j < n and neg[j]: j += 1
            comps += 1
            if (j-i)*dr >= cell: resolvable += 1
            i = j
        else:
            i += 1
    return comps, resolvable

def build_M(a, b, W, N, freqs, weights, dr=0.004):
    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    CEN = list(np.linspace(a, b, N))
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    P = np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
    cosu = np.cos(np.outer(freqs, RGRID)); prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*weights[None, :]).sum(axis=1)
    prime = 0.5*(prime+prime.T)
    return arch, prime, G

def nneg(M, G):
    gw, gv = np.linalg.eigh(G); keep = gw > gw.max()*1e-12
    B = gv[:, keep]/np.sqrt(gw[keep]); ev = np.linalg.eigvalsh(B.T @ M @ B)
    s = max(1.0, np.abs(ev).max()); return int(np.sum(ev < -1e-9*s))

rng = np.random.default_rng(11)
a, b = 14, 44
print("="*102)
print("§951 κ vs symbol negative-well topology. cell = π/(2W) (Fejér resolution). Window [14,44].")
print("="*102)
print(f"{'4W':>4} {'N':>3} {'variant':>12} | {'#neg(arch)':>10} {'#neg(prime)':>11} {'κ=#neg(M_σ)':>11} | {'σ<0 comps':>9} {'resolvable':>10}")
rfine = np.arange(a, b, 0.002)
for fourW in [3.0, 4.0, 5.0, 6.0]:
    W = fourW/4; L=b-a; Ncap = L*2*W/np.pi; N = max(4, round(0.8*Ncap)); cell = np.pi/(2*W)
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]; m = len(ln)
    variants = {
        'REAL': (ln, lw),
        'perturb1e-2': (ln*(1+1e-2*rng.uniform(-1,1,m)), lw),
        'scrambled': (rng.uniform(ln.min(), ln.max(), m), lw),
    }
    for name,(fr,wt) in variants.items():
        arch, prime, G = build_M(a,b,W,N,fr,wt)
        na, npr, k = nneg(arch,G), nneg(prime,G), nneg(arch-prime,G)
        sig = symbol(rfine, fr, wt); comps, res = wells(sig, rfine, cell)
        print(f"{fourW:>4.1f} {N:>3} {name:>12} | {na:>10} {npr:>11} {k:>11} | {comps:>9} {res:>10}")
print("="*102)
print("READ: if κ tracks 'resolvable' (#negative wells wider than the Fejér cell) — REAL symbol having ONE deep")
print("resolvable well while perturbations shatter σ into many — then κ is a MORSE/topological count of the")
print("symbol, and 'κ=1' is a provable structural fact about σ=archimedean−vonMangoldt (index-theorem realization).")
print("#neg(arch),#neg(prime): if arch⪰0 on-window but prime indefinite, the Krein-Langer split needs a smarter")
print("K₊/K₋ than arch/prime — the resolvable-well structure suggests where the rank-1 K₋ comes from.")
print("="*102)
