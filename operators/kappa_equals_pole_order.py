#!/usr/bin/env python3
"""§958 SHARP TEST OF THE FRAMING: κ(arch−prime) = ORDER OF THE POLE, across the Selberg class.

Central claim (§947-954): the Weil form's one negative square is ζ's s=1 POLE; RH ⟺ pole absorbs it. Prediction
that would CONFIRM this cleanly: for a POLE-FREE L-function (non-principal Dirichlet χ: L(s,χ) is ENTIRE, no pole),
the arch−prime form should be ALREADY PSD — κ=0, nothing to lift — whereas ζ (simple pole) has κ=1. i.e.
    κ(arch−prime for L) = ord_{s=1} of the pole of L   (1 for ζ, 0 for entire χ/Δ).
If κ tracks the pole order across characters, the 'one negative square = the pole' identification is confirmed
NON-trivially (a NEW cross-L test, distinct from §928's razor-universality). Build the Weil form for real
(quadratic) Dirichlet χ: completed Λ(s,χ)=(q/π)^{(s+a)/2}Γ((s+a)/2)L(s,χ), a=parity(χ). arch density at s=1/2+ir:
Ω_χ(r)=Re ψ((a+1/2)/2 + i r/2) + ln(q/π); prime = (1/π)Σ Λ(n)χ(n)/√n·cos(r log n); NO pole term (χ non-principal).
Compare κ(arch−prime) for ζ vs several χ. HANDLE (own code)."""
import numpy as np
from scipy.special import digamma

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam
NMAX = 8000; LAM = vonmangoldt(NMAX)
NS = np.arange(2, NMAX+1); LOGN = np.log(NS); LAMv = LAM[2:]

def legendre_char(q):
    """real primitive quadratic character mod q (Kronecker symbol (·/q)) for prime q≡? — use Jacobi symbol."""
    def jacobi(a, n):
        a %= n; result = 1
        while a != 0:
            while a % 2 == 0:
                a //= 2
                if n % 8 in (3, 5): result = -result
            a, n = n, a
            if a % 4 == 3 and n % 4 == 3: result = -result
            a %= n
        return result if n == 1 else 0
    return np.array([jacobi(n, q) for n in range(0, NMAX+1)])

def weil_form(chi, q, a_par, cen_lo, cen_hi, fourW, frac, is_zeta=False, dr=0.004):
    W = fourW/4; L = cen_hi-cen_lo; Ncap = L*2*W/np.pi; N = max(4, round(frac*Ncap))
    RGRID = np.arange(-(cen_hi+120), (cen_hi+120), dr)
    # archimedean density: Re psi((a+1/2)/2 + i r/2) + ln(q/pi)
    OMEGA = np.real(digamma((a_par+0.5)/2 + 0.5j*RGRID)) + np.log(q/np.pi)
    CEN = list(np.linspace(cen_lo, cen_hi, N))
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    P = np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]
    wts = (LAMv[mask]*(chi[NS[mask]] if not is_zeta else 1.0))/np.sqrt(NS[mask])
    cosu = np.cos(np.outer(ln, RGRID)); prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*wts[None, :]).sum(axis=1)
    prime = 0.5*(prime+prime.T)
    return arch - prime, G, N

def kappa(M, G):
    gw, gv = np.linalg.eigh(G); keep = gw > gw.max()*1e-12
    B = gv[:, keep]/np.sqrt(gw[keep]); ev = np.linalg.eigvalsh(B.T @ M @ B)
    s = max(1.0, np.abs(ev).max()); return int(np.sum(ev < -1e-9*s)), ev.min()

# targets: (name, q, parity a, is_zeta)  — quadratic characters are real; parity a: χ(-1)=(-1)^a
# χ_3 = (·/3): χ(-1)=χ(2)=-1 → odd (a=1). χ_4: odd (a=1). χ_5=(·/5): χ(-1)=1 → even (a=0).
# χ_8=(2/·): even (a=0). χ_7: χ(-1)=(-1/7)=-1 → odd (a=1).
targets = [("zeta", 1, 0, True), ("chi_3(odd)", 3, 1, False), ("chi_4(odd)", 4, 1, False),
           ("chi_5(even)", 5, 0, False), ("chi_7(odd)", 7, 1, False), ("chi_8(even)", 8, 0, False),
           ("chi_11(odd)", 11, 1, False)]
print("="*92)
print("§958 κ(arch−prime) vs pole order across the Selberg class. Predict: ζ→1 (pole), entire χ→0.")
print("="*92)
print(f"{'L-function':>14} {'q':>3} {'a':>2} {'pole ord':>8} | {'4W':>4} {'N':>3} {'κ=#neg(arch−prime)':>18} {'min eig':>11}")
for fourW in [4.0, 5.0]:
    for (name, q, a_par, isz) in targets:
        chi = None if isz else legendre_char(q)
        M, G, N = weil_form(chi, q, a_par, 14, 44, fourW, 0.8, is_zeta=isz)
        k, mn = kappa(M, G)
        pole = 1 if isz else 0
        flag = "✓" if k == pole else "✗ MISMATCH"
        print(f"{name:>14} {q:>3} {a_par:>2} {pole:>8} | {fourW:>4.1f} {N:>3} {k:>18} {mn:>+11.3e}  {flag}")
    print("-"*92)
print("READ: if κ(arch−prime)=1 for ζ and =0 for every entire χ, then 'the one negative square IS the s=1 pole'")
print("is confirmed across the Selberg class — a NEW test (κ = pole order), independent of §928 razor-universality.")
print("A χ with κ=1 (no pole to lift) would BREAK the framing — so this is a real falsification test. Honest either way.")
print("="*92)
