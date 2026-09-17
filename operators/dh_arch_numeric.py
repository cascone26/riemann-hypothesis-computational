#!/usr/bin/env python3
"""
Numerical archimedean block for the truncated Weil form, from Groskin eq (arch-source):
   psi_R(x) = (1/2pi^2) INT_{-inf}^{inf} h_+(r) S(r,x,L) dr,
   S(r,x,L) = INT_0^L sin(2 pi x (1 - y/L)) cos(r y) dy = [cos(rL) - cos(2 pi x)] * B/(B^2 - r^2), B=2 pi x/L,
   h_+(r)  = Re psi_Gamma(nu + i r/2) - log(pi / cond),   nu=1/4 (zeta, even), nu=3/4 (D-H odd char mod 5), cond=1 (zeta) or 5 (D-H).
(Q_arch)_{mn} = divided difference of psi_R:  offdiag (psi_R(m)-psi_R(n))/(m-n),  diag psi_R'(n).

WHY numeric: for entire D-H (no pole -> DROP the W02=Q_pole block) the archimedean CONSTANTS kappa,J change with
nu+conductor; rather than hand-derive them (error-prone) we evaluate psi_R directly and VALIDATE against Groskin's
certified closed-form WR for zeta (nu=1/4). If it matches, the same machinery at nu=3/4 gives the D-H arch block.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "external", "groskin_2607.02828", "code"))
import mpmath as mp
mp.mp.dps = 30
from flint import arb, ctx
import arb_ldlt_certify as G


def make_psiR(Lval, nu, cond):
    """Return psi_R(x) as an mpmath function, for gamma-shift nu and conductor cond."""
    L = mp.mpf(Lval)
    logconst = mp.log(mp.pi / cond)      # h_+ constant term  - log(pi/cond)
    nu = mp.mpf(nu)

    def hplus(r):
        return mp.re(mp.digamma(nu + 1j*r/2)) - logconst

    def Sr(r, x):
        B = 2*mp.pi*x/L
        # limit-safe near r^2 = B^2
        denom = B*B - r*r
        cosrL = mp.cos(r*L)
        num = (cosrL - mp.cos(2*mp.pi*x)) * B
        if abs(denom) < mp.mpf(10)**(-20):
            # l'Hopital in r at r->|B|: use series; rare on the grid, nudge r
            r = r + mp.mpf(10)**(-15)
            denom = B*B - r*r
            cosrL = mp.cos(r*L); num = (cosrL - mp.cos(2*mp.pi*x))*B
        return num/denom

    period = mp.pi/L                      # cos(rL) oscillation period in r
    def psiR(x):
        x = mp.mpf(x)
        f = lambda r: hplus(r)*Sr(r, x)
        B = 2*mp.pi*abs(x)/L if x != 0 else mp.mpf(1)
        R = mp.mpf(50)                     # integrand ~ log(r)/r^2, tail past 50 ~1e-3 (fine for validation)
        pts = [mp.mpf(0)]
        k = 1
        while k*period < R:                # breakpoint every half-period of cos(rL)
            pts.append(k*period)
            k += 1
        pts.append(R)
        allpts = sorted(set([float(p) for p in pts] + [float(B)]))
        val = 2*mp.quad(f, allpts)
        # analytic tail INT_R^inf: S ~ B[cos(2pi x)-cos(rL)]/r^2; the non-oscillatory part ~ B cos(2pi x)/r^2
        # gives a LINEAR-in-x term that must NOT be dropped (else a*x error -> constant added to every entry).
        tail_int = mp.quad(lambda r: hplus(r)/r**2, [R, 2*R, 10*R, mp.inf])
        val += 2*B*mp.cos(2*mp.pi*x)*tail_int    # oscillatory cos(rL) tail is O(1/R^3), negligible
        return -val/(2*mp.pi**2)           # sign: tau = W02 - WR; psi_R convention = -WR
    return psiR


def zeta_closed_WR(n, m, N, c, prec=400):
    """Groskin's certified WR entry (archimedean block) for zeta, for validation."""
    ctx.prec = prec
    S, CC, XC, L = G.arb_closed_forms(N, c, prec)
    kappa = G.arb_kappa(L); J = G.arb_J(L); PI = arb.pi()
    def Ss(k): return S[k] if k >= 0 else -S[-k]
    if n == m:
        return float((kappa + 2*CC[abs(n)] + J - (2/L)*XC[abs(n)]).mid())
    return float(((Ss(m)-Ss(n))/(PI*(n-m))).mid())


if __name__ == "__main__":
    c = 14; N = 6
    Lval = float(arb(c).log().mid()) if False else float(mp.log(c))
    print(f"VALIDATE numeric psi_R against Groskin certified WR (zeta nu=1/4, cond=1), c={c}, L={Lval:.5f}")
    psiR = make_psiR(Lval, mp.mpf(1)/4, 1)
    # off-diagonal check: (psi_R(m)-psi_R(n))/(m-n) vs closed WR(m,n)
    print("  OFF-DIAGONAL (numeric divided-diff  vs  Groskin closed WR):")
    vals = {}
    for k in range(0, N+1):
        vals[k] = psiR(k)
    for (m, n) in [(1, 0), (2, 1), (3, 1), (4, 2), (5, 3)]:
        num = float((vals[m]-vals[n])/(m-n))
        clo = zeta_closed_WR(n, m, N, c)
        print(f"    (m={m},n={n}) numeric={num:+.6e}  closed={clo:+.6e}  rel.err={abs(num-clo)/(abs(clo)+1e-30):.2e}")
    # diagonal check via central difference of psi_R
    print("  DIAGONAL (numeric psi_R'(n)  vs  Groskin closed WR(n,n)):")
    h = mp.mpf(10)**(-6)
    for n in [1, 2, 3]:
        d = float((psiR(n+h)-psiR(n-h))/(2*h))
        clo = zeta_closed_WR(n, n, N, c)
        print(f"    n={n}  numeric psi_R'={d:+.6e}  closed diag={clo:+.6e}  rel.err={abs(d-clo)/(abs(clo)+1e-30):.2e}")
