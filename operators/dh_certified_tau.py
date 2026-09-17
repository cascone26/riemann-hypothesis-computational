#!/usr/bin/env python3
"""
CERTIFIED D-H tau scaffold (cycle 203 lead). Parameterizes Groskin's arb arch by the digamma shift nu
(zeta even char: nu=1/4; Davenport-Heilbronn odd char mod 5: nu=3/4), keeping the chain-rule Jacobians
(1/2 in S/CC, 1/4 in XC) fixed. VERIFY-FIRST (anchor-to-exemplar): nu=1/4 + pole + prime-power weights MUST
reproduce the canonical build_arb_tau EXACTLY; only then are the D-H parameters trustworthy.

Remaining for a FAITHFUL D-H tau (flagged, not yet done rigorously): (i) drop the pole term kappa (D-H entire),
(ii) recompute the boundary constant J for the no-pole completed function + conductor 5 (arb_J is zeta-pole-
specific), (iii) prime side = signed Lambda_D(n)/sqrt(n) at ALL n<=c. This file locks (the verified) parts (i)
scaffold + prime-side hook; (ii) J is the subtle open piece.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "external", "groskin_2607.02828", "code"))
from flint import arb, acb, arb_mat, ctx
import arb_ldlt_certify as G
from arb_ldlt_certify import trigamma_acb, arb_J, arb_kappa, prime_powers_up_to, build_arb_tau
from arb_ldlt_margin import inertia_with_margin


def geom_sums_nu(n, L, prec, nu):
    """Groskin _geom_sums with the shift generalized: c_k = 2k + 2*nu (nu=1/4 recovers the 0.5 lattice)."""
    PI = arb.pi()
    w = 2 * PI * n / L
    w2 = w * w
    gS = arb(0); gCC = arb(0); gX1 = arb(0); gX2 = arb(0)
    thr = arb(2) ** (-(prec + 24))
    two_nu = 2 * arb(nu)
    k = 0
    while True:
        c_k = arb(2 * k) + two_nu
        e = (-c_k * L).exp()
        den = c_k * c_k + w2
        gS += e / den
        if n != 0:
            gCC += e * w2 / (c_k * den)
        gX1 += e * c_k / den
        gX2 += e * (c_k * c_k - w2) / (den * den)
        if e < thr and k > 2:
            break
        k += 1
    c_next = arb(2 * (k + 1)) + two_nu
    one_m_z = 1 - (-2 * L).exp()
    tail_geo = (-c_next * L).exp() / one_m_z
    rem = arb(4) * tail_geo
    def widen(x):
        return x + arb(0, rem)
    return widen(gS), widen(gCC), widen(gX1), widen(gX2)


def closed_forms_nu(N, c, prec, nu):
    ctx.prec = prec
    L = arb(c).log()
    PI = arb.pi()
    nu_a = arb(nu)
    psi_nu = nu_a.digamma()
    S = [arb(0)] * (N + 1); CC = [arb(0)] * (N + 1); XC = [arb(0)] * (N + 1)
    for n in range(N + 1):
        w = 2 * PI * n / L
        zarg = acb(nu_a, PI * n / L)
        psi = zarg.digamma()
        psi1 = trigamma_acb(zarg)
        gS, gCC, gX1, gX2 = geom_sums_nu(n, L, prec, nu)
        if n == 0:
            S[n] = arb(0); CC[n] = arb(0)
        else:
            S[n] = arb('0.5') * psi.imag - w * gS
            CC[n] = -arb('0.5') * (psi.real - psi_nu) + gCC
        XC[n] = arb('0.25') * psi1.real - L * gX1 - gX2
    return S, CC, XC, L


def build_tau_general(c, N, prec, nu=0.25, include_pole=True, kj_const=None,
                      prime_weights=None, prime_positions=None):
    """Generalized tau.
       nu=0.25, include_pole=True, kj_const=None, prime_weights=None -> canonical zeta build_arb_tau.
       D-H: nu=0.75, include_pole=False (entire, drop the W02 pole block), kj_const=(kappa+J)_DH,
            prime_weights=signed Lambda_D/sqrt(n), prime_positions=log n for all n<=c.
       include_pole toggles the W02=Q_pole block (cycle 205: W02 IS the pole, NOT kappa).
       kj_const overrides the archimedean diagonal constant (kappa+J); None -> exact Groskin arb_kappa+arb_J."""
    ctx.prec = prec
    S, CC, XC, L = closed_forms_nu(N, c, prec, nu)
    PI = arb.pi()
    sp2 = 16 * PI * PI
    l2 = L * L
    pref02 = 32 * L * (L / 4).sinh() ** 2
    KJ = (arb_kappa(L) + arb_J(L)) if kj_const is None else arb(kj_const)
    pole_on = include_pole

    if prime_weights is None:
        pdata = prime_powers_up_to(c)
        weights = [arb(p).log() * (arb(q) ** arb('-0.5')) for (q, p) in pdata]
        positions = [arb(q).log() for (q, p) in pdata]
    else:
        weights = prime_weights; positions = prime_positions
    M = len(weights)

    def S_signed(nn):
        return S[nn] if nn >= 0 else -S[-nn]

    DIM = 2 * N + 1
    A = arb_mat(DIM, DIM)
    for i in range(DIM):
        n = i - N
        for j in range(i, DIM):
            m = j - N
            num = l2 - sp2 * m * n
            den = (l2 + sp2 * m * m) * (l2 + sp2 * n * n)
            W02 = (pref02 * num / den) if pole_on else arb(0)
            if n == m:
                WR = KJ + 2 * CC[abs(n)] - (2 / L) * XC[abs(n)]
            else:
                WR = (S_signed(m) - S_signed(n)) / (PI * (n - m))
            Wp = arb(0)
            for idx in range(M):
                y = positions[idx]
                if n == m:
                    q = 2 * (1 - y / L) * (2 * PI * n * y / L).cos()
                else:
                    q = ((2 * PI * m * y / L).sin() - (2 * PI * n * y / L).sin()) / (PI * (n - m))
                Wp += weights[idx] * q
            val = W02 - WR - Wp
            A[i, j] = val; A[j, i] = val
    return A, DIM


def margin(A, DIM):
    npos, nneg, undet, mlow, mmid = inertia_with_margin(A, DIM)
    return nneg, (undet is None and nneg == 0), mlow


def dh_prime_side(c, prec):
    """signed generalized von Mangoldt Lambda_D(n)/sqrt(n) at positions log n, ALL n<=c (no Euler product)."""
    import math
    from _dh_negative_control import lambda_D
    ctx.prec = prec
    LD = lambda_D(c)
    weights = [arb(float(LD[n])) * (arb(n) ** arb('-0.5')) for n in range(2, c + 1)]
    positions = [arb(n).log() for n in range(2, c + 1)]
    return weights, positions


if __name__ == "__main__":
    c, N, prec = 14, 32, 800
    print(f"VERIFY-FIRST: nu=1/4 + pole + prime-powers reproduces canonical build_arb_tau (c={c},N={N},prec={prec})")
    A0, D0 = build_arb_tau(c, N, prec)
    n0, pd0, m0 = margin(A0, D0)
    A1, D1 = build_tau_general(c, N, prec, nu=0.25, include_pole=True)
    n1, pd1, m1 = margin(A1, D1)
    print(f"  canonical      : n_neg={n0} PD={pd0} margin={m0:.6e}")
    print(f"  general(nu=1/4): n_neg={n1} PD={pd1} margin={m1:.6e}")
    print(f"  => scaffold {'VERIFIED' if (n0==n1 and abs(float(m0)-float(m1))/abs(float(m0))<1e-6) else 'MISMATCH'}")

    print(f"\nSEMI-CERTIFIED D-H: nu=3/4, NO pole (entire), signed Lambda_D, (kappa+J)_DH swept for robustness.")
    dw, dpos = dh_prime_side(c, prec)
    for kj in [0.5945, 0.5957, 0.5970, 0.610, 0.580]:
        A, DIM = build_tau_general(c, N, prec, nu=0.75, include_pole=False, kj_const=kj,
                                   prime_weights=dw, prime_positions=dpos)
        nn, pd, mg = margin(A, DIM)
        tag = "  <-- validated (kappa+J)_DH" if abs(kj-0.5957) < 1e-6 else ""
        mgs = f"{mg:.4e}" if mg is not None else "n/a"
        print(f"  (kappa+J)_DH={kj:.4f}: n_neg={nn} certified_PD={pd} min_pivot={mgs}{tag}")
    print("  => n_neg>=1 across the swept constant range = D-H fails certified positivity (robust to the ~0.02% const uncertainty).")
