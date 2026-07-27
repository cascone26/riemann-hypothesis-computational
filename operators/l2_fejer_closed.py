#!/usr/bin/env python3
"""§381: CLOSED-FORM Fejer prime g_pair (triangle-FT convolution) — the §380 enabler.
fej(x)=(sin Wx/Wx)^2 has FT T(v)=(pi/W)(1-|v|/2W)_+ on |v|<2W (even triangle).
  int fej(r-a)fej(r-b) e^{iur} dr = (1/2pi) e^{iub} J(u, a-b),  J(u,d)=int T(v)T(u+v) e^{-ivd} dv.
gpair_closed(cj,ck,u) = (1/2pi) Re sum_{a in {cj,-cj}, b in {ck,-ck}} e^{iub} J(u,a-b).
J computed EXACTLY by piecewise (linear x linear x exp) integration over v where both triangles live.
VALIDATE: build the prime matrix from gpair_closed and compare to l2_fejer's NUMERICAL prime matrix.
If they match, the closed form is correct and the rigorous cert (§380) becomes feasible.
"""
import numpy as np
import l2_fejer as F

W = F.W

def T(v):
    return np.where(np.abs(v) < 2 * W, (np.pi / W) * (1 - np.abs(v) / (2 * W)), 0.0)

def J(u, d, NG=20000):
    # exact-ish via fine grid on the exact support (validation target; "closed form" = this piecewise
    # integral, evaluated to high precision). Support: v in [-2W,2W] AND u+v in [-2W,2W].
    lo = max(-2 * W, -2 * W - u); hi = min(2 * W, 2 * W - u)
    if hi <= lo: return 0.0 + 0.0j
    v = np.linspace(lo, hi, NG)
    integ = T(v) * T(u + v) * np.exp(-1j * v * d)
    return np.trapezoid(integ, v)

def gpair_closed(cj, ck, u):
    tot = 0.0 + 0.0j
    for a in (cj, -cj):
        for b in (ck, -ck):
            tot += np.exp(1j * u * b) * J(u, a - b)
    return np.real(tot) / (2 * np.pi) ** 2

def Jexact(u, d):
    """truly closed-form J via per-piece antiderivative of (alpha+beta v)(gamma+delta v)e^{-ivd}.
    breakpoints in v: -2W, 2W, -u-2W, -u+2W, and kinks of |v| (0) and |u+v| (-u). Sort within support."""
    lo = max(-2 * W, -2 * W - u); hi = min(2 * W, 2 * W - u)
    if hi <= lo: return 0.0 + 0.0j
    bpts = sorted(set([lo, hi] + [x for x in (0.0, -u) if lo < x < hi]))
    c = np.pi / W
    def antideriv(v, A, B, d):
        # int (A + B v) e^{-i d v} dv  (indefinite), d may be 0
        if abs(d) < 1e-14:
            return A * v + B * v * v / 2
        e = np.exp(-1j * d * v)
        # int (A+Bv)e^{-idv} = e^{-idv}[ (A+Bv)/(i d)?? ] -- do it stably:
        # int e^{-idv}dv = e^{-idv}/(-id); int v e^{-idv}dv = e^{-idv}(-v/(id) - 1/(d^2))... derive:
        # d/dv[e^{-idv}(a+bv)] = e^{-idv}(b - id(a+bv)). Solve for int (A+Bv)e^{-idv}:
        # try F=e^{-idv}(p+qv): F' = e^{-idv}(q - id(p+qv)) = e^{-idv}(-idq v + (q-idp))
        # match: -id q = B => q = B/(-id) = iB/d ; q - id p = A => p=(q-A)/(id)= (iB/d - A)/(id)
        q = 1j * B / d
        p = (q - A) / (1j * d)
        return e * (p + q * v)
    tot = 0.0 + 0.0j
    for lft, rgt in zip(bpts[:-1], bpts[1:]):
        mid = (lft + rgt) / 2
        # on this piece T(v)=c(1 - |v|/2W): sign of v fixed
        sv = 1.0 if mid > 0 else -1.0
        A1 = c * 1.0; B1 = -c * sv / (2 * W)   # T(v) = A1 + B1 v
        # T(u+v)=c(1-|u+v|/2W): sign of u+v fixed
        su = 1.0 if (u + mid) > 0 else -1.0
        A2 = c * (1 - su * u / (2 * W)); B2 = -c * su / (2 * W)  # T(u+v)=A2 + B2 v
        # product (A1+B1 v)(A2+B2 v) = P0 + P1 v + P2 v^2
        P0 = A1 * A2; P1 = A1 * B2 + B1 * A2; P2 = B1 * B2
        # int (P0+P1 v+P2 v^2) e^{-idv} dv over [lft,rgt]: split; reuse linear helper + quadratic term
        def antiquad(v, d):
            if abs(d) < 1e-14:
                return P0 * v + P1 * v * v / 2 + P2 * v ** 3 / 3
            e = np.exp(-1j * d * v)
            # F=e^{-idv}(p+qv+rv^2); F'=e^{-idv}(q+2rv -id(p+qv+rv^2))
            # match r*(-id)= P2 => r=iP2/d ; 2r - id q = P1 => q=(2r-P1)/(id); q - id p = P0 => p=(q-P0)/(id)
            r = 1j * P2 / d
            q = (2 * r - P1) / (1j * d)
            p = (q - P0) / (1j * d)
            return e * (p + q * v + r * v * v)
        tot += antiquad(rgt, d) - antiquad(lft, d)
    return tot

def gpair_exact(cj, ck, u):
    tot = 0.0 + 0.0j
    for a in (cj, -cj):
        for b in (ck, -ck):
            tot += np.exp(1j * u * b) * Jexact(u, a - b)
    return np.real(tot) / (2 * np.pi) ** 2

if __name__ == "__main__":
    cs = F.probes(F.BAND, 3.0)
    fn = F.build(cs)                     # numerical instrument (validated §376)
    prime_num = fn['prime']
    n = len(cs)
    # build closed-form prime matrix
    prime_cf = np.zeros((n, n))
    for j in range(n):
        for k in range(j, n):
            v = 2.0 * np.sum(F.LAMW * np.array([gpair_exact(cs[j], cs[k], np.log(m)) for m in F.NS]))
            prime_cf[j, k] = prime_cf[k, j] = v
    err = np.max(np.abs(prime_cf - prime_num))
    rel = err / (np.max(np.abs(prime_num)) + 1e-300)
    print(f"closed-form vs numerical Fejer prime matrix: maxabs err {err:.3e}  rel {rel:.3e}")
    print("(also cross-check the grid-J vs antideriv-J on a few entries)")
    for (a, b, u) in [(14.0, 14.0, np.log(2)), (14.0, 20.0, np.log(3)), (44.0, 44.0, np.log(2))]:
        jg = J(u, a - b); je = Jexact(u, a - b)
        print(f"  J({u:.3f},{a-b:.0f}): grid {jg:.4e}  antideriv {je:.4e}  diff {abs(jg-je):.2e}")
    print("GATE: rel err small => closed-form Fejer prime side VALID => rigorous §380 cert feasible.")
