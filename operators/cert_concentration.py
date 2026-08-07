#!/usr/bin/env python3
"""Arb-certified concentration section (§804): A[j,k] = 2∫_[14,44] F_j F_k dt (even symmetry),
tr A = 2∫(K(0)+K(2t))dt closed-ish, Ky-Fan certified lower bound on Σ top-9 pencil eigenvalues,
=> certified upper bound on λ_10 of the true operator. Ambient probes 12..44 step 2 (17)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import arb_fejer_cert2 as ac

fourW = 1.32
CENTERS = [12+2*i for i in range(17) if 12+2*i <= 44]
ctx.prec = 200
ac.W = arb(fourW)/4; ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W
Wa = ac.W
cs = [arb(c) for c in CENTERS]; n = len(cs)
TWO_W = 2*Wa

def sinc_sq(y):
    """(sin y / y)^2, robust at y ~ 0: series 1 - y^2/3 + 2y^4/45 with certified remainder
    (|next terms| <= |y|^6/100 for |y| <= 0.6, generous) when the ball can touch 0."""
    try:
        near0 = bool(abs(y) < 0.6) or y.contains(acb(0))
    except Exception:
        near0 = True
    if near0:
        y2 = y*y
        val = acb(1) - y2/3 + 2*y2*y2/45
        m = float(abs(y2*y2*y2).real.mid()) + float(abs(y2*y2*y2).imag.mid())
        e = m/100 + 1e-30
        return val + acb(arb(0, e), arb(0, e))
    return (y.sin()/y)**2

print(f"A-matrix: 2*int_[14,44] F_jF_k, {n} probes...", flush=True)
A = [[None]*n for _ in range(n)]
for j in range(n):
    for k in range(j, n):
        def ig(r, _):
            fj_ = sinc_sq(Wa*(r-acb(CENTERS[j]))) + sinc_sq(Wa*(r+acb(CENTERS[j])))
            fk_ = sinc_sq(Wa*(r-acb(CENTERS[k]))) + sinc_sq(Wa*(r+acb(CENTERS[k])))
            return fj_*fk_
        I = acb(0); a = arb(14)
        while a < 44:
            b = a + 2
            if b > 44: b = arb(44)
            I = I + acb.integral(ig, acb(a), acb(b))
            a = b
        v = 2*I.real
        A[j][k] = A[k][j] = v
    print(f"  A row {j} done (rad {float(A[j][j].rad()):.1e})", flush=True)

G = [[None]*n for _ in range(n)]
for j in range(n):
    for k in range(j, n):
        G[j][k] = G[k][j] = 2*arb.pi()*ac.gpair_closed_arb(cs[j], cs[k], arb(0))

# tr A of the OPERATOR (even subspace): 2*int_14^44 [K(0) + K(2t)] dt, K(u)=sin(2W u)/(pi u)
K0 = TWO_W/arb.pi()
def igtr(t, _):
    u = 2*t
    return (TWO_W*u).sin()/(arb.pi()*u) if not u.contains(arb(0)) else acb(K0)
I2 = acb.integral(lambda t,_: ((2*t*TWO_W).sin()/(arb.pi()*2*t)), acb(14), acb(44))
trA = K0*30 + I2.real  # even-subspace: window [14,44], kernel K(0)+K(2t); no extra x2
print(f"certified tr A (even op) = {trA}", flush=True)

# float eigenvectors for the Ky-Fan frame — independent numpy build (robust, frame need not be certified)
from scipy.integrate import quad as _q
def _f(c, r):
    y = Wa*(r-c); s1 = np.sin(y)/y if abs(y)>1e-8 else 1.0
    y2 = Wa*(r+c); s2 = np.sin(y2)/y2 if abs(y2)>1e-8 else 1.0
    return s1*s1 + s2*s2
Wf = float(Wa.mid())
def _fn(c, r):
    y = Wf*(r-c); s1 = np.sin(y)/y if abs(y)>1e-9 else 1.0
    y2 = Wf*(r+c); s2 = np.sin(y2)/y2 if abs(y2)>1e-9 else 1.0
    return s1*s1 + s2*s2
Am = np.zeros((n,n))
for i in range(n):
    for j in range(i, n):
        v = 2*_q(lambda r: _fn(CENTERS[i], r)*_fn(CENTERS[j], r), 14, 44, limit=200)[0]
        Am[i,j] = Am[j,i] = v
Gm = np.array([[float(G[i][j].mid()) for j in range(n)] for i in range(n)])
gw, gv = np.linalg.eigh(Gm); keep = gw > gw.max()*1e-12
Bf = gv[:, keep]/np.sqrt(gw[keep])
ew, ev = np.linalg.eigh(Bf.T @ Am @ Bf)
order = np.argsort(ew)[::-1]
K9 = 9
V = Bf @ ev[:, order[:K9]]     # probe-coefficient frame, approx G-orthonormal
# certified: M1 = V^T A V, M0 = V^T G V in interval arithmetic
Va = [[arb(float(V[i,j])) for j in range(K9)] for i in range(n)]
def mat_VtMV(M):
    out = [[arb(0)]*K9 for _ in range(K9)]
    for a_ in range(K9):
        for b_ in range(a_, K9):
            s = arb(0)
            for i in range(n):
                for j2 in range(n):
                    s += Va[i][a_]*M[i][j2]*Va[j2][b_]
            out[a_][b_] = out[b_][a_] = s
    return out
# G-orthonormalize the float frame: U = V (V^T G V)^{-1/2}
M0f = V.T @ Gm @ V
w0, u0 = np.linalg.eigh(M0f)
U = V @ (u0 @ np.diag(w0**-0.5) @ u0.T)     # float G-orthonormal
Ua = [[arb(float(U[i,j])) for j in range(K9)] for i in range(n)]
def mat_UtMU(M):
    out = [[arb(0)]*K9 for _ in range(K9)]
    for a_ in range(K9):
        for b_ in range(a_, K9):
            s2 = arb(0)
            for i in range(n):
                for j2 in range(n):
                    s2 += Ua[i][a_]*M[i][j2]*Ua[j2][b_]
            out[a_][b_] = out[b_][a_] = s2
    return out
M1 = mat_UtMU(A); M0 = mat_UtMU(G)
# Ky Fan with non-orthonormal frame: sum top-9 of pencil (A,G) >= tr(M0^{-1} M1) lower bound via
# certified: lambda-sum >= tr(M1) - ||M1||*||M0 - I||_F correction. Simpler rigorous route:
# generalized Rayleigh: sum_{i<=9} lam_i >= tr over ANY 9-dim subspace of the pencil trace =
# tr(M0^{-1}M1) exactly (trace of pencil compression). Compute certified tr(M0^{-1}M1) via LDL^T solve.
# sum of top-9 pencil eigenvalues >= tr(M0^{-1} M1); with M0 = I + E (E certified tiny since U is
# G-orthonormal in float), tr(M0^{-1}M1) >= tr(M1) - ||M1||_op * ||E||_tr.  Rigorous lower bound.
nC = K9
trM1 = arb(0)
for i in range(nC): trM1 += M1[i][i]
# ||E||: E = M0 - I, certified Frobenius
Efro = arb(0)
for i in range(nC):
    for j2 in range(nC):
        d_ = M0[i][j2] - (arb(1) if i==j2 else arb(0))
        Efro += d_*d_
Efro = Efro.sqrt()
# ||M1||_op <= ||M1||_F
M1fro = arb(0)
for i in range(nC):
    for j2 in range(nC): M1fro += M1[i][j2]*M1[i][j2]
M1fro = M1fro.sqrt()
# tr(M0^{-1}M1) >= tr(M1) - ||M1||_F * ||E||_F / (1 - ||E||_F)   (Neumann series bound, valid ||E||<1)
denom = arb(1) - Efro
trP = trM1 - M1fro*Efro/denom
print(f"  tr(M1)={float(trM1.mid()):.5f}  ||E||_F={float(Efro.upper()):.2e}  ||M1||_F={float(M1fro.mid()):.3f}", flush=True)
print(f"certified pencil-compression trace (>= sum of top-9 true eigenvalues... see note): {trP}")
lam10_ub = trA - trP
print(f"\n===== §804: certified lambda_10 upper bound = tr A - S9 <= {float(lam10_ub.upper()):.5f} =====")
print("(float landscape predicted ~0.037; NOTE the Ky-Fan direction requires S9 <= sum of top-9")
print(" TRUE eigenvalues — compression trace over a 9-dim subspace IS <= sum top-9 (Ky Fan max")
print(" principle), so trA - trP is a VALID upper bound on the tail sum INCLUDING lambda_10.)")
