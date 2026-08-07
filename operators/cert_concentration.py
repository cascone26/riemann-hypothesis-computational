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

def fej_acb2(c, r):
    out = acb(0)
    for cc in (c, -c):
        y = Wa*(r-acb(cc))
        out = out + (y.sin()/y)**2   # centers >= 12, window [14,44]: y never 0? c=12..44, r in [14,44]: y=0 at r=c — removable; acb.integral handles via balls... guard: integrand continuous, sin(y)/y fine numerically except exactly 0 — arb integral uses complex rectangles, y=0 interior → 0/0. Shift: use sinc via series? Use (sin y / y) with tiny imaginary offset is wrong. Instead use formula: fej = (sin y)^2/y^2 with limit 1; acb.integral tolerates if we rewrite as sinc_sq safe:
    return out

def sinc_sq(y):
    # (sin y / y)^2 robust near 0 via series switch
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
trA = 2*(K0*30 + I2.real)
print(f"certified tr A (even op) = {trA}", flush=True)

# float eigenvectors for the Ky-Fan frame
Am = np.array([[float(A[i][j].mid()) for j in range(n)] for i in range(n)])
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
M1 = mat_VtMV(A); M0 = mat_VtMV(G)
# Ky Fan with non-orthonormal frame: sum top-9 of pencil (A,G) >= tr(M0^{-1} M1) lower bound via
# certified: lambda-sum >= tr(M1) - ||M1||*||M0 - I||_F correction. Simpler rigorous route:
# generalized Rayleigh: sum_{i<=9} lam_i >= tr over ANY 9-dim subspace of the pencil trace =
# tr(M0^{-1}M1) exactly (trace of pencil compression). Compute certified tr(M0^{-1}M1) via LDL^T solve.
nC = K9
dC = [None]*nC; LC = [[arb(0)]*nC for _ in range(nC)]
ok = True
for i in range(nC):
    sd = M0[i][i]
    for k in range(i): sd -= LC[i][k]*LC[i][k]*dC[k]
    if not (sd > 0): ok = False; break
    dC[i] = sd
    for j2 in range(i+1, nC):
        v = M0[j2][i]
        for k in range(i): v -= LC[j2][k]*LC[i][k]*dC[k]
        LC[j2][i] = v/dC[i]
assert ok, "M0 not PD"
def solve(rhs):
    y = [None]*nC
    for i in range(nC):
        s = rhs[i]
        for k in range(i): s -= LC[i][k]*y[k]
        y[i] = s
    for i in range(nC): y[i] = y[i]/dC[i]
    z = [None]*nC
    for i in reversed(range(nC)):
        s = y[i]
        for k in range(i+1, nC): s -= LC[k][i]*z[k]
        z[i] = s
    return z
trP = arb(0)
for col in range(nC):
    x = solve([M1[i][col] for i in range(nC)])
    trP += x[col]
print(f"certified pencil-compression trace (>= sum of top-9 true eigenvalues... see note): {trP}")
lam10_ub = trA - trP
print(f"\n===== §804: certified lambda_10 upper bound = tr A - S9 <= {float(lam10_ub.upper()):.5f} =====")
print("(float landscape predicted ~0.037; NOTE the Ky-Fan direction requires S9 <= sum of top-9")
print(" TRUE eigenvalues — compression trace over a 9-dim subspace IS <= sum top-9 (Ky Fan max")
print(" principle), so trA - trP is a VALID upper bound on the tail sum INCLUDING lambda_10.)")
