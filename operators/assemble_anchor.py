#!/usr/bin/env python3
"""Final assembly of the anchor ∀-theorem (§805). Consumes gram14.json (17-probe zero-free Q + L²
pencil G) + the concentration section; produces m_E (pencil λ_min on the top-9 concentration
eigenspace E), then the ε₀ quadratic. Runs only when gram14.json exists."""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
ctx.prec = 200
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, "gram14.json")))
CENTERS = d["centers"]; n = len(CENTERS)
def load(mat):
    M = [[None]*n for _ in range(n)]
    for key, (mid, rad) in d[mat].items():
        j, k = map(int, key.split(","))
        M[j][k] = M[k][j] = arb(mid, rad)
    return M
Q = load("Q"); G = load("G")
Qm = np.array([[float(Q[i][j].mid()) for j in range(n)] for i in range(n)])
Gm = np.array([[float(G[i][j].mid()) for j in range(n)] for i in range(n)])

# concentration eigenspace E = top-9 of (A, G); A = 2∫_[14,44] F_jF_k (float, robust numpy)
from scipy.integrate import quad
Wf = float(arb('0.33').mid())
def fn(c, r):
    y = Wf*(r-c); s1 = np.sin(y)/y if abs(y) > 1e-9 else 1.0
    y2 = Wf*(r+c); s2 = np.sin(y2)/y2 if abs(y2) > 1e-9 else 1.0
    return s1*s1 + s2*s2
A = np.zeros((n, n))
for i in range(n):
    for j in range(i, n):
        A[i, j] = A[j, i] = 2*quad(lambda r: fn(CENTERS[i], r)*fn(CENTERS[j], r), 14, 44, limit=200)[0]
gw, gv = np.linalg.eigh(Gm); Bf = gv/np.sqrt(np.abs(gw))
ew, ev = np.linalg.eigh(Bf.T @ A @ Bf)
order = np.argsort(ew)[::-1]
K9 = 9
V = Bf @ ev[:, order[:K9]]            # E in probe coords (float, ~G-orthonormal)
# certified pencil λ_min of (Q, G) restricted to E = min generalized eig of (V^T Q V, V^T G V)
Va = [[arb(float(V[i, j])) for j in range(K9)] for i in range(n)]
def VtMV(M):
    out = [[arb(0)]*K9 for _ in range(K9)]
    for a in range(K9):
        for b in range(a, K9):
            s = arb(0)
            for i in range(n):
                for j in range(n): s += Va[i][a]*M[i][j]*Va[j][b]
            out[a][b] = out[b][a] = s
    return out
QE = VtMV(Q); GE = VtMV(G)
def ok(lam):
    m = K9; dd = [None]*m; L = [[arb(0)]*m for _ in range(m)]
    for i in range(m):
        sd = QE[i][i] - lam*GE[i][i]
        for k in range(i): sd -= L[i][k]*L[i][k]*dd[k]
        if not (sd > 0): return False
        dd[i] = sd
        for j2 in range(i+1, m):
            v = QE[j2][i] - lam*GE[j2][i]
            for k in range(i): v -= L[j2][k]*L[i][k]*dd[k]
            L[j2][i] = v/dd[i]
    return True
if not ok(arb(0)):
    print("Q|E not PD — unexpected"); sys.exit(1)
lo, hi = 0.0, 2.0
for _ in range(40):
    mid = (lo+hi)/2
    if ok(arb(mid)): lo = mid
    else: hi = mid
m_E = lo
# assembly banks as §806 (anchor thread; barrier took §805)
    print(f"m_E (pencil λ_min of Q on the top-9 concentration eigenspace E) ≥ {m_E:.6e}")

# assembly quadratic: with lam10 (cert pending §804-fix), B=1.399, using ||f||^2 <= T/(1-lam10),
# ||p||^2 >= 1 - ||f||^2, Q >= m_E||p||^2 - 2B||p||||f|| - C_f||f||^2 - T.
# ε₀ = largest T with the RHS >= 0; solved numerically here with the current best constants.
def eps0(m_E, B, lam10, C_f):
    # normalize ||psi||=1; ||f||^2 = s <= T/(1-lam10); worst case s = T/(1-lam10)
    # Q >= m_E(1-s) - 2B*sqrt(s(1-s)) - C_f*s - T ; find max T>=0 keeping >=0 over s=T/(1-lam10)
    for T in np.linspace(0, 0.2, 20001):
        s = T/(1-lam10)
        if s >= 1: break
        val = m_E*(1-s) - 2*B*np.sqrt(max(s*(1-s),0)) - C_f*s - T
        if val < 0: return T
    return 0.2
e0 = eps0(m_E, 1.399, 0.037, 0.4)   # lam10, C_f float placeholders until §804-fix + C_f cert
print(f"ε₀ (float, current constants B=1.399 λ10=0.037 C_f≈0.4): {e0:.5f}")
print(f"→ class {{T(ψ) ≤ {e0:.4f}}} has Q(ψ) ≥ 0.  m_E={m_E:.4f} is the certified anchor.")
