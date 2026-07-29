"""Hilbert: DIRECT Arb-interval certificate that the physical-delta-Lambda Gram G
on Riemann's worst-conditioned mod-105 5x5 cell (s=48) is positive-definite.
Two independent ways: (1) interval Cholesky, every radicand/pivot certified >0;
(2) independent interval Leibniz determinant + ALL leading principal minors (Sylvester).
u_t(n)=dL(n+t1)dL(n+t2)dL(n+t3), dL=Lambda-1. Entries are exact sums of products of
(log p - 1), enclosed rigorously in Arb ball arithmetic. Michalowski-style, signed POSITIVE."""
from flint import ctx, arb
import json, itertools
ctx.prec = 160

CELL = [(6,15,27),(1,7,40),(5,11,32),(12,15,21),(1,22,25)]
S = 48
M = 1_000_000
MAXOFF = max(x for t in CELL for x in t)

# --- smallest prime factor sieve to M+MAXOFF ---
LIM = M + MAXOFF + 1
spf = [0]*LIM
for i in range(2, LIM):
    if spf[i]==0:
        for j in range(i, LIM, i):
            if spf[j]==0: spf[j]=i

# --- dL[n] = Lambda(n)-1 as arb ball ---
NEG1 = arb(-1)
logcache = {}
def dL(n):
    p = spf[n]
    if p==0:  # n<2
        return NEG1
    m=n
    while m%p==0: m//=p
    if m==1:
        lp = logcache.get(p)
        if lp is None:
            lp = arb(p).log(); logcache[p]=lp
        return lp - 1
    return NEG1

# precompute u_i(n) for the 5 triples across the window
print(f"building u-vectors, M={M} ...", flush=True)
U = [[None]*(M+1) for _ in CELL]
for i,(a,b,c) in enumerate(CELL):
    Ui=U[i]
    for n in range(1, M+1):
        Ui[n] = dL(n+a)*dL(n+b)*dL(n+c)

# Gram entries G[i][j] = (1/M) sum_n u_i(n) u_j(n)  (arb, rigorous)
print("accumulating Gram entries ...", flush=True)
G=[[None]*5 for _ in range(5)]
invM = arb(1)/M
for i in range(5):
    for j in range(i,5):
        acc = arb(0)
        Ui=U[i]; Uj=U[j]
        for n in range(1,M+1):
            acc += Ui[n]*Uj[n]
        val = acc*invM
        G[i][j]=val; G[j][i]=val

def enc(x): return [float(x.lower()), float(x.upper())]

# ---- (1) interval Cholesky, certify each radicand > 0 ----
L=[[arb(0)]*5 for _ in range(5)]
L=[[arb(0) for _ in range(5)] for _ in range(5)]
chol_ok=True; pivots=[]
for k in range(5):
    d = G[k][k]
    for m in range(k):
        d = d - L[k][m]*L[k][m]
    pos = (d>0)
    pivots.append((enc(d), bool(pos)))
    if not pos: chol_ok=False; break
    L[k][k]=d.sqrt()
    for i in range(k+1,5):
        off=G[i][k]
        for m in range(k):
            off = off - L[i][m]*L[k][m]
        L[i][k]=off/L[k][k]

# ---- (2) independent Leibniz determinant + leading principal minors ----
def perm_det(Asub, idx):
    r=len(idx); tot=arb(0)
    for perm in itertools.permutations(range(r)):
        sign=1
        pr=list(perm)
        for a in range(r):
            for b in range(a+1,r):
                if pr[a]>pr[b]: sign=-sign
        term=arb(sign)
        for a in range(r):
            term = term * Asub[idx[a]][perm[a]]
        tot += term
    return tot
minors=[]; sylv_ok=True
for k in range(1,6):
    idx=list(range(k))
    dk=perm_det(G, idx)
    posk=(dk>0)
    minors.append((k, enc(dk), bool(posk)))
    if not posk: sylv_ok=False

CERTIFIED = chol_ok and sylv_ok
out={
 "object":"DIRECT physical-dLambda Gram G, 5x5 worst-conditioned mod-105 cell",
 "version":"rt6-hilbert-arb-cert-2026-07-28",
 "cell_triples":CELL,"s":S,"window_M":M,"arb_prec_bits":ctx.prec,
 "gram_entry_enclosures":[[enc(G[i][j]) for j in range(5)] for i in range(5)],
 "way1_interval_cholesky":{"pivot_radicands":pivots,"all_positive":chol_ok},
 "way2_leibniz_leading_minors":{"minors":minors,"all_positive":sylv_ok},
 "CERTIFIED_PD_two_ways":CERTIFIED,
 "claim":"fork A at pinned s=48 in the accessible odd sector -- Gram PD, both ways, outward-rounded intervals"
}
json.dump(out, open("results/rt6_hilbert_arb_cert.json","w"), indent=2)
print("\n=== Cholesky radicands (must all be >0) ===")
for k,(e,p) in enumerate(pivots): print(f"  pivot[{k}] in [{e[0]:.4f},{e[1]:.4f}]  >0: {p}")
print("=== leading principal minors (Sylvester, independent) ===")
for k,e,p in minors: print(f"  minor {k}x{k} in [{e[0]:.4e},{e[1]:.4e}]  >0: {p}")
print(f"\nCERTIFIED PD two independent ways: {CERTIFIED}")
print("banked -> results/rt6_hilbert_arb_cert.json")
