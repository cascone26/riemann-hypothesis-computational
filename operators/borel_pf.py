#!/usr/bin/env python3
"""§323b: CORRECTED test — B∈LP ⟺ {c_k=γ_k/k!} is a PÓLYA FREQUENCY sequence (Aissen-Schoenberg-Whitney-Edrei):
generating function Σc_k y^k entire with only negative real zeros ⟺ the (lower-triangular banded) Toeplitz
matrix T=[c_{i-j}] is TOTALLY POSITIVE (all minors ≥0). §323's Hankel test was the wrong object (that tests a
moment property, not zero-negativity) — honest self-catch.
Unconditional structural angle: c_k=γ_k/k!=M_k/((2k)! k!), M_k a Stieltjes moment (Φ≥0). Test whether the
1/((2k)!k!) Bessel weighting yields a PF sequence.
Tests: (a) 2x2 Toeplitz minors = log-concavity c_k²≥c_{k-1}c_{k+1} (already ✓ §323); (b) 3x3 and higher Toeplitz
minors ≥0 (PF_r for r=3,4); (c) direct: is Σc_k y^k's partial sums' zeros all negative-real?"""
import mpmath as mp, pickle
mp.mp.dps=60
gam=[mp.mpf(g) for g in pickle.load(open("operators/_gam130.pkl","rb"))]
K=30
c=[gam[k]/mp.factorial(k) for k in range(K)]
def toeplitz_minor(rows,cols):
    # T[i,j]=c_{i-j} (0 if i-j<0 or >K); minor over given row set, col set (equal length)
    n=len(rows)
    Mx=mp.matrix(n,n)
    for a in range(n):
        for b in range(n):
            d=rows[a]-cols[b]
            Mx[a,b]=c[d] if 0<=d<K else mp.mpf(0)
    return mp.det(Mx)
import itertools
print("PF test: minimum r×r Toeplitz minor of {c_k} (PF ⟺ all ≥0):",flush=True)
for r in (2,3,4):
    worst=None
    # sample contiguous+spread row/col index sets in [0,20]
    idxpool=list(range(0,20))
    cnt=0
    for rows in itertools.combinations(idxpool,r):
        for cols in itertools.combinations(idxpool,r):
            if sum(rows)<sum(cols): continue   # need rows≥cols elementwise-ish for nonzero; prune
            m=toeplitz_minor(list(rows),list(cols))
            if worst is None or m<worst[0]: worst=(m,rows,cols)
            cnt+=1
            if cnt>4000: break
        if cnt>4000: break
    print(f"  r={r}: min minor sampled = {mp.nstr(worst[0],4)}  {'≥0 (PF_r OK)' if worst[0]>-mp.mpf('1e-45') else 'NEGATIVE'}  at rows{worst[1]} cols{worst[2]}",flush=True)
print("",flush=True)
print("(c) partial-sum zero-reality: zeros of Σ_{k≤K} c_k y^k (expect all negative-real if PF):",flush=True)
co=[c[k] for k in range(24)]  # low->high
rts=mp.polyroots(list(reversed(co)),maxsteps=400,extraprec=250)
allneg_real=all(abs(mp.im(r))/abs(r)<mp.mpf('1e-12') and mp.re(r)<0 for r in rts)
maximag=max(abs(mp.im(r))/abs(r) for r in rts)
posreal=[r for r in rts if mp.re(r)>0 and abs(mp.im(r))/abs(r)<1e-12]
print(f"   degree-23 partial sum: max|Im|/|z|={mp.nstr(maximag,3)}  all zeros negative-real? {allneg_real}  (#pos-real={len(posreal)})",flush=True)
print("READ: all Toeplitz minors ≥0 across r ⇒ {γ_k/k!} PF ⇒ B∈LP UNCONDITIONALLY (Φ≥0 only, no RH).",flush=True)
print("DONE",flush=True)
