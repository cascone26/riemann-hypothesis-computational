#!/usr/bin/env python3
"""§418 (ATTACK): effective dimension of the crux. residual push along critical v* = (1/π)Σ_n Λ(n)/√n ĥ²(log n),
h=Σ v*_j F_j, ĥ²(u)=∫h(r)²cos(ur)dr. Decompose per prime-power n; rank |contribution|; cumulative. If a FEW
low primes carry ~all of it, the crux is finite-dim (a bounded arithmetic condition) — the geometric angle
bound (§417) reduces to finitely many prime-frequency alignments. If spread across all n, it's genuinely
infinite-dim (the full Guinand-Weil wall). Also: sign pattern of ĥ²(log n) — which primes erode vs help."""
import numpy as np, l2_edge_sweep as L
BAND=(14.0,44.0); w=3.0; R=L.RGRID; dr=R[1]-R[0]
for s in (4.2,3.9,3.6):
    cs=L.probes(BAND,s); f=L.build_forms(cs,w)
    AP=f['archpole']; va,vv=np.linalg.eigh(AP); vstar=vv[:,0]
    P=L.Fmat(cs,w,R); h=vstar@P; h2=h*h                       # h(r) on grid, h^2
    contrib=[]
    for n,lw in zip(L.NS,L.LAMW):
        if lw==0: continue
        hhat2=np.sum(h2*np.cos(R*np.log(n)))*dr               # ĥ²(log n)
        contrib.append((n,lw*hhat2/np.pi))
    push=sum(c for _,c in contrib)
    archm=float(va[0])
    order=sorted(contrib,key=lambda t:-abs(t[1]))
    tot=sum(abs(c) for _,c in contrib)
    cum=0; effdim=0
    for _,c in order:
        cum+=abs(c)
        effdim+=1
        if cum>=0.9*tot: break
    top=[(n,round(c,4)) for n,c in order[:8]]
    print(f"s={s} n_probes={len(cs)}: push={push:+.4e} archm={archm:+.4e} push/archm={push/archm:.3f}")
    print(f"   L1 mass tot={tot:.4e}  effdim(90% mass)={effdim} primes  top contribs (n,contrib)={top}")
print("""
reading: effdim small (a handful) => the crux is finite-dim: RH-on-this-band = finitely many prime-frequency
  alignments of the critical direction. Big net |push| from a few primes with CANCELLING signs (L1 >> |push|)
  => the razor-tightness is exact cancellation among a few low-prime contributions (matches §417 0.3% residual).
  effdim large => genuinely infinite-dim Guinand-Weil, no finite reduction.
""")
