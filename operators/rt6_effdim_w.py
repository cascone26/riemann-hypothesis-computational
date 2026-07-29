#!/usr/bin/env python3
"""§418b: does effdim of the crux survive w-variation, or is it the §412 bandwidth artifact? Sweep w."""
import numpy as np, l2_edge_sweep as L
BAND=(14.0,44.0); R=L.RGRID; dr=R[1]-R[0]
print(f"{'w':>4} {'s':>4} {'np':>3} {'push/archm':>10} {'effdim90%':>9} {'L1/|push|':>9}  top-primes")
for w in (4.0,3.0,2.0,1.5,1.0):
    # pick s near the edge for this w (probe count ~ zero count); use s that gives ~8-9 probes
    for s in (3.8,):
        cs=L.probes(BAND,s); f=L.build_forms(cs,w)
        AP=f['archpole']; va,vv=np.linalg.eigh(AP); vstar=vv[:,0]
        P=L.Fmat(cs,w,R); h2=(vstar@P)**2
        contrib=[(n,lw*np.sum(h2*np.cos(R*np.log(n)))*dr/np.pi) for n,lw in zip(L.NS,L.LAMW) if lw>0]
        push=sum(c for _,c in contrib); archm=float(va[0])
        order=sorted(contrib,key=lambda t:-abs(t[1])); tot=sum(abs(c) for _,c in contrib)
        cum=0; ed=0
        for _,c in order:
            cum+=abs(c); ed+=1
            if cum>=0.9*tot: break
        tp=[int(n) for n,_ in order[:6]]
        print(f"{w:4.1f} {s:4.1f} {len(cs):3d} {push/archm:10.3f} {ed:9d} {tot/abs(push) if abs(push)>1e-12 else np.nan:9.1f}  {tp}",flush=True)
print("\nreading: effdim GROWS as w sharpens => finite-dim was a bandwidth artifact (§412 redux). effdim STABLE")
print("across w => genuinely finite-dim crux. L1/|push| large & growing => cancellation is the essence.")
