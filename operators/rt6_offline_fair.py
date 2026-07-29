#!/usr/bin/env python3
"""§431c (fixed): intrinsic off-line sensitivity. ALL 150 zeros in the sum (complete, PSD), band CENTERED on the
moved zero, probes WELL BELOW edge (spacing=2×local => ~half zero count => well-conditioned λmin(β=0)>0). β* vs γ
= intrinsic detectability of off-line displacement at height γ. This is the clean §422 intrinsic test."""
import numpy as np, l2_fejer as LF
Z=np.load('zeros150.npy'); W=LF.W
def fej_c(x): y=W*x; return (np.sin(y)/y)**2
def Frow_c(c,r): return fej_c(r-c)+fej_c(r+c)
def zeroside(cs,kmove,beta):
    n=len(cs); M=np.zeros((n,n))
    Fon=np.array([[Frow_c(c,Z[i]).real for c in cs] for i in range(len(Z)) if i!=kmove])
    M=2.0*Fon.T@Fon
    g=Z[kmove]
    for arg in (g-1j*beta,g+1j*beta):
        Fv=np.array([Frow_c(c,arg) for c in cs]); M=M+np.real(np.outer(Fv,Fv))
    return 0.5*(M+M.T)
print("§431c intrinsic off-line sensitivity, band centered, all zeros, sub-edge probes. β* vs γ.")
print(f"{'zero γ':>9} {'k':>3} {'n':>3} {'localΔγ':>8} {'λmin(β=0)':>11} {'β*':>9}")
res=[]
for k in [3,7,12,20,30,45,60,80,100,120,140]:
    if k>=len(Z): continue
    g=Z[k]; dloc=(Z[min(k+1,149)]-Z[max(k-1,0)])/2
    band=(g-18,g+18); cs=LF.probes(band, 2.0*dloc)   # ~half density => sub-edge
    if len(cs)<4: continue
    l0=float(np.linalg.eigvalsh(zeroside(cs,k,0.0))[0])
    if l0<1e-9:   # not well-conditioned, skip
        print(f"{g:9.3f} {k:3d} {len(cs):3d} {dloc:8.3f} {l0:+11.3e}   (ill-cond, skip)"); continue
    lo,hi=0.0,4.0
    if float(np.linalg.eigvalsh(zeroside(cs,k,hi))[0])>=0: bstar=np.inf
    else:
        for _ in range(30):
            mid=(lo+hi)/2
            if float(np.linalg.eigvalsh(zeroside(cs,k,mid))[0])>=0: lo=mid
            else: hi=mid
        bstar=(lo+hi)/2
    res.append((g,bstar))
    print(f"{g:9.3f} {k:3d} {len(cs):3d} {dloc:8.3f} {l0:+11.3e} {bstar:9.4f}",flush=True)
if len(res)>=3:
    gs=np.array([r[0] for r in res]); bs=np.array([r[1] for r in res])
    fin=np.isfinite(bs)
    if fin.sum()>=3:
        slope=np.polyfit(np.log(gs[fin]),bs[fin],1)[0]
        print(f"\nβ* vs log γ slope = {slope:+.3f}  ({'GROWS with γ (high-γ harder, confirms §422 intrinsically)' if slope>0.05 else 'FLAT (height-independent; §422 wall is effective/global not local)' if abs(slope)<=0.05 else 'shrinks'})")
