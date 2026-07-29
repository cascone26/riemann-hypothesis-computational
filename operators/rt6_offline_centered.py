#!/usr/bin/env python3
"""§431b: fair intrinsic test — center the band on each moved zero (local density), so β* measures the INTRINSIC
sensitivity to off-line displacement at height γ, not in-band/out-of-band effects. β* vs γ: growing => high-γ
off-line genuinely weaker signal (confirms §422 intrinsically); flat => height-independent (the §431 growth was
just fixed-band geometry)."""
import numpy as np, l2_fejer as LF
Z=np.load('zeros150.npy'); W=LF.W
LF.RGRID=np.arange(-360.0,360.0,0.006); from scipy.special import digamma
LF.OMEGA=np.real(digamma(0.25+0.5j*LF.RGRID))-np.log(np.pi)
def fej_c(x): y=W*x; return (np.sin(y)/y)**2
def Frow_c(c,r): return fej_c(r-c)+fej_c(r+c)
def zeroside(cs,kmove,beta,zin):
    n=len(cs); M=np.zeros((n,n))
    for i in zin:
        if i==kmove: continue
        Fv=np.array([Frow_c(c,Z[i]).real for c in cs]); M+=2.0*np.outer(Fv,Fv)
    g=Z[kmove]
    for arg in (g-1j*beta,g+1j*beta):
        Fv=np.array([Frow_c(c,arg) for c in cs]); M+=np.real(np.outer(Fv,Fv))
    return 0.5*(M+M.T)
print(f"§431b intrinsic off-line sensitivity (band centered on each zero). β* vs γ.")
print(f"{'zero γ':>9} {'k':>3} {'half-wid':>8} {'n':>3} {'local Δγ':>8} {'β*':>9}")
for k in [3,7,12,20,30,50,80,110,140]:
    if k>=len(Z): continue
    g=Z[k]
    hw=15.0
    band=(g-hw,g+hw)
    zin=[i for i in range(len(Z)) if band[0]<=Z[i]<=band[1]]
    # local spacing
    dloc=(Z[min(k+1,len(Z)-1)]-Z[max(k-1,0)])/2
    cs=LF.probes(band, dloc)   # probe at local zero spacing
    if len(cs)<4: continue
    l0=float(np.linalg.eigvalsh(zeroside(cs,k,0.0,zin))[0])
    lo,hi=0.0,3.0
    if float(np.linalg.eigvalsh(zeroside(cs,k,hi,zin))[0])>=0: bstar=np.inf
    else:
        for _ in range(28):
            mid=(lo+hi)/2
            if float(np.linalg.eigvalsh(zeroside(cs,k,mid,zin))[0])>=0: lo=mid
            else: hi=mid
        bstar=(lo+hi)/2
    print(f"{g:9.3f} {k:3d} {hw:8.1f} {len(cs):3d} {dloc:8.3f} {bstar:9.4f}",flush=True)
print("""
reading: β* vs γ with band centered on the zero. GROWING β* => high-γ off-line displacement intrinsically harder
to detect (weaker signal per unit β) — confirms §422's wall as an intrinsic feature. FLAT β* => the §431 growth
was fixed-band geometry, and intrinsic off-line sensitivity is height-independent (then §422's wall is only about
the EFFECTIVE/global argument, not local detectability).
""")
