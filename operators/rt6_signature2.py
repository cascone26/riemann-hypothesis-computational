#!/usr/bin/env python3
"""§399 (§392 follow-up a+b): characterize arch-prime's ONE negative direction.
(a) Is its eigenvector the §375 zero-interpolation mode (phi vanishes on zeros)?
(b) Which prime terms build the negative direction? (the one-sided sieve target)
(c) Robustness: is it ALWAYS exactly 1 negative across bands/widths/fills, and what's its magnitude law?"""
import numpy as np
import l2_edge_sweep as L
w = 3.0
def analyze(cs, band_desc):
    f = L.build_forms(cs, w)
    z = 0.5j
    fp = np.exp(-((z-cs)**2)/(2*w*w)) + np.exp(-((z+cs)**2)/(2*w*w))
    pole = 2.0*np.real(np.outer(fp, fp))
    K = f['full'] - pole                     # arch - prime
    ev, V = np.linalg.eigh(K)
    nneg = int((ev < -1e-11*max(1,abs(ev).max())).sum())
    if nneg == 0:
        return f"{band_desc}: n={len(cs)} nneg=0 (K already PSD)"
    vneg = V[:, 0]                            # most-negative eigenvector
    # (a) does phi=sum vneg F vanish on zeros?
    zin = L.ZEROS[(L.ZEROS>=cs.min()-2*w)&(L.ZEROS<=cs.max()+2*w)]
    def phi(vv, r): return sum(vv[j]*(np.exp(-((r-c)**2)/(2*w*w))+np.exp(-((r+c)**2)/(2*w*w))) for j,c in enumerate(cs))
    rg = np.linspace(cs.min()-4, cs.max()+4, 1200); sc = max(abs(phi(vneg,x)) for x in rg)+1e-30
    zvals = np.mean([abs(phi(vneg,g))/sc for g in zin]) if len(zin) else float('nan')
    mids = (zin[:-1]+zin[1:])/2 if len(zin)>1 else zin
    mvals = np.mean([abs(phi(vneg,m))/sc for m in mids]) if len(mids) else float('nan')
    # (b) prime content of the negative direction: <vneg, P_n vneg>
    def pform(n):
        u=np.log(n); M=np.zeros((len(cs),)*2)
        for i in range(len(cs)):
            for j in range(i,len(cs)): M[i,j]=M[j,i]=2.0*L.g_pair(cs[i],cs[j],w,u)
        return M
    pc={}
    for i,n in enumerate(L.NS):
        if L.LAMW[i]==0 or L.LAMW[i]*np.exp(-(np.log(n)**2)*w*w/4)<1e-13: continue
        pc[int(n)]=L.LAMW[i]*float(vneg@pform(int(n))@vneg)
    top=sorted(pc,key=lambda k:-abs(pc[k]))[:4]
    # <vneg,arch vneg> and <vneg,prime vneg>
    arch=f['archpole']-pole
    va=float(vneg@arch@vneg); vp_=float(vneg@f['prime']@vneg)
    return (f"{band_desc}: n={len(cs)} nneg={nneg} negeig={ev[0]:+.3e} | phi(zeros)/phi(mid)={zvals/(mvals+1e-12):.3f} "
            f"| <v,arch v>={va:+.3f} <v,prime v>={vp_:+.3f} | topP {[(n,round(pc[n],3)) for n in top]}")
print("arch-prime NEGATIVE-DIRECTION characterization:")
for band,B in [("[0.5,30]",30),("[0.5,40]",40)]:
    for n in (8,12,15):
        print(analyze(np.linspace(0.5,B,n), band))
print("\nreading: phi(zeros)/phi(mid)<<1 => neg direction IS the §375 zero-interpolation mode (arch-prime")
print("goes negative exactly along the zero-interpolant). <v,prime v> > <v,arch v> => the negativity is")
print("prime-DRIVEN (prime overwhelms arch on that one mode) => one-sided target = bound that ONE prime form.")
