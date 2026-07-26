#!/usr/bin/env python3
"""storm#1 recognition: reconstruct the saturating Weil function f(x) and fingerprint it (zeros, decay, oscillation)."""
import mpmath as mp
mp.mp.dps=22
def Omega(r): return mp.re(mp.digamma(mp.mpf(1)/4+1j*r/2))-mp.log(mp.pi)
w=mp.mpf(3)
def F(c,r): return mp.e**(-(r-c)**2/(2*w**2))+mp.e**(-(r+c)**2/(2*w**2))
def Fc(c,z): return mp.e**(-(z-c)**2/(2*w**2))+mp.e**(-(z+c)**2/(2*w**2))
def bpts(a,b):
    p=[mp.mpf(0)]
    for c in (a,b):
        for x in (c-4*w,c,c+4*w):
            if x>0:p.append(x)
    return sorted(set(p))+[mp.inf]
Lam={}
for pp in range(2,50):
    if all(pp%d for d in range(2,int(pp**.5)+1)):
        q=pp
        while q<50: Lam[q]=mp.log(pp);q*=pp
centers=[mp.mpf(14)+3*i for i in range(7)]  # [14,32], 7 probes > ~5 zeros
M=len(centers)
print(f"reconstruct saturating f: {M} probes on [14,32]",flush=True)
B=mp.zeros(M)
for j in range(M):
    for k in range(M):
        pts=bpts(centers[j],centers[k])
        arch=(1/(2*mp.pi))*2*mp.quad(lambda r:F(centers[j],r)*F(centers[k],r)*Omega(r),pts)
        pole=mp.re(Fc(centers[j],1j/2)*Fc(centers[k],1j/2))*2
        pr=-2*sum(Lam[n]/mp.sqrt(n)*(1/(2*mp.pi))*2*mp.quad(lambda r:F(centers[j],r)*F(centers[k],r)*mp.cos(mp.log(n)*r),pts) for n in sorted(Lam))
        B[j,k]=pole+arch+pr
ev,V=mp.eigsy(B)
a=[V[i,0] for i in range(M)]
print(f"  min-eig={mp.nstr(ev[0],4)}; saturating coeffs a = {[mp.nstr(x,4) for x in a]}",flush=True)
def f(x): return sum(a[j]*F(centers[j],x) for j in range(M))
# FINGERPRINT: zeros of f on [10,36], decay, and compare zero-spacing to zeta-zero spacing
print("  f(x) sampled on [8,38] step 1:",flush=True)
xs=[mp.mpf(8)+i for i in range(31)]
fv=[f(x) for x in xs]
print("   "+", ".join(mp.nstr(v,2) for v in fv),flush=True)
# find sign changes (zeros)
zeros_f=[]
for i in range(len(xs)-1):
    if fv[i]*fv[i+1]<0:
        z=mp.findroot(f,(xs[i]+xs[i+1])/2)
        zeros_f.append(z)
print(f"  ZEROS of f: {[mp.nstr(z,6) for z in zeros_f]}",flush=True)
zz=[mp.mpf(z) for z in ['14.134725','21.022040','25.010858','30.424876','32.935062']]
print(f"  zeta zeros in band: {[mp.nstr(z,6) for z in zz]}",flush=True)
print("  => are f's zeros AT / interlacing the zeta zeros? (interlacing = de Branges HB structure)",flush=True)
if len(zeros_f)>1:
    print(f"  f zero-spacings: {[mp.nstr(zeros_f[i+1]-zeros_f[i],4) for i in range(len(zeros_f)-1)]}",flush=True)
print("DONE",flush=True)
