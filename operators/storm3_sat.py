#!/usr/bin/env python3
"""storm#3/#4: is the near-saturating (min-eig) Weil test function's f(gamma) ~0 at in-band zeros? (§290 claim)."""
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
ZEROS=[mp.mpf(z) for z in ['14.134725','21.022040','25.010858','30.424876','32.935062','37.586178','40.918719','43.327073','48.005151','49.773832']]
centers=[mp.mpf(14)+3*i for i in range(13)]  # [14,50], 13 probes > 10 zeros -> saturating direction
M=len(centers)
print(f"storm3: build {M}x{M} Weil form, band [14,50], {M} probes vs {len(ZEROS)} in-band zeros",flush=True)
B=mp.zeros(M)
for j in range(M):
    for k in range(M):
        pts=bpts(centers[j],centers[k])
        arch=(1/(2*mp.pi))*2*mp.quad(lambda r:F(centers[j],r)*F(centers[k],r)*Omega(r),pts)
        pole=mp.re(Fc(centers[j],1j/2)*Fc(centers[k],1j/2))*2
        pr=-2*sum(Lam[n]/mp.sqrt(n)*(1/(2*mp.pi))*2*mp.quad(lambda r:F(centers[j],r)*F(centers[k],r)*mp.cos(mp.log(n)*r),pts) for n in sorted(Lam))
        B[j,k]=pole+arch+pr
    print(f"  row {j} done",flush=True)
ev,V=mp.eigsy(B)
print(f"  eigenvalues: min={mp.nstr(ev[0],4)} max={mp.nstr(ev[-1],4)}  (min~0 => saturating direction exists)",flush=True)
a=[V[i,0] for i in range(M)]
norm=mp.sqrt(sum(abs(x)**2 for x in a))
vals=[abs(sum(a[j]*F(centers[j],g) for j in range(M)))/norm for g in ZEROS]
print("  |f(gamma)|/||a|| at in-band zeros:",flush=True)
print("   "+", ".join(mp.nstr(v,3) for v in vals),flush=True)
print(f"  MAX = {mp.nstr(max(vals),5)}  => {'CONFIRMED §290 (saturating f vanishes at zeros)' if max(vals)<0.05 else 'NOT vanishing'}",flush=True)
# also: does f(gamma) vanish at OUT-of-band zeros too? (it should NOT need to)
OUT=[mp.mpf('52.970321'),mp.mpf('59.347044')]
ovals=[abs(sum(a[j]*F(centers[j],g) for j in range(M)))/norm for g in OUT]
print(f"  |f| at out-of-band zeros (52.97,59.35): {', '.join(mp.nstr(v,3) for v in ovals)} (expect NOT ~0: f localizes to its band)",flush=True)
print("DONE",flush=True)
