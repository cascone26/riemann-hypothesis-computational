#!/usr/bin/env python3
"""§331: test the NON-TAUTOLOGICAL content of the C-C-arch = primitive-height identification (§330).
The kernel match (arch Weil kernel Ω = Gamma_R log-derivative = arch height kernel) is tautological — both share
the archimedean L-factor. The REAL claim: the arch Weil form is INDEFINITE with its NEGATIVE part = a low-dimensional
"ample/fiber" piece (the Omega<0 region, r<6.29), and removing it (C-C's Sonin compression = Hodge-index's primitive
restriction) leaves the POSITIVE "primitive" part. Test: build arch form on Gaussians spanning below+above the 6.29
crossover; compute signature (#+,#-); check negative eigenvectors localize in the low-r region; check removing them
gives positive-definite. Charisma: HANDLE (signature is a computed fact)."""
import mpmath as mp
mp.mp.dps=25
def Omega(r): return mp.re(mp.digamma(mp.mpf(1)/4+1j*r/2))-mp.log(mp.pi)
w=mp.mpf(2)
def F(c,r): return mp.e**(-(r-c)**2/(2*w**2))+mp.e**(-(r+c)**2/(2*w**2))
def bpts(a,b):
    p=[mp.mpf(0)]
    for c in (a,b):
        for x in (c-4*w,c,c+4*w):
            if x>0:p.append(x)
    return sorted(set(p))+[mp.inf]
centers=[mp.mpf(2*i) for i in range(0,16)]  # r=0,2,...,30 -> spans Omega<0 (r<6.29) and Omega>0
M=len(centers)
CROSS=mp.mpf('6.2898')
print(f"arch Weil form on {M} Gaussians r=0..30 (crossover Omega=0 at r={mp.nstr(CROSS,5)}):",flush=True)
n_lowband=sum(1 for c in centers if c<CROSS)
print(f"  centers in Omega<0 band [0,{mp.nstr(CROSS,4)}): {n_lowband}  (predicted # negative directions ~ this)",flush=True)
A=mp.zeros(M)
for j in range(M):
    for k in range(M):
        pts=bpts(centers[j],centers[k])
        A[j,k]=(1/(2*mp.pi))*2*mp.quad(lambda r:F(centers[j],r)*F(centers[k],r)*Omega(r),pts)
ev,V=mp.eigsy((A+A.T)/2)
evl=[ev[i] for i in range(M)]
npos=sum(1 for e in evl if e>mp.mpf('1e-12'))
nneg=sum(1 for e in evl if e<-mp.mpf('1e-12'))
print(f"  eigenvalues: {[mp.nstr(e,3) for e in evl]}",flush=True)
print(f"  SIGNATURE: {npos} positive, {nneg} negative  (Hodge-index-shaped iff nneg is small/controlled)",flush=True)
# localize each NEGATIVE eigenvector: its 'center of mass' in r (weighted by |coeff|^2 * center)
print(f"  negative eigenvectors' frequency localization (center-of-mass in r):",flush=True)
for i in range(M):
    if evl[i]<-mp.mpf('1e-12'):
        vv=[V[t,i] for t in range(M)]
        nrm=sum(abs(x)**2 for x in vv)
        com=sum(abs(vv[t])**2*centers[t] for t in range(M))/nrm
        print(f"    eig={mp.nstr(evl[i],3)}: freq center-of-mass r={mp.nstr(com,4)}  {'LOW-BAND (fiber)' if com<CROSS+2 else 'HIGH?!'}",flush=True)
# Now: remove the negative subspace, check the complement (primitive) is positive-definite
# project A onto the span of the POSITIVE eigenvectors and confirm positivity (trivially true) -- instead do the
# meaningful test: is the # of negative directions STABLE / matching the low-band count, and does a high-band-only
# form come out positive-definite (pure primitive)?
print("",flush=True)
hb=[c for c in centers if c>CROSS]
Mh=len(hb)
Ah=mp.zeros(Mh)
for j in range(Mh):
    for k in range(Mh):
        pts=bpts(hb[j],hb[k])
        Ah[j,k]=(1/(2*mp.pi))*2*mp.quad(lambda r:F(hb[j],r)*F(hb[k],r)*Omega(r),pts)
evh,_=mp.eigsy((Ah+Ah.T)/2)
evhl=[evh[i] for i in range(Mh)]
print(f"  HIGH-BAND-ONLY (r>{mp.nstr(CROSS,4)}, pure 'primitive') form: min eig={mp.nstr(min(evhl),4)}  => {'POSITIVE-DEFINITE (primitive part is definite, as Hodge-index predicts)' if min(evhl)>0 else 'not pos-def'}",flush=True)
print("READ: nneg small & low-band-localized + high-band pos-def ⇒ arch form = (low-r fiber, indefinite) ⊕ (high-r primitive, positive) = the Hodge-index decomposition at the arch place; C-C's compression removes the fiber part.",flush=True)
print("DONE",flush=True)
