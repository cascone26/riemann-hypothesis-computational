#!/usr/bin/env python3
"""§300 forward move: restrict the prime form P to the de Branges / annihilator subspace (where Q=A-P ~ 0,
i.e. the near-saturating subspace, dim ~ N - #in-band-zeros) and inspect its structure. RH-in-band lives
ENTIRELY here (§299 gap). Question: does P|_sub have recognizable arithmetic structure (a small kernel whose
PSD-ness = RH-in-band)? Charisma: REPORT (numerics)."""
import mpmath as mp
mp.mp.dps=18
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
for pp in range(2,60):
    if all(pp%d for d in range(2,int(pp**.5)+1)):
        q=pp
        while q<60: Lam[q]=mp.log(pp);q*=pp
lo,hi,M=14,38,9
centers=[mp.mpf(lo)+i*(mp.mpf(hi-lo)/(M-1)) for i in range(M)]
print(f"build A(pole+arch),P(prime) band[{lo},{hi}] {M} probes",flush=True)
A=mp.zeros(M);P=mp.zeros(M)
for j in range(M):
    for k in range(M):
        pts=bpts(centers[j],centers[k])
        arch=(1/(2*mp.pi))*2*mp.quad(lambda r:F(centers[j],r)*F(centers[k],r)*Omega(r),pts)
        pole=mp.re(Fc(centers[j],1j/2)*Fc(centers[k],1j/2))*2
        pr=2*sum(Lam[n]/mp.sqrt(n)*(1/(2*mp.pi))*2*mp.quad(lambda r:F(centers[j],r)*F(centers[k],r)*mp.cos(mp.log(n)*r),pts) for n in sorted(Lam))
        A[j,k]=pole+arch;P[j,k]=pr
    print(f"  row {j}",flush=True)
Q=A-P
eQ,VQ=mp.eigsy((Q+Q.T)/2)
eQl=[eQ[i] for i in range(M)]
print(f"Q eigs: {[mp.nstr(x,4) for x in eQl]}",flush=True)
# annihilator/de Branges subspace = eigenvectors with |Q-eig| below the gap (§299: gap after 2-3 tiny ones)
thr=mp.mpf('0.05')
idx=[i for i in range(M) if abs(eQl[i])<thr]
print(f"de Branges subspace dim = {len(idx)} (indices {idx}); band has ~6 zeta zeros, N={M} => expect ~{M}-6",flush=True)
U=mp.zeros(M,len(idx))
for c,i in enumerate(idx):
    for r_ in range(M): U[r_,c]=VQ[r_,i]
# restrict P and A to this subspace
Psub=U.T*P*U; Asub=U.T*A*U; Qsub=U.T*Q*U
d=len(idx)
ePs=[mp.eigsy((Psub+Psub.T)/2)[0][i] for i in range(d)]
eQs=[mp.eigsy((Qsub+Qsub.T)/2)[0][i] for i in range(d)]
print(f"  P restricted to de Branges subspace, eigenvalues: {[mp.nstr(x,5) for x in ePs]}",flush=True)
print(f"  Q restricted (RH-in-band = these >=0): {[mp.nstr(x,4) for x in eQs]}  min={mp.nstr(min(eQs),5)}",flush=True)
print(f"  Psub matrix (the localized crux, {d}x{d}):",flush=True)
for i in range(d):
    print("    "+", ".join(mp.nstr(Psub[i,j],5) for j in range(d)),flush=True)
# structure hunt: is Psub ~ a rank-structured / Gram matrix of {evaluations at zeros}?
# ratio of Psub to Asub entrywise (should be ~1 if P~A on subspace i.e. saturation)
print(f"  Psub/Asub entrywise (diag): {[mp.nstr(Psub[i,i]/Asub[i,i],5) for i in range(d)]} (->1 = saturated)",flush=True)
print("DONE",flush=True)
