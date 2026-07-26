#!/usr/bin/env python3
"""FORWARD LEVER (§299): the reality<->correspondence link probed as a generalized eigenproblem.

Weil: Q(f) = pole(f) + arch(f) - prime(f) >= 0  <=>  RH.
Split the metric: A := pole+arch (the POSITIVE background the archimedean+pole terms set up),
                  P := prime (the oscillatory arithmetic correction).
Then Q>=0  <=>  P <= A in PSD order  <=>  lambda_max of the pencil (P; A) <= 1.

§290 established sup = EXACTLY 1 (zero uniform slack). The NEW question this run asks:
  does the generalized spectrum lambda_i(P;A) have a GAP below 1 (arch dominates except in the
  de Branges E-direction) or does it ACCUMULATE at 1 (marginal in many directions)?
A gap => the arch background genuinely forces positivity, only ONE direction is tight (=> a
lever: bound that one direction). Accumulation => the '0 slack' is pervasive (harder).
Also: track the TOP eigenvector as the band widens -> is it stably the E-direction?

Charisma: REPORT (numerical spectral exploration). No proof claimed."""
import mpmath as mp
mp.mp.dps=16
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

def build(centers):
    M=len(centers)
    A=mp.zeros(M); P=mp.zeros(M)   # A = pole+arch (positive background), P = prime
    for j in range(M):
        for k in range(M):
            pts=bpts(centers[j],centers[k])
            arch=(1/(2*mp.pi))*2*mp.quad(lambda r:F(centers[j],r)*F(centers[k],r)*Omega(r),pts)
            pole=mp.re(Fc(centers[j],1j/2)*Fc(centers[k],1j/2))*2
            pr=2*sum(Lam[n]/mp.sqrt(n)*(1/(2*mp.pi))*2*mp.quad(lambda r:F(centers[j],r)*F(centers[k],r)*mp.cos(mp.log(n)*r),pts) for n in sorted(Lam))
            A[j,k]=pole+arch
            P[j,k]=pr
        print(f"    row {j}/{M-1}",flush=True)
    return A,P

def gen_spectrum(A,P):
    # symmetric generalized eig P v = lam A v  ->  Cholesky A=L L^T, standard eig of L^-1 P L^-T
    L=mp.cholesky(A)
    Linv=L**-1
    Mmat=Linv*P*Linv.T
    Mmat=(Mmat+Mmat.T)/2
    ev,V=mp.eigsy(Mmat)
    return ev,V,L

def qform(Mtx,v):
    M=len(v); s=mp.mpf(0)
    for i in range(M):
        for j in range(M): s+=v[i]*Mtx[i,j]*v[j]
    return s

for band in [(14,38,9)]:
    lo,hi,M=band
    step=(hi-lo)/(M-1)
    centers=[mp.mpf(lo)+i*step for i in range(M)]
    print(f"BAND [{lo},{hi}], {M} probes, step {mp.nstr(step,3)}:",flush=True)
    A,P=build(centers)
    # --- SCALE sanity: diagonal magnitudes of A (arch+pole) vs P (prime) ---
    print(f"  scale: mean|A_ii|={mp.nstr(sum(abs(A[i,i]) for i in range(M))/M,4)}  mean|P_ii|={mp.nstr(sum(abs(P[i,i]) for i in range(M))/M,4)}",flush=True)
    # --- direct symmetric eigenspectra of A and P separately (in the plain probe basis) ---
    eA,_=mp.eigsy((A+A.T)/2); eP,_=mp.eigsy((P+P.T)/2)
    print(f"  A eigs (arch+pole, want >0): {[mp.nstr(x,4) for x in eA]}",flush=True)
    print(f"  P eigs (prime form):        {[mp.nstr(x,4) for x in eP]}",flush=True)
    # --- Weil form Q = A - P : its spectrum is the real RH object ---
    Q=A-P; eQ,VQ=mp.eigsy((Q+Q.T)/2)
    print(f"  Q=A-P eigs (Weil form, RH-relevant, want >=0): {[mp.nstr(x,4) for x in eQ]}",flush=True)
    print(f"    min Weil eig = {mp.nstr(eQ[0],6)}  (near 0 => near-saturating direction present)",flush=True)
    # --- generalized pencil lambda(P;A) via Cholesky ---
    ev,V,L=gen_spectrum(A,P)
    evl=[ev[i] for i in range(len(list(ev)))]   # eigsy returns a matrix; list it explicitly (neg-index was buggy)
    print(f"  generalized eigs lambda(P;A), RH<=> all<=1, §290 sup->1:",flush=True)
    print(f"    {[mp.nstr(x,6) for x in evl]}",flush=True)
    lam_max=evl[-1]
    gap=lam_max-evl[-3]   # gap from top pair down to the bulk
    print(f"  lambda_max = {mp.nstr(lam_max,10)}   spectral GAP (top2 -> bulk) = {mp.nstr(gap,5)}",flush=True)
    # --- the DIRECT Rayleigh ratio P(f)/A(f) in the min-Weil (saturating) direction ---
    vsat=[VQ[i,0] for i in range(M)]
    Pf=qform(P,vsat); Af=qform(A,vsat); Qf=qform(Q,vsat)
    print(f"  min-Weil direction: A(f)={mp.nstr(Af,5)}  P(f)={mp.nstr(Pf,5)}  Q(f)={mp.nstr(Qf,4)}  ratio P/A={mp.nstr(Pf/Af,5)} (->1 = saturation)",flush=True)
    def f(x): return sum(vsat[j]*F(centers[j],x) for j in range(M))
    zz=[mp.mpf(z) for z in ['14.134725','21.022040','25.010858','30.424876','32.935062','37.586178'] if mp.mpf(z)<hi-2]
    nrm=mp.sqrt(sum(abs(x)**2 for x in vsat))
    print(f"  min-Weil eigvec |f(gamma)|/||a|| at zeta zeros: {[mp.nstr(abs(f(g))/nrm,2) for g in zz]} (small=>zero-vanishing=de Branges dir)",flush=True)
    print("",flush=True)
print("DONE",flush=True)
