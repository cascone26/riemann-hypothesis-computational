#!/usr/bin/env python3
"""§302: which arithmetic data forces the Weil positivity? Build each prime-power's contribution matrix P_p ONCE,
so Q(s) = A - sum_p s_p P_p is cheap. Baseline s_p=1 => true zeta. Measure the spectral GAP (smallest POSITIVE
Weil eigenvalue, ~0.0177 from §300) and its sensitivity d(gap)/d(s_p) per prime => fingerprint of the forcing.
Also: does perturbing the primes drive the bulk NEGATIVE (positivity is arithmetic-tuned, not generic)? REPORT."""
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
prime_powers=sorted(Lam)
lo,hi,M=14,38,9
centers=[mp.mpf(lo)+i*(mp.mpf(hi-lo)/(M-1)) for i in range(M)]
print(f"build A + per-prime-power matrices, band[{lo},{hi}] {M} probes, {len(prime_powers)} prime powers",flush=True)
A=mp.zeros(M)
Pp={n:mp.zeros(M) for n in prime_powers}   # contribution of prime power n (before -sign)
for j in range(M):
    for k in range(M):
        pts=bpts(centers[j],centers[k])
        arch=(1/(2*mp.pi))*2*mp.quad(lambda r:F(centers[j],r)*F(centers[k],r)*Omega(r),pts)
        pole=mp.re(Fc(centers[j],1j/2)*Fc(centers[k],1j/2))*2
        A[j,k]=pole+arch
        for n in prime_powers:
            Pp[n][j,k]=2*Lam[n]/mp.sqrt(n)*(1/(2*mp.pi))*2*mp.quad(lambda r:F(centers[j],r)*F(centers[k],r)*mp.cos(mp.log(n)*r),pts)
    print(f"  row {j}",flush=True)
def Qmat(scale):  # scale: dict n->s_n
    Q=A.copy()
    for n in prime_powers: Q-=scale.get(n,mp.mpf(1))*Pp[n]
    return Q
def gap(Q):  # smallest positive eigenvalue (the RH-in-band gap) + smallest eigenvalue overall
    e,_=mp.eigsy((Q+Q.T)/2); el=[e[i] for i in range(M)]
    pos=[x for x in el if x>mp.mpf('1e-9')]
    return (min(pos) if pos else mp.mpf(0)), el[0]
base={n:mp.mpf(1) for n in prime_powers}
g0,min0=gap(Qmat(base))
print(f"BASELINE (true zeta): spectral gap (smallest +eig) = {mp.nstr(g0,6)}, min eig = {mp.nstr(min0,4)}",flush=True)
# sensitivity: bump each prime's weight by +5% and see gap shift
print("per-prime sensitivity d(gap) for +5% weight bump (which primes HOLD UP positivity):",flush=True)
sens=[]
for p in [n for n in prime_powers if n<30]:
    sc=dict(base); sc[p]=mp.mpf('1.05')
    gp,mnp=gap(Qmat(sc))
    sens.append((p,gp-g0,mnp))
    print(f"  n={p}: d(gap)={mp.nstr(gp-g0,3)}  new_min_eig={mp.nstr(mnp,3)}",flush=True)
# stress: scale ALL primes by (1+eps), find eps where bulk goes negative
print("stress: uniform prime overweight (1+eps) -> when does a bulk eig go negative?",flush=True)
for eps in ['0.02','0.05','0.10','0.20']:
    sc={n:1+mp.mpf(eps) for n in prime_powers}
    g,mn=gap(Qmat(sc))
    print(f"  eps={eps}: gap={mp.nstr(g,4)} min_eig={mp.nstr(mn,4)} {'<-- POSITIVITY BROKEN' if mn<mp.mpf('-1e-6') else ''}",flush=True)
print("DONE",flush=True)
