#!/usr/bin/env python3
"""
ldlt_structure_mine.py -- ATTACK A1: mine the certified LDL^T / Cholesky of tau for a
POSITIVITY-PROOF pattern.

tau = A + B - P is certified PSD (Groskin). Its LDL^T pivots d_i are all >0 (that IS the
certificate). Question nobody has asked: do the pivots d_i(N), or the Cholesky factor, have
a RECOGNIZABLE STRUCTURE -- a closed form, a clean recurrence, a manifest-positivity form --
that would explain WHY they stay positive uniformly in N? A structural law here is a
constructive proof seed, not a confirmation.

We build tau at high precision via the Groskin closed forms (mpmath mirror -- same entries
the Arb tool certifies), do an exact-ish LDL^T, and study:
  1. the pivot spectrum d_i: how many are ~0 (near-kernel), the gap, the bulk law.
  2. pivot STABILITY across N: does d_i (for fixed i from the top) converge as N grows?
  3. the near-kernel pivots: their magnitude vs N -- bounded below? clean power law?
  4. whether the Cholesky factor L has decaying/structured off-diagonals (banded? Toeplitz-ish?).
"""
import mpmath as mp
import numpy as np
import json, sys, time

mp.mp.dps = 60
TWO_PI = 2*mp.pi

def prime_powers_up_to(c):
    primes, x = [], 2
    while x <= c:
        if all(x % p for p in primes): primes.append(x)
        x += 1
    out = []
    for p in primes:
        q = p
        while q <= c: out.append((q, p)); q *= p
    return out

def geom_sums(n, L):
    w = 2*mp.pi*n/L; w2 = w*w
    gS=gCC=gX1=gX2=mp.mpf(0); k=0
    while True:
        ck = mp.mpf(2*k)+mp.mpf('0.5'); e = mp.e**(-ck*L); den = ck*ck+w2
        gS += e/den
        if n != 0: gCC += e*w2/(ck*den)
        gX1 += e*ck/den; gX2 += e*(ck*ck-w2)/(den*den)
        if e < mp.mpf(10)**(-70) and k>2: break
        k += 1
    return gS,gCC,gX1,gX2

def build_tau(N, c):
    L = mp.log(c); PI = mp.pi
    psi_q = mp.digamma(mp.mpf(1)/4)
    S=[mp.mpf(0)]*(N+1); CC=[mp.mpf(0)]*(N+1); XC=[mp.mpf(0)]*(N+1)
    for n in range(N+1):
        w = 2*PI*n/L; z = mp.mpf(1)/4 + 1j*PI*n/L
        psi = mp.digamma(z); psi1 = mp.polygamma(1, z)
        gS,gCC,gX1,gX2 = geom_sums(n, L)
        if n==0: S[n]=mp.mpf(0); CC[n]=mp.mpf(0)
        else: S[n]=mp.im(psi)/2 - w*gS; CC[n]=-(mp.re(psi)-psi_q)/2 + gCC
        XC[n]=mp.re(psi1)/4 - L*gX1 - gX2
    sp2 = 16*PI*PI; l2=L*L
    pref02 = 32*L*mp.sinh(L/4)**2
    eL=mp.e**L; kappa = mp.log(4*PI*(eL-1)/(eL+1)) + mp.euler
    U=mp.e**(L/2); J = -2*mp.log(U+1)+mp.log(U*U+1)+2*mp.atan(U)+mp.log(2)-PI/2
    pdata = prime_powers_up_to(c)
    wts=[mp.log(p)*mp.mpf(q)**mp.mpf('-0.5') for (q,p) in pdata]
    pos=[mp.log(q) for (q,p) in pdata]
    def Ssg(nn): return S[nn] if nn>=0 else -S[-nn]
    D=2*N+1; tau=mp.matrix(D,D)
    for i in range(D):
        n=i-N
        for j in range(i,D):
            m=j-N
            num=l2-sp2*m*n; den=(l2+sp2*m*m)*(l2+sp2*n*n)
            W02=pref02*num/den
            if n==m: WR = kappa+2*CC[abs(n)]+J-(2/L)*XC[abs(n)]
            else: WR=(Ssg(m)-Ssg(n))/(PI*(n-m))
            Wp=mp.mpf(0)
            for idx in range(len(wts)):
                y=pos[idx]
                if n==m: qv=2*(1-y/L)*mp.cos(2*PI*n*y/L)
                else: qv=(mp.sin(2*PI*m*y/L)-mp.sin(2*PI*n*y/L))/(PI*(n-m))
                Wp += wts[idx]*qv
            v=W02-WR-Wp; tau[i,j]=v; tau[j,i]=v
    return tau, D

def ldlt_pivots(tau, D):
    """LDL^T; return pivots d[0..D-1] (mpf)."""
    d=[mp.mpf(0)]*D; Lf=[[mp.mpf(0)]*D for _ in range(D)]
    for i in range(D):
        s=tau[i,i]
        for k in range(i): s-=Lf[i][k]*Lf[i][k]*d[k]
        d[i]=s
        for j in range(i+1,D):
            t=tau[j,i]
            for k in range(i): t-=Lf[j][k]*Lf[i][k]*d[k]
            Lf[j][i]=t/d[i]
    return d, Lf

if __name__=="__main__":
    c = int(sys.argv[1]) if len(sys.argv)>1 else 14
    Ns = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv)>2 else [16,24,32,48,64]
    out={"c":c,"dps":mp.mp.dps,"runs":{}}
    print(f"[ldlt-mine] c={c}, N={Ns}, dps={mp.mp.dps}", flush=True)
    for N in Ns:
        t0=time.time()
        tau,D=build_tau(N,c)
        d,Lf=ldlt_pivots(tau,D)
        dv=[float(x) for x in d]
        dvpos=[x for x in dv if x>0]; dvneg=[x for x in dv if x<=0]
        sd=sorted(dv)
        # near-zero pivots (the RH content) vs bulk
        nnz=sum(1 for x in dv if abs(x)<1e-6)
        print(f"\n=== N={N} (D={D}, {time.time()-t0:.0f}s) ===", flush=True)
        print(f"  pivots: {len(dvpos)} pos, {len(dvneg)} <=0   min={min(dv):.3e}  max={max(dv):.3e}")
        print(f"  smallest 6 pivots: {[f'{x:.3e}' for x in sd[:6]]}")
        print(f"  #|pivot|<1e-6 = {nnz}   (near-kernel size)")
        # bulk pivot law: median and a few percentiles
        print(f"  pivot percentiles [10,50,90]: {[f'{np.percentile(dv,p):.3e}' for p in (10,50,90)]}")
        out["runs"][N]={"D":D,"n_pos":len(dvpos),"n_nonpos":len(dvneg),
                        "min":min(dv),"max":max(dv),"smallest6":sd[:6],
                        "near_kernel_count":nnz,
                        "pct10_50_90":[float(np.percentile(dv,p)) for p in (10,50,90)]}
    # cross-N: does the k-th smallest pivot converge?
    print("\n" + "="*70)
    print("k-th SMALLEST pivot vs N (convergence => structural limit):")
    for k in range(6):
        row=[f"N={N}:{out['runs'][N]['smallest6'][k]:.3e}" for N in Ns]
        print(f"  k={k}: " + "  ".join(row))
    p="/Users/scones/projects/riemann/results/ldlt_structure_mine.json"
    json.dump(out,open(p,"w"),indent=2)
    print(f"\nsaved -> {p}")
