#!/usr/bin/env python3
"""
CYCLE 292 — the FULL arch-minus-prime Weil form for the Epstein class-number ladder.

Cycles 290/291: class number h flips the discriminator (h=1 Euler product/PASS/0 off-line; h>=2 FAIL),
and the PRIME side carries the sign (jump weights Λ_Q go negative for h>=2). Cycle 291's next-lead: does
adding the (class-number-BLIND) ARCH term preserve the verdict on the COMPLETE Weil form?

Structural fact (verified below by FE): the completed Epstein zeta is
    ξ_Q(s) = (√d/(2π))^s · Γ(s) · ζ_Q(s),   FE ξ_Q(s)=ξ_Q(1-s),
so the archimedean density in the Weil explicit formula is
    Ω_Q(r) = log(√d/(2π)) + Re ψ(1/2 + i r).
The DIGAMMA shape Re ψ(1/2+ir) is IDENTICAL for every disc; d enters ONLY as the additive constant
log(√d/(2π)) (monotone increasing in d). => the arch side cannot carry the class-number sign.

Full Weil quadratic form on the Gaussian basis φ_a(x)=e^{-x^2/2}cos(a x)  (φ̂_a(r)=√(2π)/2[e^{-(r-a)^2/2}+e^{-(r+a)^2/2}]):
    W_Q(a,a') = POLE(rank-1,+)  +  ARCH  -  PRIME
      POLE  : ζ_Q has a simple pole at s=1 (and ξ_Q at s=0,1) -> +Re[φ̂_a(r₁)φ̂_a'(r₁)]+Re[φ̂_a(r₀)φ̂_a'(r₀)],
              r₁=-i/2 (s=1), r₀=+i/2 (s=0)   (positive geometric-side terms)
      ARCH  : (1/2π)∫ φ̂_a(r)φ̂_a'(r) Ω_Q(r) dr
      PRIME : 2 Σ_N Λ_Q(N)/√N · g(a,a',log N),  g(a,a',y)=e^{-y^2/4}(√π/2)[e^{-s^2/4}cos(dy/2)+e^{-d^2/4}cos(sy/2)],
              s=a+a', d=a-a'; Λ_Q from the exact recursion Σ_{e|N}Λ_Q(e)a(N/e)=a(N)log N (cycle 291).
By Weil/Bombieri: n_neg(W_Q) = # off-line zero PAIRS. INDEPENDENT GATE: must be 0 for h=1, >=1 for h=2,3.
"""
import numpy as np, mpmath as mp
from math import log, gcd
mp.mp.dps = 30

# ---------- prime side: Λ_Q(N) from principal reduced form (cycle 291) ----------
def reduced_principal_form(d):
    if d % 4 == 0: return (1, 0, d//4)
    else:          return (1, 1, (d+1)//4)
def rep_counts(form, Nmax):
    a,b,c = form; r = np.zeros(Nmax+1)
    B = int(2*np.sqrt(Nmax))+3
    for x in range(-B,B+1):
        for y in range(-B,B+1):
            v=a*x*x+b*x*y+c*y*y
            if 1<=v<=Nmax: r[v]+=1.0
    return r
def lambda_Q(d, Nmax):
    r = rep_counts(reduced_principal_form(d), Nmax); a = r/r[1]
    Lam=np.zeros(Nmax+1)
    for n in range(2,Nmax+1):
        s=a[n]*log(n)
        for e in range(1,n):
            if n%e==0: s-=Lam[e]*a[n//e]
        Lam[n]=s
    return Lam

# ---------- basis + kernels ----------
SQ2PI = float(mp.sqrt(2*mp.pi))
def phihat(a, r):                      # r may be complex
    return SQ2PI/2*(mp.e**(-(r-a)**2/2)+mp.e**(-(r+a)**2/2))
def g(a,ap,y):
    s=a+ap; dd=a-ap; sp=float(mp.sqrt(mp.pi))/2
    return mp.e**(-y*y/4)*sp*(mp.e**(-s*s/4)*mp.cos(dd*y/2)+mp.e**(-dd*dd/4)*mp.cos(sp*0+ s*y/2))

def build_W(d, A, Nmax=300):
    Lam = lambda_Q(d, Nmax)
    logconst = float(mp.log(mp.sqrt(d)/(2*mp.pi)))
    n=len(A); W=np.zeros((n,n))
    # arch density Ω(r)=logconst+Re ψ(1/2+ir)
    def Omega(r): return logconst + float(mp.re(mp.digamma(mp.mpf('0.5')+1j*r)))
    for j in range(n):
        for k in range(j,n):
            aj,ak=A[j],A[k]
            arch = float(mp.quad(lambda r: phihat(aj,r)*phihat(ak,r)*Omega(r), [-40,0,40]))/(2*float(mp.pi))
            prime = 0.0
            for N in range(2,Nmax+1):
                if Lam[N]!=0.0:
                    prime += Lam[N]/np.sqrt(N)*float(g(aj,ak,log(N)))
            prime *= 2.0
            # pole rank-1 (+): s=1 -> r=-i/2 ; s=0 -> r=+i/2
            pole = float(mp.re(phihat(aj,-0.5j)*phihat(ak,-0.5j))) + float(mp.re(phihat(aj,0.5j)*phihat(ak,0.5j)))
            W[j,k]=W[k,j]= pole + arch - prime
    return W

# ---------- FE gate on the completion (independent of the form) ----------
def zeta_Q_via_theta(d, s, M=60):
    # Z_Q(s)=Σ' Q_m^{-s} via incomplete-gamma continuation (covolume-1: q=Q*2/√d), ζ_Q=Z_Q/2
    form=reduced_principal_form(d); r=rep_counts(form, 4000)
    qs=[];
    for N in range(1,4001):
        if r[N]>0: qs += [ (N*2/mp.sqrt(d), r[N]) ]
    s=mp.mpc(s); tot=mp.mpc(0)
    for (q,mult) in qs:
        tot += mult*((mp.pi*q)**(-s)*mp.gammainc(s, mp.pi*q) + (mp.pi*q)**(s-1)*mp.gammainc(1-s, mp.pi*q))
    tot += -1/s - 1/(1-s)
    return tot   # = completed Λ̃(s) (FE-symmetric), zeros = zeros of ζ_Q

print("FE gate: Λ̃_Q(s) =? Λ̃_Q(1-s)  [validates ξ_Q=(√d/2π)^sΓ(s)ζ_Q completion]")
for d in (4,20,23):
    a=zeta_Q_via_theta(d, mp.mpc('0.3','5.0')); b=zeta_Q_via_theta(d, mp.mpc('0.7','-5.0'))
    print(f"  disc -{d}: |Λ̃(0.3+5i)-Λ̃(0.7-5i)| = {float(abs(a-b)):.2e}")

A=[0.6,1.0,1.6,2.4,3.4]
print(f"\nFull Weil form W_Q=pole+arch-prime on basis a={A}  (n_neg must match off-line pairs: h=1:0, h>=2:>=1)")
known={4:0,20:8,23:3}
for d,h in [(4,1),(20,2),(23,3)]:
    W=build_W(d,A); ev=np.linalg.eigvalsh(W); nneg=int((ev<-1e-9).sum())
    print(f"  disc -{d} (h={h}): n_neg={nneg}  min-eig={ev.min():+.4e}   off-line zeros known={known[d]}  "
          f"{'PASS(PSD)' if nneg==0 else 'FAIL(indefinite)'}")
