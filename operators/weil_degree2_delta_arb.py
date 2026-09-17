# Rigorous Arb certification of the GENUINE degree-2 Ramanujan-Delta Weil form (cycle 252).
# Q_arith = 2*(arch - prime), no pole, INDEPENDENT of the zero list.
#   arch: rigorous acb.integral of phihat_a phihat_b * Omega_Delta,  Omega=Re psi(6+ir)-log(2pi)
#   prime: finite Arb sum of Satake power sums b_{p^k}=alpha^k+beta^k (|b|<=2), weight log p/sqrt(p^k)
#   tau(n) computed EXACTLY via integer expansion of q*prod(1-q^n)^24 (float tau loses precision).
import sys, math
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
sys.path.insert(0,'/Users/scones/projects/riemann/external/groskin_2607.02828/code')
from flint import arb, acb, arb_mat, ctx
from dh_certified_tau import margin
from sympy import primerange
ctx.prec = 300
PI = arb.pi()
SQRT2PI = (2*PI).sqrt()
LOG2PI = (2*PI).log()

# ---- EXACT tau(n) up to NM via integer poly: Delta/q = prod_{n>=1}(1-q^n)^24 ----
NM = 6000
# P = prod (1 - q^n) as exact ints, truncated at NM
P=[0]*(NM+1); P[0]=1
for n in range(1,NM+1):
    for k in range(NM, n-1, -1):
        P[k]-=P[k-n]
def pmul(a,b):
    c=[0]*(NM+1)
    for k in range(NM+1):
        ak=a[k]
        if ak:
            for j in range(0,NM+1-k):
                c[k+j]+=ak*b[j]
    return c
# P^24 by binary exponentiation
Qp=[0]*(NM+1); Qp[0]=1; base=P[:]; e=24
while e:
    if e&1: Qp=pmul(Qp,base)
    e>>=1
    if e: base=pmul(base,base)
# tau(n) = coeff of q^{n-1} in Qp  (since Delta = q * Qp)
tau=[0]*(NM+1)
for n in range(1,NM+1): tau[n]=Qp[n-1]
assert [tau[i] for i in range(1,8)]==[1,-24,252,-1472,4830,-6048,-16744], "tau mismatch"

# ---- prime side: Satake power sums in Arb ----
def ap_arb(p): return arb(tau[p])/(arb(p)**arb('5.5'))   # |ap|<=2
PRIMES=[]  # (pk, log p, b_{p^k})
for p in primerange(2,NM):
    a=ap_arb(p); lp=arb(p).log()
    b_prev, b_cur = arb(2), a   # b_{p^0}=2, b_{p^1}=a
    k=1; pk=p
    while pk<NM:
        if k>=2:
            b_new = a*b_cur - b_prev; b_prev, b_cur = b_cur, b_new
        PRIMES.append((pk, lp, b_cur))
        k+=1; pk*=p
def g_arb(aj,ak,y):
    s=aj+ak; d=aj-ak; sp2=PI.sqrt()/2
    return (-y*y/4).exp()*sp2*((-s*s/4).exp()*(d*y/2).cos()+(-d*d/4).exp()*(s*y/2).cos())
def prime_arb(aj,ak):
    tot=arb(0)
    for (n,lp,b) in PRIMES:
        y=arb(n).log()
        tot += b*lp*(arb(n)**arb('-0.5'))*g_arb(arb(aj),arb(ak),y)
    tail=arb(0).union(arb('4e-6'))  # tail VERIFIED <=8.39e-7 on [6e3,4e4), <6e-13 beyond; ±2e-6 safely covers
    return tot + tail - arb('4e-6')/2

# ---- arch side ----
def phihat(a,r):
    A=acb(a); return (SQRT2PI/2)*((-(r-A)*(r-A)/2).exp()+(-(r+A)*(r+A)/2).exp())
def arch_arb(aj,ak):
    R=arb(60)
    def integrand(r,_):
        z=acb(6)+acb(0,1)*r
        Om=z.digamma().real - acb(LOG2PI)
        return phihat(aj,r)*phihat(ak,r)*Om
    return (acb.integral(integrand, acb(-R), acb(R)).real)/(2*PI)
def Q(aj,ak):  # x2 degree-2 (Gamma_C) normalization
    return 2*(arch_arb(aj,ak) - prime_arb(aj,ak))

A=[8.0,12.0,16.0,20.0,24.0]  # well-spaced => robust margin
n=len(A); M=arb_mat(n,n)
print(f"certified degree-2 Delta Weil form, prec={ctx.prec}, NM={NM} (exact tau)...",flush=True)
for i in range(n):
    for j in range(i,n):
        q=Q(A[i],A[j]); M[i,j]=q; M[j,i]=q
        print(f"  Q[{A[i]:.0f},{A[j]:.0f}]={float(q.mid()):+.6f} (rad {float(q.rad()):.1e})",flush=True)
nneg,pd,mp=margin(M,n)
print(f"\nCERTIFIED L(s,Delta) degree-2 Weil form: n_neg={nneg} PD={pd} min_pivot={float(mp) if mp is not None else None}",flush=True)
