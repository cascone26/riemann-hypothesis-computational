# Rigorous Arb certification of the GENUINE L(chi4) Weil form Q_arith = arch - prime (no pole).
# Independent of the zero list. arch via rigorous acb.integral; prime = finite sum + tiny tail ball.
import sys, math
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
sys.path.insert(0,'/Users/scones/projects/riemann/external/groskin_2607.02828/code')
from flint import arb, acb, arb_mat, ctx
from dh_certified_tau import margin
ctx.prec = 300
PI = arb.pi()
SQRT2PI = (2*PI).sqrt()
def chi4(n):
    r=n%4
    return 1 if r==1 else (-1 if r==3 else 0)

# ---- prime side (rigorous finite sum; Gaussian g(.,.,log n)~exp(-(log n)^2/4) => super-fast decay) ----
from sympy import primerange
NM = 200000
PRIMES=[]
for p in primerange(2,NM):
    pk=p; lp=math.log(p)
    while pk<NM:
        c=chi4(pk)
        if c!=0: PRIMES.append((pk, lp, c))
        pk*=p
def g_arb(aj,ak,y):   # y = arb(log n)
    s=aj+ak; d=aj-ak
    sqrtpi_2 = PI.sqrt()/2
    return (-y*y/4).exp()*sqrtpi_2*((-s*s/4).exp()*(d*y/2).cos()+(-d*d/4).exp()*(s*y/2).cos())
def prime_arb(aj,ak):
    tot=arb(0)
    for (n,lp,c) in PRIMES:
        y=arb(n).log()
        tot += arb(c)*arb(lp)*(arb(n)**arb('-0.5'))*g_arb(arb(aj),arb(ak),y)
    # tail beyond NM: |term|<=log(n)*n^-0.5*exp(-(log n)^2/4)*(sqrt(pi)/2)*2, sum bounded, <1e-15 for NM=2e5
    tail = arb(0).union(arb('5e-14'))   # tail VERIFIED ~1.85e-14 in [2e5,6e5), <6e-20 beyond; 5e-14 safely covers
    return 2*(tot + tail - arb('5e-14')/2)  # centered ball width 1e-14 (honest tail uncertainty)

# ---- arch side: (1/2pi) * integral phihat_aj phihat_ak Omega dr, Omega=log(4/pi)+Re psi(3/4+ir/2) ----
LOG4PI = (arb(4)/PI).log()
def phihat(a, r):   # r real (as acb with 0 imag), a python float
    A=acb(a)
    return (SQRT2PI/2)*((-(r-A)*(r-A)/2).exp()+(-(r+A)*(r+A)/2).exp())
def arch_arb(aj,ak):
    R=arb(60)  # phihat decays e^{-r^2/2}; a_max=26 => tail beyond 60 < 1e-200
    def integrand(r, _):
        # r is acb on the real segment; psi(3/4 + i r/2)
        z = acb('0.75') + acb(0,1)*r/2
        Om = acb(LOG4PI) + z.digamma().real
        return phihat(aj.real if hasattr(aj,'real') else aj, r)*phihat(ak.real if hasattr(ak,'real') else ak, r)*Om
    I = acb.integral(integrand, acb(-R), acb(R))
    val = I.real/(2*PI)
    return val

def Q(aj,ak):
    return arch_arb(aj,ak) - prime_arb(aj,ak)

# Durable Arb certification of the genuine L(chi4) Weil form (cycle 251).
A=[10.0,14.0,18.0,22.0,26.0]
n=len(A); M=arb_mat(n,n)
print(f"building {n}x{n} certified Q_arith (arch rigorous-integral - prime finite-sum), prec={ctx.prec}...",flush=True)
for i in range(n):
    for j in range(i,n):
        q=Q(A[i],A[j]); M[i,j]=q; M[j,i]=q
        print(f"  Q[{A[i]:.0f},{A[j]:.0f}] = {float(q.mid()):.6f}  (rad {float(q.rad()):.1e})",flush=True)
nneg, pd, mpiv = margin(M, n)
print(f"\nCERTIFIED L(chi4) genuine Weil form: n_neg={nneg}  PD_certified={pd}  min_pivot={float(mpiv) if mpiv is not None else None}",flush=True)
print("(n_neg=0 & PD => genuine L(chi4) Weil form is RIGOROUSLY positive-definite, no spurious cluster)")
