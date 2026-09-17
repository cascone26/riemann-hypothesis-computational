"""Katz-Sarnak symplectic-family test: lowest-zero distribution of quadratic Dirichlet L-functions.
Symplectic family => density of the normalized lowest zero is SUPPRESSED near 0 (repulsion from center);
unitary => finite at 0. Decides where class-C/symplectic structure lives (a family, not ζ)."""
import mpmath as mp, numpy as np
from sympy import jacobi_symbol
mp.mp.dps = 15

def sqfree(n):
    n=abs(n); k=2
    while k*k<=n:
        if n%(k*k)==0: return False
        k+=1
    return n!=0

def is_fundamental(d):
    if d==1: return False
    if d%4==1: return sqfree(d)
    if d%4==0:
        m=d//4
        return (m%4 in (2,3)) and sqfree(m)
    return False

def kronecker(a,n):
    # full Kronecker symbol (a|n), n>0
    if n==0: return 1 if abs(a)==1 else 0
    res=1; 
    while n%2==0:
        n//=2
        if a%2==0: return 0
        if a%8 in (3,5): res=-res
    a%=n
    if n==1: return res
    return res*int(jacobi_symbol(a,n))

def chi_vec(d):
    D=abs(d)
    return [kronecker(d,n) for n in range(D)]  # χ_d(n)=(d|n)

def Zfun(t,d,chiv):
    D=abs(d); s=mp.mpf('0.5')+1j*t
    L=mp.mpf(0)
    for a in range(1,D):
        c=chiv[a]
        if c: L+=c*mp.zeta(s,mp.mpf(a)/D)
    L*=D**(-s)
    even=(chiv[D-1]==1)   # χ(-1)=(d|D-1)=(d|-1); even char if +1
    ae=0 if even else 1
    g=(D/mp.pi)**(s/2)*mp.gamma((s+ae)/2)
    return float(mp.re(g*L/mp.fabs(g)))

def first_zero(d,chiv,tmax=14.0,step=0.2):
    pt=0.02; pv=Zfun(pt,d,chiv); t=step
    while t<tmax:
        cv=Zfun(t,d,chiv)
        if pv*cv<0:
            lo,hi=pt,t
            for _ in range(30):
                m=(lo+hi)/2; vm=Zfun(m,d,chiv)
                if pv*vm<0: hi=m
                else: lo,pv=m,vm
            return (lo+hi)/2
        pt,pv,t=t,cv,t+step
    return None

if __name__=="__main__":
    import sys
    # smoke test
    for d in [-3,5,-4,8,-7,12,13]:
        if is_fundamental(d):
            cv=chi_vec(d); g=first_zero(d,cv)
            print(f"  d={d:4d} fundamental, χ(-1)={cv[abs(d)-1]:+d}, first zero γ₁={g}")

def run_batch(dmax=300, target=200, outfile='spectral/symplectic_first_zeros.npy'):
    fds=[d for d in range(5,dmax) if is_fundamental(d)] + [d for d in range(-3,-dmax,-1) if is_fundamental(d)]
    norm=[]; raw=[]
    for d in fds:
        cv=chi_vec(d); g=first_zero(d,cv)
        if g is not None:
            dens=float(mp.log(abs(d)/mp.pi)/(2*mp.pi))
            norm.append(g*dens); raw.append((d,g))
        if len(norm)>=target: break
    norm=np.array(norm); np.save(outfile,norm)
    print(f"BATCH DONE: {len(norm)} normalized first zeros (fundamental discriminants |d|<{dmax})")
    print(f"  mean normalized γ₁ = {norm.mean():.3f}")
    print(f"  P(γ̃₁<0.1)={np.mean(norm<0.1):.3f}  P(<0.2)={np.mean(norm<0.2):.3f}  P(<0.3)={np.mean(norm<0.3):.3f}  P(<0.5)={np.mean(norm<0.5):.3f}")
    print("  Katz-Sarnak: SYMPLECTIC one-level density W(x)=1-sinc(2x) VANISHES at 0 => strong repulsion,")
    print("  P(γ̃₁ small)~0 and mean pushed UP. UNITARY W(x)=1 => finite near 0. ORTHOGONAL(even) => attracted to 0.")
    h,e=np.histogram(norm,bins=np.linspace(0,2.5,11))
    for i in range(len(h)): print(f"    [{e[i]:.2f},{e[i+1]:.2f}): {'#'*h[i]} {h[i]}")
    return norm


# DISPATCH_TAIL
if __name__=="__main__":
    import sys
    if len(sys.argv)>1 and sys.argv[1]=="batch":
        run_batch()
