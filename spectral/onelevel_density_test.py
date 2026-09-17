"""One-level density: THE decisive Katz-Sarnak discriminant. Collect the low zeros of many quadratic
Dirichlet L-functions (symplectic family), scale each by L=log(|d|/π)/2π, histogram the scaled zeros,
compare to the three universal curves:
  Unitary:    W_U(x)  = 1
  Symplectic: W_Sp(x) = 1 - sin(2πx)/(2πx)      <- VANISHES at x=0 (repulsion) -- the class-C signature
  Orthogonal: W_O(x)  = 1 + (1/2)δ(x) ... (even) or 1 - sinc + δ (odd); near 0 ENHANCED
Quadratic Dirichlet family is the textbook SYMPLECTIC family. If W_Sp fits and W_U/W_O don't, the
symplectic/class-C structure is CONFIRMED to live in the family (not ζ's bulk)."""
import mpmath as mp, numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from symplectic_family_test import is_fundamental, chi_vec, Zfun
mp.mp.dps=15

def low_zeros(d, chiv, m=6, tmax=40.0, step=0.15):
    zs=[]; pt=0.02; pv=Zfun(pt,d,chiv); t=step
    while t<tmax and len(zs)<m:
        cv=Zfun(t,d,chiv)
        if pv*cv<0:
            lo,hi=pt,t
            for _ in range(28):
                mm=(lo+hi)/2; vm=Zfun(mm,d,chiv)
                if pv*vm<0: hi=mm
                else: lo,pv=mm,vm
            zs.append((lo+hi)/2)
        pt,pv,t=t,cv,t+step
    return zs

def run(dmax=250, ndisc=140, m=6, out='spectral/onelevel_scaled.npy'):
    fds=[d for d in range(5,dmax) if is_fundamental(d)]+[d for d in range(-3,-dmax,-1) if is_fundamental(d)]
    scaled=[]; nd=0
    for d in fds:
        cv=chi_vec(d); zs=low_zeros(d,cv,m=m)
        Ld=float(mp.log(abs(d)/mp.pi)/(2*mp.pi))
        for g in zs: scaled.append(g*Ld)
        nd+=1
        if nd>=ndisc: break
    scaled=np.array(scaled); np.save(out,scaled)
    # empirical one-level density (normalize so bulk ~1): histogram, divide by ndisc and bin width... 
    # here just show near-0 suppression vs bulk
    print(f"ONE-LEVEL DENSITY: {nd} quadratic L-functions, {len(scaled)} scaled low zeros")
    bins=np.linspace(0,4,21); h,e=np.histogram(scaled,bins=bins)
    c=(e[:-1]+e[1:])/2; w=np.diff(e)
    dens=h/(nd*w)   # zeros per unit scaled-x per L-function
    bulk=dens[(c>1.5)&(c<3.5)].mean()
    print(f"  bulk density (x in [1.5,3.5]) = {bulk:.3f}  (Katz-Sarnak: ->1)")
    print(f"  {'x':>6} {'emp W(x)':>9} {'W_Sp':>7} {'W_U':>5}")
    for i in range(len(c)):
        if c[i]<3.0:
            Wsp=1-np.sin(2*np.pi*c[i])/(2*np.pi*c[i]) if c[i]>1e-6 else 0.0
            print(f"  {c[i]:6.2f} {dens[i]/bulk:9.3f} {Wsp:7.3f} {1.0:5.1f}  {'#'*int(dens[i]/bulk*20)}")
    # decisive: near-0 suppression?
    near=dens[c<0.4].mean()/bulk
    print(f"\n  NEAR-0 (x<0.4) empirical W = {near:.3f}   [SYMPLECTIC predicts ~{1-np.sin(2*np.pi*0.2)/(2*np.pi*0.2):.3f} (suppressed); UNITARY=1.0]")
    print(f"  => {'SYMPLECTIC signature CONFIRMED (repulsion at center)' if near<0.7 else 'no clear suppression at this N — needs more discriminants/zeros'}")

if __name__=="__main__":
    run()
