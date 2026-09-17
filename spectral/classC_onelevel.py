"""Close the RMT<->family loop: does the AZ class-C (BdG) ensemble reproduce the SAME symplectic
one-level density W_Sp(x)=1-sinc(2πx) near E=0 that the quadratic Dirichlet family showed?
If yes, AZ class C IS the correct RMT target for the Hilbert-Polya construction."""
import numpy as np, math
rng=np.random.default_rng(11)

def classC(n):  # 2n x 2n BdG class C: H=[[a,b],[b†,-aᵀ]], a herm, b symmetric
    a=(rng.standard_normal((n,n))+1j*rng.standard_normal((n,n))); a=(a+a.conj().T)/2
    b=(rng.standard_normal((n,n))+1j*rng.standard_normal((n,n))); b=(b+b.T)/2
    return np.vstack([np.hstack([a,b]),np.hstack([b.conj().T,-a.T])])

n=150; T=400
pos=[]
for _ in range(T):
    ev=np.linalg.eigvalsh(classC(n))
    ev=ev[ev>0]                     # positive branch (spectrum symmetric about 0)
    # local mean density near 0: use the smallest ~n/3 eigenvalues; density = count/interval
    ev=np.sort(ev)
    k=len(ev)//3
    dens0=k/ev[k-1]                 # mean density near origin per trial
    pos.extend((ev[:k]*dens0).tolist())   # scale so mean spacing ->1
pos=np.array(pos)
print(f"AZ class-C ensemble: {T} matrices (2x{n}), {len(pos)} scaled near-0 eigenvalues")
bins=np.linspace(0,4,33); h,e=np.histogram(pos,bins=bins); c=(e[:-1]+e[1:])/2; w=np.diff(e)
d=h/(T*w); bulk=d[(c>1.5)&(c<3.5)].mean()
print(f"  {'x':>6} {'C-ens/bulk':>10} {'W_Sp':>7} {'W_U':>5}")
for i in range(len(c)):
    if c[i]<3.0:
        x=c[i]; Wsp=1-math.sin(2*math.pi*x)/(2*math.pi*x)
        print(f"  {x:6.2f} {d[i]/bulk:10.3f} {Wsp:7.3f} {1.0:5.1f}  {'#'*int(max(0,d[i]/bulk)*18)}")
near=d[c<0.3].mean()/bulk
print(f"\n  class-C near-0 (x<0.3) W/bulk = {near:.3f}   [symplectic W_Sp(0.15)=0.14; unitary 1.0]")
print(f"  LOOP: {'CLOSED — AZ class C reproduces the symplectic-family W_Sp (center repulsion)' if near<0.6 else 'partial — check scaling'}")
