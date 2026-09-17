"""One-level density of the quadratic Dirichlet (SYMPLECTIC) family via fast AFE. Decisive:
does the scaled low-zero density VANISH near x=0 (symplectic W_Sp=1-sinc(2πx)) vs flat (unitary)?"""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastL import is_fundamental, low_zeros_fast, kronecker
import math

def run(dmax=500, ndisc=300, m=8):
    fds=[d for d in range(5,dmax) if is_fundamental(d)]+[d for d in range(-3,-dmax,-1) if is_fundamental(d)]
    scaled=[]; nd=0
    for d in fds:
        try: zs=low_zeros_fast(d,m=m)
        except Exception: continue
        Ld=math.log(abs(d)/math.pi)/(2*math.pi)
        for g in zs: scaled.append(g*Ld)
        nd+=1
        if nd>=ndisc: break
    scaled=np.array(scaled); np.save('spectral/onelevel_scaled.npy',scaled)
    print(f"ONE-LEVEL DENSITY (SYMPLECTIC family): {nd} quadratic L-functions, {len(scaled)} scaled zeros")
    bins=np.linspace(0,4,33); h,e=np.histogram(scaled,bins=bins)
    c=(e[:-1]+e[1:])/2; w=np.diff(e); dens=h/(nd*w)
    bulk=dens[(c>1.5)&(c<3.5)].mean()
    print(f"  bulk density [1.5,3.5] = {bulk:.3f} (KS normalizes ->1); rescaling by it.")
    print(f"  {'x':>6} {'emp/bulk':>9} {'W_Sp':>7} {'W_U':>5} {'W_O(even)':>9}")
    for i in range(len(c)):
        if c[i]<3.0:
            x=c[i]
            Wsp=1-math.sin(2*math.pi*x)/(2*math.pi*x)
            Wo =1+math.sin(2*math.pi*x)/(2*math.pi*x)
            bar='#'*int(max(0,dens[i]/bulk)*18)
            print(f"  {x:6.2f} {dens[i]/bulk:9.3f} {Wsp:7.3f} {1.0:5.1f} {Wo:9.3f}  {bar}")
    near=dens[c<0.3].mean()/bulk
    Wsp02=1-math.sin(2*math.pi*0.15)/(2*math.pi*0.15)
    print(f"\n  NEAR-0 (x<0.3) empirical W/bulk = {near:.3f}")
    print(f"  Symplectic predicts W_Sp(0.15) = {Wsp02:.3f} (SUPPRESSED); Unitary=1.0; Orthogonal(even) ENHANCED>1")
    verdict = ('SYMPLECTIC signature CONFIRMED — center repulsion' if near<0.75 else
               'no clear suppression at this N' )
    print(f"  VERDICT: {verdict}")

if __name__=="__main__": run()
