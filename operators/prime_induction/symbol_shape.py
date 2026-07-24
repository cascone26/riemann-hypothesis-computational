# Characterize the first-step symbol S_a(xi) (a=0.40, only prime 2, NO zeta'/zeta poles -> BOUNDED).
# Where/how deep/how wide is the negative dip? Narrow => finite-support concentration bound has a chance; broad => fatal.
import numpy as np, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
import suzuki_screw_lambda_a as S
from scipy.special import digamma
def build_S(a, xis):
    Om=np.real(digamma(0.25+0.5j*xis))-np.log(np.pi)
    # truncated pole kernel on [0,2a], C^2 taper to 2a+1
    T=2*a+1.0; ts=np.linspace(0,T,4000); dt=ts[1]-ts[0]; w=np.full_like(ts,dt); w[0]=w[-1]=dt/2
    k=-4*(np.exp(ts/2)+np.exp(-ts/2)-2); tap=np.ones_like(ts); m=ts>2*a; u=(ts[m]-2*a); tap[m]=1-u**3*(10-15*u+6*u**2)
    k=k*tap; C=np.cos(np.outer(xis,ts)); pole=2*(xis**2)*(C@(w*k))
    # prime 2 (only prime with log n <= 2a=0.8 < log3): -2 (log2/sqrt2) cos(xi log2)
    prime=-2*(np.log(2)/np.sqrt(2))*np.cos(xis*np.log(2))
    return Om+pole+prime
a=0.40; xis=np.linspace(0.05,80,8000); Sv=build_S(a,xis)
mn=Sv.min(); imn=np.argmin(Sv); xistar=xis[imn]
# width of the dip: where S < mn/2 (half-depth), and where S first turns >=0
below=Sv<0
neg_xis=xis[below]
half=Sv< mn/2
halfwidth=(xis[half].max()-xis[half].min()) if half.any() else 0
xi_pos = xis[~below][0] if (~below).any() else None
print(f"a={a}: symbol S_a(xi) bounded (no poles). min S = {mn:.3f} at xi*={xistar:.2f}")
print(f"  negative on xi in [{neg_xis.min():.2f}, {neg_xis.max():.2f}] (width {neg_xis.max()-neg_xis.min():.1f})")
print(f"  dip half-depth width (S<minS/2) = {halfwidth:.2f}")
print(f"  S_a>=0 for xi >= {xis[Sv>=0][ -1 if False else 0]:.1f}? first nonneg at xi={xis[Sv>=0][0]:.2f}; stays>=0 beyond xi={ (xis[np.where(np.diff((Sv>=0).astype(int))==1)][-1] if (Sv>=0).any() else 0):.1f}")
# time-bandwidth for the NEGATIVE band: c = a * (band half-width). If band up to Xi0, c=a*Xi0.
Xi0=neg_xis.max(); c=a*Xi0
print(f"  bad band up to Xi0={Xi0:.1f} => time-bandwidth c=a*Xi0={c:.1f} (large => v_hat CAN concentrate in band => concentration bound FAILS)")
print(f"  BUT dip is {'NARROW' if halfwidth<3 else 'BROAD'} (half-depth width {halfwidth:.1f}) => {'narrow-band concentration bound may work' if halfwidth<3 else 'broad dip, finite-support cannot avoid it -> hard'}")
