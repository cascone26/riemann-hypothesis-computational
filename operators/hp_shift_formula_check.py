"""S2875: numeric check of the S2874 closed form for the first-order node shift (lam2=14, N=8, K=8)."""
import numpy as np, mpmath as mp
mp.mp.dps=40
lam2=14; N=8; L=mp.log(lam2)
nu=[2*mp.pi*n/L for n in range(-N,N+1)]
G=[mp.zetazero(k).imag for k in range(1,N+1)]
P=lambda t: mp.fprod(t*t-g*g for g in G)
Q=lambda t: mp.fprod(t+v for v in nu)
w=lambda t: 4*mp.sin(t*L/2)**2/L
# tail zeros (positive), 100k list beyond gamma_N
zs=np.loadtxt('data/zeros_100k.txt',usecols=-1)[N:]
Lf=float(L); Gf=np.array([float(g) for g in G]); nuf=np.array([float(v) for v in nu])
def Wf(t):  # float weights w P^2/Q^2 for array t (log-space for stability)
    lp=np.sum(np.log(np.abs(t[:,None]**2-Gf[None,:]**2)),1); lq=np.sum(np.log(np.abs(t[:,None]+nuf[None,:])),1)
    return 4*np.sin(t*Lf/2)**2/Lf*np.exp(2*lp-2*lq)
Wz=Wf(zs)
pred=[]
for k in range(1,6):
    gk=G[k-1]; gkf=float(gk)
    qk=Q(gk)/mp.diff(P,gk)
    gz=1/(zs**2-gkf**2)
    SW=Wz.sum(); SWg=(Wz*gz).sum()
    # density tail beyond last zero: W(t)~(2/L)*t^-2*..., g~t^-2 -> negligible; add avg-sin^2 estimate for SW
    T=zs[-1]; SW+= float(mp.quad(lambda t:(2/L)*mp.log(t/(2*mp.pi))/(2*mp.pi)*t**(4*N)/t**(4*N+2),[T,mp.inf]))
    f=lambda t: w(t)*P(t)**2/Q(t)**2
    pts=[0]+[g for g in G]+[60,200,mp.inf]
    m0=mp.quad(f,pts); m1=mp.quad(lambda t:f(t)/(t*t-gk*gk),pts)
    gm=m1/m0
    dt=2*qk**2*gk/w(gk)*(SWg-SW*gm)
    pred.append(dt)
    print(f"k={k}: <g>_tail={SWg/SW:.6e}  <g>_m~={mp.nstr(gm,6)}  sign={'+' if SWg/SW>gm else '-'}  predicted dt/s={mp.nstr(dt,5)}")
meas=[1.69e-13,2.05e-4,1.72,4.35e4,1.14e6]
print("ratio measured/predicted:",[mp.nstr(m/p,5) for m,p in zip(meas,pred)])
