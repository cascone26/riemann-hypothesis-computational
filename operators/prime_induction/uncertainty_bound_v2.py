# v2 CORRECTED symbol assembly. Exact multiplier identity for supp v in (-a,a):
#   2pi*Q(v) = int S_a(xi) |vhat(xi)|^2 dxi,   S_a = OMEGA(xi) [arch, exact classical]
#            + POLE_a(xi) [smooth-truncated pole kernel, xi^2-FT, C^2 taper on [2a,2a+1]]
#            + PRIME(xi)  [= -2 sum_{log n<=2a} Lam(n)/sqrt(n) cos(xi log n), exact §81]
# VALIDATE S_a against direct lambda_a on the minimizer, then bad set + Donoho-Stark bound.
import numpy as np, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
import suzuki_screw_lambda_a as S
from scipy.special import digamma
def Omega(xi): return np.real(digamma(0.25+0.5j*xi)) - np.log(np.pi)
def pole_mult(a, xis):
    # kernel k(t)=-4(e^{t/2}+e^{-t/2}-2) on |t|<=2a, tapered C^2 to 0 on [2a,2a+1]; multiplier = xi^2 * khat
    T=2*a+1.0; ts=np.linspace(0,T,6000); dt=ts[1]-ts[0]; w=np.full_like(ts,dt); w[0]=w[-1]=dt/2
    k=-4*(np.exp(ts/2)+np.exp(-ts/2)-2)
    taper=np.ones_like(ts); m=ts>2*a; u=(ts[m]-2*a)/1.0
    taper[m]=1-u**3*(10-15*u+6*u**2)   # C^2 smoothstep down
    k=k*taper
    out=np.empty_like(xis)
    for i,xi in enumerate(xis):
        out[i]=2*(xi**2)*np.sum(w*k*np.cos(xi*ts))
    return out
def prime_mult(a, xis):
    out=np.zeros_like(xis); Nmax=int(np.floor(np.exp(2*a)))
    if Nmax>=2:
        Lam=S.vonmangoldt_upto(max(Nmax,3))
        for n in range(2,Nmax+1):
            if Lam[n]>0: out += -2*(Lam[n]/np.sqrt(n))*np.cos(xis*np.log(n))
    return out
def build_S(a, xis): return Omega(xis)+pole_mult(a,xis)+prime_mult(a,xis)
# --- VALIDATION: S_a-quotient on the minimizer vs direct lambda_a ---
def validate(a):
    M=max(7000,int(550/a)); K=80
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=S.vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; diffs=np.arange(-M,M+1)*h
    gv=np.array([S.screw_g(d,Lam,Nmax) for d in diffs]); idx=np.arange(M+1); G=gv[(idx[:,None]-idx[None,:])+M]
    kk=np.arange(1,K+1); Vp=(kk*np.pi/(2*a))*np.cos(np.pi*np.outer((x+a),kk)/(2*a))
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h; Wp=Vp*w[:,None]; Q=Wp.T@G@Wp; Q=(Q+Q.T)/2
    ev,V=np.linalg.eigh(Q); c=V[:,0]; lam_direct=ev[0]/a
    vx=np.sin(np.pi*np.outer((x+a),kk)/(2*a))@c
    xis=np.linspace(-400,400,32000); dxi=xis[1]-xis[0]
    vhat=np.array([np.sum(vx*np.exp(-1j*xi*x)*w) for xi in xis])
    Sv=build_S(a,xis)
    quot=np.sum(Sv*np.abs(vhat)**2)*dxi/(2*np.pi)/np.sum(vx*vx*w)
    print(f"  validate a={a}: direct lambda_a={lam_direct:+.6f}  symbol-quotient={quot:+.6f}  (match => S_a correct)")
    return abs(quot-lam_direct)<5e-3
ok=validate(0.40)
# --- bad set + DS bound ---
for a in [0.30,0.37,0.42,0.50]:
    xis=np.linspace(0,400,16000); Sv=build_S(a,xis)
    bad=Sv<0
    Bmeas=2*(xis[1]-xis[0])*bad.sum()
    if bad.any():
        infB=Sv[bad].min(); ginf=Sv[~bad].min(); bmax=xis[bad].max()
        lower=(ginf*2*np.pi+(infB-ginf)*2*a*Bmeas)/(2*np.pi)
        print(f"a={a:.2f}: |B|={Bmeas:6.2f} (bad xi up to {bmax:6.1f}), inf_B={infB:8.2f}, inf_good={ginf:6.3f} => DS bound lambda_a >= {lower:+9.3f}")
    else:
        print(f"a={a:.2f}: S_a>=0 everywhere")
