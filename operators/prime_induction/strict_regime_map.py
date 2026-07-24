# Map the strict->tight transition in the first prime step: arch+pole MARGIN vs PRIME-2 magnitude, fine a-grid
# just above a_2=0.347. If margin >> prime near a_2 => genuine provable strict sub-regime (modest theorem).
import numpy as np, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
import suzuki_screw_lambda_a as S
def build_Q(a,gfunc,M,K):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=S.vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; diffs=np.arange(-M,M+1)*h
    gv=np.array([gfunc(d,Lam,Nmax) for d in diffs]); idx=np.arange(M+1); G=gv[(idx[:,None]-idx[None,:])+M]
    kk=np.arange(1,K+1); Vp=(kk*np.pi/(2*a))*np.cos(np.pi*np.outer((x+a),kk)/(2*a))
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h; Wp=Vp*w[:,None]; Q=Wp.T@G@Wp; return (Q+Q.T)/2
def g_archpole(t,Lam,Nmax): return S.screw_g_no_primes(t,Lam,Nmax)
def g_prime(t,Lam,Nmax): return S.screw_g(t,Lam,Nmax)-S.screw_g_no_primes(t,Lam,Nmax)
print("a      lambda_a    archpole-margin   |prime2|    ratio margin/|prime|  (>>1 strict, ~1 tight)")
for a in [0.355,0.37,0.39,0.41,0.43,0.45,0.48,0.51]:
    M=max(7000,int(550/a)); K=80
    Qf=build_Q(a,S.screw_g,M,K); ev,V=np.linalg.eigh(Qf); v=V[:,0]; lam=ev[0]/a
    Qap=build_Q(a,g_archpole,M,K); Qp=build_Q(a,g_prime,M,K)
    m=float(v@Qap@v)/a; p=float(v@Qp@v)/a
    ratio=m/abs(p) if abs(p)>1e-12 else float('inf')
    print(f"{a:.3f}  {lam:+.6f}   {m:+.6f}      {abs(p):.6f}    {ratio:8.2f}", flush=True)
print("strict sub-regime (margin>>|prime|) near a_2 => provable; tight (ratio->1) as a->0.5 => the hard onset.")
