# Verify: Q_p2(v) = -(Lambda2/sqrt2)(1/pi) int cos(xi*log2)|vhat|^2 dxi  (prime form = cosine multiplier).
# Compare formula vs direct-kernel Q_p2 on the minimizer v_a. If match => clean rigorous formula confirmed.
import numpy as np, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
import suzuki_screw_lambda_a as S
LOG2=np.log(2); W2=LOG2/np.sqrt(2)
def build_Q(a,gfunc,M,K):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=S.vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; diffs=np.arange(-M,M+1)*h
    gv=np.array([gfunc(d,Lam,Nmax) for d in diffs]); idx=np.arange(M+1); G=gv[(idx[:,None]-idx[None,:])+M]
    kk=np.arange(1,K+1); Vp=(kk*np.pi/(2*a))*np.cos(np.pi*np.outer((x+a),kk)/(2*a))
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h; Wp=Vp*w[:,None]; Q=Wp.T@G@Wp; return (Q+Q.T)/2, x, w, kk
def g_prime(t,Lam,Nmax): return S.screw_g(t,Lam,Nmax)-S.screw_g_no_primes(t,Lam,Nmax)
for a in [0.40,0.45,0.50]:
    M=max(7000,int(550/a)); K=80
    Qf,x,w,kk=build_Q(a,S.screw_g,M,K); ev,V=np.linalg.eigh(Qf); c=V[:,0]
    Qp,_,_,_=build_Q(a,g_prime,M,K)
    direct=float(c@Qp@c)   # direct-kernel Q_p2(v) (un-normalized)
    # v(x)=sum c_k sin(k pi (x+a)/2a); vhat(xi)=int_{-a}^a v e^{-i xi x} dx  -> compute on xi-grid
    vx=np.sin(np.pi*np.outer((x+a),kk)/(2*a))@c
    xis=np.linspace(-120,120,24000); dxi=xis[1]-xis[0]
    vhat=np.array([np.sum(vx*np.exp(-1j*xi*x)*w) for xi in xis])
    formula=-(W2/np.pi)*np.sum(np.cos(xis*LOG2)*np.abs(vhat)**2)*dxi
    nrm=(1/(2*np.pi))*np.sum(np.abs(vhat)**2)*dxi   # ||v||^2 via Parseval (check vs direct)
    nrm_direct=np.sum(vx*vx*w)
    print(f"a={a:.2f}: Q_p2 direct={direct:+.6e}  formula={formula:+.6e}  ratio={direct/formula:.4f}   (||v||^2 parseval={nrm:.4f} vs direct={nrm_direct:.4f})", flush=True)
print("ratio~1 => prime form = cosine multiplier -(L2/sqrt2)(1/pi)cos(xi log2) CONFIRMED (=> full prime = Re[zeta'/zeta]).")
