# Phase 3 move 1: at a in the first prime step (a_2,a_3)=(0.347,0.549), decompose Q(v_a) on the FULL
# minimizer into arch-core / arch-tail / POLE / PRIME-2 parts. Reveals the actual sign structure to prove.
import numpy as np, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
import suzuki_screw_lambda_a as S

def build_Q(a, gfunc, M, K):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=S.vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gv=np.array([gfunc(d,Lam,Nmax) for d in diffs])
    idx=np.arange(M+1); G=gv[(idx[:,None]-idx[None,:])+M]
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h; Wp=Vp*w[:,None]
    Q=Wp.T@G@Wp; return (Q+Q.T)/2

# component g-functions (all take (t,Lam,Nmax))
def g_pole(t,Lam,Nmax):
    at=abs(t); return -4.0*(np.exp(t/2)+np.exp(-t/2)-2.0)
def g_archcore(t,Lam,Nmax):
    at=abs(t); return -(at/2.0)*S.PSI_QUARTER_MINUS_LOGPI
def g_archtail(t,Lam,Nmax):
    at=abs(t); return -0.25*(S.PHI_1_2_QUARTER - np.exp(-at/2)*S.phi_tail(t))
def g_prime(t,Lam,Nmax):
    return S.screw_g(t,Lam,Nmax)-S.screw_g_no_primes(t,Lam,Nmax)

for a in [0.40,0.45,0.50]:
    M=max(6000,int(500/a)); K=70
    Qf=build_Q(a,S.screw_g,M,K)
    ev,V=np.linalg.eigh(Qf); v=V[:,0]           # full minimizer (coeffs), lam_a=ev[0]/a
    lam=ev[0]/a
    def comp(gf): 
        Q=build_Q(a,gf,M,K); return float(v@Q@v)/a
    cp=comp(g_pole); cc=comp(g_archcore); ct=comp(g_archtail); cpr=comp(g_prime)
    tot=cp+cc+ct+cpr
    print(f"a={a:.2f}  lambda_a={lam:+.5f}  = POLE {cp:+.4f} + ARCHcore {cc:+.4f} + ARCHtail {ct:+.4f} + PRIME2 {cpr:+.5f}  (sum {tot:+.5f})", flush=True)
print("\nReveals which parts are +/- on the minimizer => the inequality to prove for the one-prime step.")
