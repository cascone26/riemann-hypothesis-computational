# Connes program crux (LEDGER §94/95), stage 1-2: reconstruct k_lambda = E(h_lambda) from the prolate combo and
# verify Connes' PROVEN convergence k_lambda^ -> Xi (Fact hermfact1). Validates the k_lambda machinery before the
# open theta_x comparison (stage 3). Convention-safe: each stage is a unit test against a Connes formula.
#
# Connes eq (1080): h(u) = (pi/2) u^2 (2 pi u^2 - 3) e^{-pi u^2}       [Hermite combo, limit object]
# Connes eq (emap): E(f)(u) = u^{1/2} sum_{n>=1} f(n u)
# Connes eq (1068): k(u) = E(h)(u) = u^{1/2}(pi/2) sum_n n^2 u^2 (2 pi n^2 u^2 - 3) e^{-pi n^2 u^2}   [UNIT TEST 1]
# Connes eq (1072): Xi(t) = int_0^inf k(u) u^{it} du/u   (K(s)=k(e^s) even in s by Poisson k(u)=k(1/u))  [TEST 2]
# k_lambda: replace h by h_lambda = combo of prolate h_{0,lambda},h_{4,lambda} with vanishing integral, sum n<=lambda^2.
import numpy as np

# ---------- prolate PSWF (from prolate_operator.py, Legendre-Galerkin) ----------
def y2_matrix(Nleg):
    Y=np.zeros((Nleg+2,Nleg+2))
    for n in range(Nleg+1):
        if n+1<=Nleg+1: Y[n+1,n]=(n+1)/(2*n+1)*np.sqrt((2*n+1)/(2*n+3))
        if n-1>=0:      Y[n-1,n]=n/(2*n+1)*np.sqrt((2*n+1)/(2*n-1))
    Y=Y[:Nleg,:Nleg]; return Y@Y
def prolate(lmbda,Nleg):
    c=2*np.pi*lmbda**2; n=np.arange(Nleg)
    A=np.diag((n*(n+1)).astype(float))+c**2*y2_matrix(Nleg); A=(A+A.T)/2
    w,V=np.linalg.eigh(A); return V,n
from numpy.polynomial.legendre import legval
def pswf_eval(V_col, y, nmax):
    # h_{.,lambda}(x=lambda y) = sum_m V[m]*P~_m(y),  P~_m = sqrt((2m+1)/2) P_m
    coeff=np.array([V_col[m]*np.sqrt((2*m+1)/2) for m in range(nmax)])
    return legval(y, coeff)

# ---------- Hermite limit object h, k (unit tests) ----------
def h_lim(u): return (np.pi/2)*u**2*(2*np.pi*u**2-3)*np.exp(-np.pi*u**2)
def k_lim(u, nmax=200):
    u=np.atleast_1d(u); s=np.zeros_like(u,float)
    for n in range(1,nmax): s+=h_lim(n*u)
    return np.sqrt(u)*s

# UNIT TEST 1: E(h) == k  (Connes 1068). Direct: k(u)=u^{1/2}(pi/2)sum n^2u^2(2pi n^2u^2-3)e^{-pi n^2u^2}
def k_direct(u,nmax=200):
    u=np.atleast_1d(u); s=np.zeros_like(u,float)
    for n in range(1,nmax): s+=(np.pi/2)*(n*u)**2*(2*np.pi*(n*u)**2-3)*np.exp(-np.pi*(n*u)**2)
    return np.sqrt(u)*s
uu=np.linspace(0.05,3,50)
err1=np.max(np.abs(k_lim(uu)-k_direct(uu)))
print(f"UNIT TEST 1 (E(h)=k, Connes 1068): max|E(h)-k| = {err1:.2e}  {'PASS' if err1<1e-10 else 'FAIL'}")

# TEST 2: Xi(t) = int k(u) u^{it} du/u has first zero at t=14.13 (=gamma_1)
def Xi_from_k(t, kfun, smax=6.0, Ns=4000):
    s=np.linspace(-smax,smax,Ns); ds=s[1]-s[0]; K=kfun(np.exp(s))
    return np.array([np.sum(K*np.cos(tt*s))*ds for tt in np.atleast_1d(t)])
tg=np.linspace(1,40,400); F=Xi_from_k(tg,k_lim)
zeros=[tg[i]-F[i]*(tg[i+1]-tg[i])/(F[i+1]-F[i]) for i in range(len(tg)-1) if F[i]*F[i+1]<0]
print(f"TEST 2 (zeros of Xi from k_lim): first 3 = {[round(z,3) for z in zeros[:3]]}  (true gamma=14.135,21.022,25.011)")

# ---------- k_lambda from prolate combo, convergence to Xi ----------
print("\nSTAGE 2: k_lambda = E(h_lambda), h_lambda = vanishing-integral combo of prolate h_{0,l},h_{4,l}. Zeros of Xi_lambda:")
print(f"{'lambda':>7} {'lambda^2':>8} {'first 3 zeros of Xi_lambda':>34}")
for lmbda in [2.0, 3.0, 4.0, 5.0]:
    Nleg=max(200,int(6*2*np.pi*lmbda**2)); V,n=prolate(lmbda,Nleg)
    # even eigenfunctions: parity by dominant Legendre component even/odd
    par=np.array([(-1)**int(np.argmax(np.abs(V[:,k]))) for k in range(V.shape[1])])
    ev_idx=np.where(par>0)[0]            # ascending-eigenvalue even eigenfns
    i0,i4=ev_idx[0],ev_idx[2]           # h_{0,l}=even#0, h_{4,l}=even#2
    # vanishing integral: int h_{n,l} propto V[0,n]  => combo b*h4 + a*h0 with b*V[0,i4]+a*V[0,i0]=0
    a_coef=1.0; b_coef=-V[0,i0]/V[0,i4] if abs(V[0,i4])>1e-30 else 0.0
    def h_lam(x):   # supported [-lambda,lambda]
        y=x/lmbda; out=np.zeros_like(x,float); m=(np.abs(y)<=1)
        out[m]=a_coef*pswf_eval(V[:,i0],y[m],Nleg)+b_coef*pswf_eval(V[:,i4],y[m],Nleg)
        return out
    def k_lam(u):
        u=np.atleast_1d(u); s=np.zeros_like(u,float)
        nmax=int(np.floor(lmbda**2))+1
        for nn in range(1,nmax): s+=h_lam(nn*u)
        return np.sqrt(u)*s
    a=np.log(lmbda)
    s=np.linspace(-a,a,4000); ds=s[1]-s[0]; K=k_lam(np.exp(s))
    Xl=np.array([np.sum(K*np.cos(tt*s))*ds for tt in tg])
    zl=[tg[i]-Xl[i]*(tg[i+1]-tg[i])/(Xl[i+1]-Xl[i]) for i in range(len(tg)-1) if Xl[i]*Xl[i+1]<0]
    print(f"{lmbda:>7.2f} {lmbda**2:>8.1f} {str([round(z,3) for z in zl[:3]]):>34}", flush=True)
print("\nCONNES Fact hermfact1: Xi_lambda -> Xi as lambda->inf. If first zeros march toward 14.135,21.022,25.011,")
print("the k_lambda machinery is VALIDATED (reproduces Connes' proven convergence) -> ready for the theta_x comparison.")
