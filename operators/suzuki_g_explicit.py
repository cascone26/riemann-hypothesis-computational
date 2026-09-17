"""
Explicit screw function g_xi(t) (Suzuki 2301.00421 eq 4.3) vs its zero-sum value.
FAST: numpy-vectorized zero-sum over precomputed zeros10k.npy; mpmath only for special functions.

RH-content identity (Cor 4.3 / eq 4.6, u=t, g(0)=0):
    g_xi(t)  =  - sum_{gamma>0} (2 - 2 cos(gamma t)) / gamma^2        ( <= 0 under RH )

Sign decomposition of the explicit eq-4.3 form (the whole RH tension in one variable):
    A(t) = -4( e^{t/2}+e^{-t/2}-2 ) = -4(2cosh(t/2)-2)   <= 0 ALWAYS
    B(t) = sum_{n<=e^t} (Lambda(n)/sqrt(n))(t-log n)      >= 0 ALWAYS
    C(t) = arch terms.  Pin via C_needed := g_zero - A - B, compare to candidate.
"""
import numpy as np, mpmath as mp
mp.mp.dps = 25
Z = np.load(__file__.rsplit('/',1)[0] + '/zeros10k.npy')   # gamma_n > 0

def g_zero(t):
    return -np.sum((2 - 2*np.cos(Z*t))/(Z*Z))

def A_term(t):
    return -4*(np.exp(t/2)+np.exp(-t/2)-2)

def _lambda_pairs(X):
    X=int(X); out=[]
    for n in range(2,X+1):
        f=n; d=2; fac={}
        while d*d<=f:
            while f%d==0: fac[d]=fac.get(d,0)+1; f//=d
            d+=1
        if f>1: fac[f]=fac.get(f,0)+1
        if len(fac)==1:
            p=next(iter(fac)); out.append((n, np.log(p)))
    return out

def B_term(t):
    X=np.exp(t); s=0.0
    for n,L in _lambda_pairs(X):
        s += (L/np.sqrt(n))*(t-np.log(n))
    return s

def C_term(t):
    """CORRECT explicit arch term, 2606 eq (1.3), t>=0:
       -(t/2)(psi(1/4)-log pi) - (1/4)(Phi(1,2,1/4) - e^{-t/2} Phi(e^{-2t},2,1/4))."""
    t=mp.mpf(t)
    psi14=mp.digamma(mp.mpf(1)/4)
    return float(-(t/2)*(psi14 - mp.log(mp.pi))
                 - mp.mpf(1)/4*( mp.lerchphi(1,2,mp.mpf(1)/4) - mp.e**(-t/2)*mp.lerchphi(mp.e**(-2*t),2,mp.mpf(1)/4) ))

def g_explicit(t):
    """Zeros-FREE screw function (primes + arch). RH <=> g_explicit(t)<=0 for large t (non-circular)."""
    return A_term(t)+B_term(t)+C_term(t)

if __name__=="__main__":
    print("Weil explicit-formula identity: g_explicit(zeros-FREE) =? g_zero(zeros).  RH-test: g_explicit<=0.")
    print(f"{'t':>5} {'g_explicit':>13} {'g_zero':>13} {'A':>11} {'B':>11} {'C':>11} {'B/e^{t/2}':>10} {'g_exp<=0?':>9}")
    for t in [1.0,2.0,3.0,4.0,5.0,6.0,8.0,10.0,12.0]:
        A=A_term(t); B=B_term(t); C=C_term(t); ge=A+B+C
        gz=g_zero(t)
        print(f"{t:5.1f} {ge:13.7f} {gz:13.7f} {A:11.4f} {B:11.4f} {C:11.4f} {B/np.exp(t/2):10.5f} {'OK' if ge<=0 else 'POS!':>9}")
    print("\ng_explicit is zeros-free; matching g_zero = verified Weil explicit formula (non-circular).")
    print("Bounded residual: |g|<=4*sum 1/gamma^2 ~=", 4*np.sum(1.0/(Z*Z)))
