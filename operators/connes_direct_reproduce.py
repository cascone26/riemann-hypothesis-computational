"""
Reproduce Connes' claim in our hands (Handle): the min eigenvector of the DIRECT Weil form G_a (kernel = the
validated screw kernel g, arch+pole-primes) on [-a,a] should have Fourier-transform zeros on the REAL axis
(critical line) approximating the zeta zeros. Uses primes up to e^{2a}.
"""
import numpy as np, mpmath as mp
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def build_G(a, M):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gv=np.array([screw_g(d,Lam,Nmax) for d in diffs])
    idx=np.arange(M+1); G=gv[(idx[:,None]-idx[None,:])+M]
    w=np.sqrt(np.full(M+1,h)); w[0]=w[-1]=np.sqrt(0.5*h)   # symmetric quadrature weights
    Gs=(w[:,None]*G*w[None,:]); Gs=(Gs+Gs.T)/2
    return x, h, Gs, w

def eta_hat(t, x, eta, h):
    # Fourier transform eta_hat(t)=int eta(x) e^{i t x} dx (real part for even eta)
    return np.sum(eta*np.cos(t*x))*h

if __name__=="__main__":
    a=2.565  # ~ log 13 : primes up to e^{5.13}~169 enter (>= Connes' <=13; test mechanism)
    M=600
    x,h,Gs,w=build_G(a,M)
    ev,V=np.linalg.eigh(Gs)
    eta = V[:,0]/w    # undo quadrature weighting -> eigenvector of the form
    eta=eta/np.sqrt(np.sum(eta*eta)*h)
    simple = (ev[1]-ev[0])>1e-6*max(1,abs(ev).max())
    evenerr=np.max(np.abs(eta-eta[::-1]))
    print(f"a={a}, M={M}: min eig={ev[0]:.4e}, simple={simple}, even-err={evenerr:.2e}")
    # scan eta_hat(t) for zeros in t; compare to first zeta zeros
    ts=np.linspace(0.5,40,4000); vals=np.array([eta_hat(t,x,eta,h) for t in ts])
    sign=np.sign(vals); crossings=ts[:-1][np.diff(sign)!=0]
    zz=[float(mp.im(mp.zetazero(k))) for k in range(1,9)]
    print(f"eta_hat(t) real-axis zero crossings (first 8): {np.round(crossings[:8],4)}")
    print(f"first 8 zeta zeros:                            {np.round(zz,4)}")
    if len(crossings)>=1:
        # match nearest
        for k,g in enumerate(zz[:6]):
            near=crossings[np.argmin(np.abs(crossings-g))] if len(crossings) else float('nan')
            print(f"   zeta_{k+1}={g:.4f}  nearest eta_hat zero={near:.4f}  diff={abs(near-g):.4f}")
    print("\nIf eta_hat zeros track zeta zeros: Connes' mechanism reproduced from OUR Weil kernel (Handle).")
    print("If not: the direct screw-g form differs from Connes' exact Q (honest result); note the difference.")
