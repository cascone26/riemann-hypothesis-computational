import numpy as np, os
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "suzuki_screw_lambda_a.py")).read()
cut = src.find("if __name__")
exec(src[:cut] if cut > 0 else src, globals())
def Qmat(a, M, K):
    Nmax = int(np.floor(np.exp(2*a))) + 1
    Lam = vonmangoldt_upto(min(Nmax, 30_000_000))
    x = np.linspace(-a, a, M+1); h = x[1]-x[0]
    diffs = np.arange(-M, M+1)*h
    gvals = np.array([screw_g(d, Lam, Nmax) for d in diffs])
    idx = np.arange(M+1)
    G = gvals[(idx[:,None]-idx[None,:])+M]
    kk = np.arange(1, K+1)
    arg = np.pi*np.outer((x+a), kk)/(2*a)
    Vp = (kk*np.pi/(2*a))*np.cos(arg)
    w = np.full(M+1, h); w[0]=w[-1]=0.5*h
    Wp = Vp*w[:,None]
    Q = Wp.T@G@Wp
    return (Q+Q.T)/2/a, x, kk
for a in (3.0, 4.0, 5.0, 6.0):
    K = int(20*a); M = 4000
    Q, x, kk = Qmat(a, M, K)
    ev, V = np.linalg.eigh(Q)
    # convergence spot-check at K+10
    Q2, _, _ = Qmat(a, M, K+10)
    ev2 = np.linalg.eigvalsh(Q2)
    # near-null eigenvector width (x-domain): reconstruct, measure central-lobe scale
    B = np.sin(np.pi*np.outer((x+a), kk)/(2*a))
    u0 = B@V[:,0]
    p = u0**2/np.trapezoid(u0**2, x)
    xbar = np.trapezoid(x*p, x); width = np.sqrt(np.trapezoid((x-xbar)**2*p, x))
    print(f"a={a:.1f}: min/max {ev[0]/ev[-1]:+.3e} (K+10: {ev2[0]/ev2[-1]:+.3e}) | null-vec rms width {width:.3f}", flush=True)
print("DONE")
