"""
lambda_a via sine-Galerkin (same as suzuki_screw_lambda_a) but with the STABLE bounded zero-sum
screw kernel g_zero(t) = -2*S + 2*sum cos(gamma t)/gamma^2  (S=sum 1/gamma^2), which equals the
true screw function g under RH (verified cycle-78 to 4 digits) with NO catastrophic cancellation.
This unblocks the LARGE-a lambda_a margin map (old float A+B+C code lost ~5 digits at t=2a large).
Status: RH-true (zeros) side => a MARGIN/consistency map (does lambda_a stay >0, with what margin),
NOT a proof. Cross-checked vs Thm 1.4 small-a asymptotic lambda_a = log(1/a)+mu1-log2pi+psi(2)-1.
"""
import numpy as np, mpmath as mp, time
G = np.load('zeros2000.npy').astype(float); W = 1.0/(G*G); S = W.sum()

def g_zero_arr(tarr):
    # vectorized bounded screw kernel, chunked
    out = np.empty_like(tarr)
    B=4000
    for i in range(0,len(tarr),B):
        c = np.cos(np.outer(tarr[i:i+B], G))
        out[i:i+B] = -2*S + 2*(c*W).sum(1)
    return out

def lambda_a(a, M=3000, K=70):
    x = np.linspace(-a, a, M+1); h = x[1]-x[0]
    diffs = np.arange(-M, M+1)*h
    gvals = g_zero_arr(diffs)
    idx = np.arange(M+1)
    Gm = gvals[(idx[:,None]-idx[None,:])+M]
    kk = np.arange(1,K+1)
    arg = np.pi*np.outer((x+a),kk)/(2*a)
    Vp = (kk*np.pi/(2*a))*np.cos(arg)
    w = np.full(M+1,h); w[0]=w[-1]=0.5*h
    Wp = Vp*w[:,None]
    Q = Wp.T @ Gm @ Wp; Q=(Q+Q.T)/2
    return (np.linalg.eigvalsh(Q)/a)[0]

if __name__=="__main__":
    const = float(-mp.log(2*mp.pi)+mp.digamma(2)-1)   # -2.4151
    print("STABLE lambda_a(a) via bounded zero-sum g (RH-true side); margin map into large a.")
    print(f"{'a':>7} {'lambda_a':>13} {'log(1/a)+mu1+c (Thm1.4, mu1=0.25)':>34} {'>0?':>5}")
    t0=time.time()
    for a in [0.05,0.1,0.2,0.35,0.5,0.75,1.0,1.5,2.0,3.0,4.0,6.0,8.0]:
        M = min(8000, max(3000, int(500/a)))
        la = lambda_a(a, M=M, K=70)
        thm = np.log(1/a)+0.25+const if a<0.35 else float('nan')
        print(f"{a:7.3f} {la:13.6f} {thm:34.5f} {'OK' if la>0 else 'NEG!':>5}", flush=True)
    print(f"done {time.time()-t0:.1f}s")
