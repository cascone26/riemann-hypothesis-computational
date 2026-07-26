import sys, json, numpy as np, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import suzuki_exact_Q as SQ
from scipy.linalg import solve
from scipy.interpolate import CubicSpline

# global Lerch kernel table: ker(t) = 0.5 e^{-t/2} Phi(e^{-2t},1,1/4), log-t spline
_ts = np.logspace(-9, np.log10(6.5), 3000)
_ks = np.array([0.5*np.exp(-t/2)*SQ.lerch1_quarter_t(t) for t in _ts])
_spl = CubicSpline(np.log(_ts), _ks)
def ker_fast(t): return _spl(np.log(np.maximum(t, 1e-9)))

def tail_nodes_fast(a, NG=320):
    L = 2*a
    xs, ws = np.polynomial.legendre.leggauss(NG)
    s = 0.5*(xs+1.0); wgt = 0.5*ws
    t = L*s**3; jac = 3*L*s**2
    return t, wgt*jac*ker_fast(t)

SQ.tail_nodes = tail_nodes_fast   # monkey-patch: same math, cached kernel

GAMMAS = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178, 40.918719]
ZLO, ZHI, NZ = 10.0, 44.0, 3600
THETAS = [0.0, np.pi/2, np.pi, 1.5*np.pi]

def scan_a(a):
    K = int(np.ceil((ZHI+6)*2*a/np.pi)) + 40
    kk = np.arange(1, K+1)
    Q = SQ.build_Q(a, K, NG=320)
    lam_a = np.linalg.eigvalsh(Q)[0]/a
    Ta = Q - (lam_a-0.2)*a*np.eye(K)
    fp = SQ.E_p(1.0, a, kk); fm = SQ.E_p(-1.0, a, kk)
    cp = solve(Ta, fp, assume_a='sym'); cm = solve(Ta, fm, assume_a='sym')
    cm *= np.sqrt(abs(cp@fp))/np.sqrt(abs(cm@fm))
    zg = np.linspace(ZLO, ZHI, NZ)
    om = kk*np.pi/(2*a)
    den = om[:,None]**2 - zg[None,:]**2
    den = np.where(np.abs(den)<1e-9, 1e-9, den)
    SH = np.exp(-1j*zg[None,:]*a)*om[:,None]*(1.0-((-1.0)**kk)[:,None]*np.exp(2j*zg[None,:]*a))/den  # (K,NZ)
    vp = cp @ SH; vm = cm @ SH                                     # (NZ,)
    out = {}
    for ti, th in enumerate(THETAS):
        Wz = np.abs((zg-1j)*vp + np.exp(1j*th)*(zg+1j)*vm)
        med = np.median(Wz)
        zs = []
        for i in range(1, NZ-1):
            if Wz[i]<Wz[i-1] and Wz[i]<Wz[i+1] and Wz[i]<0.05*med:
                y0,y1,y2 = Wz[i-1],Wz[i],Wz[i+1]
                zs.append(zg[i] + 0.5*(y0-y2)/(y0-2*y1+y2+1e-30)*(zg[1]-zg[0]))
        out[ti] = {str(j): float(min(zs, key=lambda z: abs(z-gm))) for j,gm in enumerate(GAMMAS) if zs}
        out[f"n{ti}"] = len(zs)
    return K, out

if __name__ == "__main__":
    lo, hi = float(sys.argv[1]), float(sys.argv[2]); tag = sys.argv[3]
    avals = [round(x,3) for x in np.arange(lo, hi, 0.01)]
    rows = []
    t0 = time.time()
    for a in avals:
        K, out = scan_a(a)
        rows.append({"a": a, "K": K, "out": {str(k): v for k,v in out.items()}})
        json.dump(rows, open(f"../results/gf2fast_{tag}.json","w"))
    print(f"{tag} DONE {len(rows)} pts in {time.time()-t0:.0f}s", flush=True)
