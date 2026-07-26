#!/usr/bin/env python3
"""suzuki_parking_force.py -- §244: Hellmann-Feynman decomposition of the parking force.

At a parked branch z_j(a) (W(a,theta;z_j)=0, Thm 1.5 zeros real): dz_j/da = -Re[dW/da / dW/dz].
dW/da splits: kinematic (explicit a in shat/f/mass/domain) + per-piece dQ/da contributions
(pole, core, tail, each prime power n). Pieces are EXACT closed forms (suzuki_exact_Q), so
dP_X/da by central difference h=1e-3 is clean; contribution of X = directional derivative of W
under Q -> Q + eps*dP_X/da at fixed a. Kinematic = total (full central diff in a) - sum(pieces).
READ: which primes hold the branch at gamma_j -- the resonance ledger of §241/§243's c(a)*gamma/a law.
"""
import numpy as np, sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import suzuki_exact_Q as SQ
from scipy.linalg import solve
from scipy.interpolate import CubicSpline

# splined Lerch kernel (same speedup as fast scan)
_ts = np.logspace(-9, np.log10(8.5), 3000)
_ks = np.array([0.5*np.exp(-t/2)*SQ.lerch1_quarter_t(t) for t in _ts])
_spl = CubicSpline(np.log(_ts), _ks)
def tail_nodes_fast(a, NG=320):
    L = 2*a
    xs, ws = np.polynomial.legendre.leggauss(NG)
    s = 0.5*(xs+1.0); t = L*s**3
    return t, 0.5*ws*3*L*s**2*_spl(np.log(np.maximum(t,1e-9)))
SQ.tail_nodes = tail_nodes_fast

GAMMAS = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062]

def W_of_Q(Qt, a, K, theta, z):
    """W at (theta,z) given an arbitrary symmetric Q-matrix (fixed a elsewhere)."""
    kk = np.arange(1, K+1)
    lam_a = np.linalg.eigvalsh(Qt)[0]/a
    Ta = Qt - (lam_a-0.2)*a*np.eye(K)
    fp = SQ.E_p(1.0, a, kk); fm = SQ.E_p(-1.0, a, kk)
    cp = solve(Ta, fp, assume_a='sym'); cm = solve(Ta, fm, assume_a='sym')
    cm *= np.sqrt(abs(cp@fp))/np.sqrt(abs(cm@fm))
    om = kk*np.pi/(2*a)
    den = om*om - z*z; den = np.where(np.abs(den)<1e-9, 1e-9, den)
    sh = np.exp(-1j*z*a)*om*(1.0-((-1.0)**kk)*np.exp(2j*z*a))/den
    return (z-1j)*(cp@sh) + np.exp(1j*theta)*(z+1j)*(cm@sh)

def find_root(a, K, theta, guess):
    Q = SQ.build_Q(a, K, NG=320)
    lo, hi = guess-0.8, guess+0.8
    zs = np.linspace(lo, hi, 161)
    Wz = np.abs([W_of_Q(Q, a, K, theta, z) for z in zs])
    i = int(np.argmin(Wz)); lo, hi = zs[max(0,i-1)], zs[min(160,i+1)]
    for _ in range(50):
        m1, m2 = lo+0.382*(hi-lo), lo+0.618*(hi-lo)
        if abs(W_of_Q(Q, a, K, theta, m1)) < abs(W_of_Q(Q, a, K, theta, m2)): hi = m2
        else: lo = m1
    return 0.5*(lo+hi), Q

def decompose(a, theta=0.0, h=1e-3, eps=1e-4):
    K = int(np.ceil(50*2*a/np.pi)) + 40
    _, piecesA = SQ.build_Q(a-h, K, NG=320, return_pieces=True)
    _, piecesB = SQ.build_Q(a+h, K, NG=320, return_pieces=True)
    dP = {}
    for key in ("pole", "core", "tail"):
        dP[key] = (piecesB[key]-piecesA[key])/(2*h)
    for n in piecesB["per_prime"]:
        Pa = piecesA["per_prime"].get(n, np.zeros_like(piecesB["per_prime"][n]))
        dP[f"n={n}"] = (piecesB["per_prime"][n]-Pa)/(2*h)
    rows = []
    for j, gm in enumerate(GAMMAS):
        z0, Q = find_root(a, K, theta, gm)
        dz = 1e-4
        dWdz = (W_of_Q(Q, a, K, theta, z0+dz) - W_of_Q(Q, a, K, theta, z0-dz))/(2*dz)
        # total velocity: full central diff in a (rebuild everything)
        Kh = K
        Qa = SQ.build_Q(a-h, Kh, NG=320); Qb = SQ.build_Q(a+h, Kh, NG=320)
        dWda_tot = (W_of_Q(Qb, a+h, Kh, theta, z0) - W_of_Q(Qa, a-h, Kh, theta, z0))/(2*h)
        v_tot = float(np.real(-dWda_tot/dWdz))
        contrib = {}
        for key, D in dP.items():
            dW = (W_of_Q(Q+eps*D, a, K, theta, z0) - W_of_Q(Q-eps*D, a, K, theta, z0))/(2*eps)
            contrib[key] = float(np.real(-dW/dWdz))
        v_pieces = sum(contrib.values())
        contrib["KINEMATIC(remainder)"] = v_tot - v_pieces
        rows.append({"gamma": gm, "z0": float(z0), "v_total": v_tot,
                     "free_scale(-z/a)": float(-z0/a), "contrib": contrib})
    return rows

if __name__ == "__main__":
    for a in [1.5, 2.2]:
        print(f"\n================ a={a}  (theta=0) ================")
        rows = decompose(a)
        for r in rows:
            print(f"\n gamma={r['gamma']:.3f}  z0={r['z0']:.4f}  v_total={r['v_total']:+.4f}  "
                  f"(free-scale would be {r['free_scale(-z/a)']:+.2f})")
            c = r["contrib"]
            groups = {}
            for k, v in c.items():
                if k.startswith("n="):
                    n = int(k[2:]); p = min([q for q in (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79) if n % q == 0 and all(n % qq for qq in (2,3,5,7,11,13) if qq < q and n % qq == 0)] or [n])
                    # base prime of prime power
                    b = 2
                    while b*b <= n:
                        if n % b == 0: break
                        b += 1
                    base = b if n % b == 0 and n != b or True else n
                    bb = None
                    for q in range(2, n+1):
                        if n % q == 0: bb = q; break
                    groups.setdefault(f"p={bb}", 0.0); groups[f"p={bb}"] += v
                else:
                    groups[k] = v
            big = sorted(((k, v) for k, v in groups.items()), key=lambda kv: -abs(kv[1]))
            print("   " + "  ".join(f"{k}:{v:+.3f}" for k, v in big if abs(v) > 5e-4))
        json.dump(rows, open(f"../results/parking_force_a{a}.json", "w"), indent=1)
    print("\n-> results/parking_force_a*.json")
