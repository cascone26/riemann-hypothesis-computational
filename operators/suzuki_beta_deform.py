#!/usr/bin/env python3
"""suzuki_beta_deform.py -- §248: beta-deformation control. Scale prime weights Lam(q)/sqrt(q) -> beta*...
beta=1: true system, zeros real (RH side), Phi'(peak) = 2pi^2 a^2 growing.
beta!=1: deformed EF -- 'zeros' generically OFF-LINE (Beurling-type) -> prediction: Phi'-peaks FLOOR
(width ~ zeta-width of deformed zero) instead of growing a^2. beta->1: recover 2pi^2.
This calibrates the RH-discriminator on a KNOWN off-line family."""
import numpy as np, sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import suzuki_exact_Q as SQ
import suzuki_phase_density as PD
from scipy.linalg import solve

BETA = 1.0
_orig_build = SQ.build_Q
def build_Q_beta(a, K, NG=320, return_pieces=False):
    Q, pieces = _orig_build(a, K, NG, return_pieces=True)
    Qb = pieces["pole"] + pieces["core"] + pieces["tail"] + BETA*pieces["prime"]
    Qb = (Qb + Qb.T)/2
    if return_pieces: return Qb, pieces
    return Qb

def peaks_at(a, zwin=(11.5, 34.0), nz=4600):
    """Return top local maxima of Phi'(z) in window (peak height, position)."""
    K = int(np.ceil(50*2*a/np.pi)) + 40
    kk = np.arange(1, K+1)
    Q = build_Q_beta(a, K)
    lam_a = np.linalg.eigvalsh(Q)[0]/a
    Ta = Q - (lam_a-0.2)*a*np.eye(K)
    fp = SQ.E_p(1.0, a, kk); fm = SQ.E_p(-1.0, a, kk)
    cp = solve(Ta, fp, assume_a='sym'); cm = solve(Ta, fm, assume_a='sym')
    cm *= np.sqrt(abs(cp@fp))/np.sqrt(abs(cm@fm))
    om = kk*np.pi/(2*a)
    zg = np.linspace(zwin[0], zwin[1], nz)
    den = om[:,None]**2 - zg[None,:]**2
    den = np.where(np.abs(den)<1e-9, 1e-9, den)
    SH = np.exp(-1j*zg[None,:]*a)*om[:,None]*(1.0-((-1.0)**kk)[:,None]*np.exp(2j*zg[None,:]*a))/den
    u = (zg-1j)*(cp@SH); w = (zg+1j)*(cm@SH)
    Phi = np.unwrap(np.angle(-u/w))
    dP = np.gradient(Phi, zg)
    pk = [(float(dP[i]), float(zg[i])) for i in range(2, nz-2)
          if dP[i]>dP[i-1] and dP[i]>dP[i+1] and dP[i] > 3*np.median(dP)]
    pk.sort(reverse=True)
    return pk[:8], float(np.median(dP))

if __name__ == "__main__":
    GAM=[14.13,21.02,25.01,30.42,32.94]
    print("beta-deformation: peak Phi' heights vs a  (true zeros gamma in window: 14.1,21.0,25.0,30.4,32.9)")
    out={}
    for beta in [1.0, 0.9, 0.7, 0.5, 1.2]:
        BETA=beta; globals()['BETA']=beta
        print(f"\n===== beta={beta} =====")
        for a in [2.0, 3.0, 4.0]:
            pk, med = peaks_at(a)
            tops=[f"{h:.0f}@{z:.2f}" for h,z in pk[:6]]
            Ks=[h/a**2 for h,z in pk[:5]]
            print(f"  a={a}: median Phi'={med:.1f} (2a={2*a})  peaks: {tops}")
            print(f"        K=peak/a^2: {[f'{k:.1f}' for k in Ks]}")
            out[f"{beta}_{a}"]={"peaks":pk,"med":med}
    json.dump(out, open("../results/beta_deform.json","w"), indent=1)
    print("\nREAD: beta=1 -> peaks at gammas, K~19.7 all a. beta!=1 -> if peaks FLOOR (K falls with a),")
    print("off-line zeros detected by the Phi'-deficit = the discriminator WORKS; if peaks still grow a^2,")
    print("the deformed zeros stayed real (Beurling rigidity) -> different, also informative.")
