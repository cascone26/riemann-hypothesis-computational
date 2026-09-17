import numpy as np, sys, os
sys.path.insert(0, "operators")
import suzuki_W_zeros as SW

GAM = np.array([14.1347,21.0220,25.0109,30.4249,32.9351])

def theta_match(a, K, M, zgrid, thetas):
    vhp, vhm, lam_a = SW.v_hat_pm(a, M=M, K=K)
    best=(9e9,None,None)
    for th in thetas:
        zs,med,offr = SW.zeros_and_reality(vhp, vhm, th, zgrid)
        if len(zs)<5: continue
        zs=np.array(sorted(zs))
        # mean nearest-match to first 5 gammas
        err=np.mean([min(abs(zs-g)) for g in GAM])
        if err<best[0]: best=(err,th,zs[:6])
    return lam_a, best

zgrid=np.linspace(0.5,45,5000)
thetas=np.linspace(0,2*np.pi,60,endpoint=False)
print(f"{'a':>4} {'K':>4} {'lam_a':>9} {'theta*':>7} {'meanErr(5γ)':>11}  zeros")
for a in [0.8, 1.2, 1.6, 2.0]:
    K = 80 if a<1.3 else 120
    M = 5000
    lam_a, (err,th,zs) = theta_match(a, K, M, zgrid, thetas)
    zshow = [round(z,2) for z in (zs if zs is not None else [])]
    print(f"{a:>4} {K:>4} {lam_a:>9.4f} {(th if th else float('nan')):>7.3f} {err:>11.4f}  {zshow}", flush=True)
