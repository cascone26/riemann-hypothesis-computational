import numpy as np, sys, os
sys.path.insert(0, os.path.abspath('operators'))
import l2_edge_sweep as m
# HEAT-FLOW REGULARIZES THE RAZOR (Path A x Path B synthesis, per LEDGER Riemann7 §500).
# de Bruijn-Newman zero flow: dgamma_k/dt = sum_{j!=k} 2/(gamma_k - gamma_j) over ALL zeros (+/- symmetric).
# t>0 SPREADS the zeros. The a=1 Weil-razor (min-eig ~ 0 at the ACTUAL arithmetic zeros, tuned basis) should
# gain POSITIVE MARGIN as the flow regularizes the zeros => links our a=1 razor to Rodgers-Tao Lambda>=0.
G0 = np.load('operators/zeros500.npy')   # positive gamma, ascending
G0 = G0[G0 < 200]                        # low zeros (well-resolved), ~ first 80
def flow_deriv(g):
    # full symmetric set {+g, -g}; d g_k/dt from repulsion by all others
    allz = np.concatenate([g, -g])
    d = np.zeros_like(g)
    for k in range(len(g)):
        diff = g[k] - allz
        diff = diff[np.abs(diff) > 1e-9]
        d[k] = np.sum(2.0/diff)
    return d
def evolve(g, t, steps=200):
    g = g.copy(); dt = t/steps
    for _ in range(steps):
        k1=flow_deriv(g); k2=flow_deriv(g+0.5*dt*k1); k3=flow_deriv(g+0.5*dt*k2); k4=flow_deriv(g+dt*k3)
        g = g + (dt/6)*(k1+2*k2+2*k3+k4)
    return g
# probe basis (the a=1 razor basis) on a well-conditioned band
band=(14,64); w=3.0
cs = m.probes(band, w/2.0)
def zeroside_mineig(gammas):
    FZ = m.Fmat(cs, w, gammas)            # [nprobe, nzero]
    Z = 2.0 * FZ @ FZ.T                    # zero-side Weil Gram
    Z = 0.5*(Z+Z.T)
    # normalize by the overlap Gram to get the a=1-style generalized margin
    G = m.build_forms(cs, w)['overlap']; G=0.5*(G+G.T)
    ge,gv=np.linalg.eigh(G); keep=ge>1e-10*ge[-1]
    B=gv[:,keep]/np.sqrt(ge[keep])[None,:]
    return float(np.linalg.eigvalsh(B.T@Z@B)[0])
print("Heat-flow regularizes the a=1 razor: zero-side Weil margin vs flow time t")
print(f"(band {band}, {len(G0)} zeros; t=0 = arithmetic zeros = the razor)")
print(f"{'t':>8} {'min-eig(zero-side)':>20} {'zero-spread change':>18}")
sp0 = np.mean(np.diff(G0[:40]))
for t in [0.0, 0.002, 0.005, 0.01, 0.02, 0.04, 0.08]:
    g = G0 if t==0 else evolve(G0, t)
    me = zeroside_mineig(g)
    sp = np.mean(np.diff(g[:40]))
    print(f"{t:>8.3f} {me:>20.6e} {(sp-sp0)/sp0*100:>17.2f}%", flush=True)
print("\nREAD: if min-eig grows from ~0 (t=0 razor) to positive as t>0 (zeros spread), the de Bruijn-Newman")
print("flow REGULARIZES the a=1 razor => our razor (§470) = Rodgers-Tao Lambda=0 boundary, made concrete.")
