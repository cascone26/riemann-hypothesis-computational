"""Weil bridge: deficiency form D_ij = [zero side - smooth side] on Gaussian windows
(centers on a grid, width sigma). Min-eig margin vs sigma; compare critical sigma to the
rigidity contraction scale delta* ~ 0.12 (t-units ramp width)."""
import numpy as np
from scipy.special import loggamma as slg
g = np.load("zeros10k.npy")
T0, T1 = 2000.0, 8000.0
gz = g[(g>T0)&(g<T1)]
tgrid = np.arange(T0+50, T1-50, 3.0)   # window centers
tt = np.arange(T0, T1, 0.01)
rho = 0.5*np.log(tt/(2*np.pi))/np.pi   # smooth density
for sigma in (2.0, 1.0, 0.5, 0.25, 0.12, 0.06):
    # phi_i(t) = exp(-(t-c_i)^2/(2 sigma^2)); zero side sum, smooth side integral
    # subsample centers to keep matrix small: spacing >= sigma*2 for near-orthogonality... keep all, matrix ~ 2000
    C = tgrid[::max(1, int(1.0/3.0))]
    C = tgrid[:600]
    Z = np.exp(-(gz[None,:]-C[:,None])**2/(2*sigma**2))      # (nc, nz)
    A = Z @ Z.T                                              # zero-side Gram
    # smooth side: integral of phi_i phi_j rho dt — Gaussian product formula
    Cd = C[:,None]-C[None,:]; Cm = 0.5*(C[:,None]+C[None,:])
    rho_m = 0.5*np.log(Cm/(2*np.pi))/np.pi
    Ssm = np.sqrt(np.pi)*sigma*np.exp(-Cd**2/(4*sigma**2))*rho_m
    D = A - Ssm
    ev = np.linalg.eigvalsh(0.5*(D+D.T))
    print(f"sigma={sigma:5.2f}: min-eig {ev[0]:+.4f}  max {ev[-1]:+.4f}  (neg count {np.sum(ev<-1e-9)}/{len(ev)})", flush=True)
print("DONE")
