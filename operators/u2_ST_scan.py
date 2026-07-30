#!/usr/bin/env python3
"""Task 2 (roundtable §425/§427 handoff): scan the U(2) self-adjoint-extension family of the
CCM prolate operator (lambda=sqrt2 build) and grade EACH extension by pearsonr(F_theta, S_true),
NOT by spectrum-match. Pre-registered verdict:
  COINCIDE : some (beta,theta) gives |r| clearly > bare (r_bare~0.10) and approaching ~1
  PARTIAL  : max|r| rises meaningfully above bare but stays well below 1
  NULL     : max|r| ~ bare everywhere -> NO CCM extension carries arithmetic S(T)
             = the real Hilbert-Polya gap, stated cleanly.

Construction (matches §427's F_op computation):
  build: FP=8pi^2, w=2pi*sqrt2, singular pt x=sqrt2, alpha=-(gamma/2)^2.
  u(x)=x*xi(x) ~ A sin(w x) + B cos(w x) as x->inf  (both solutions L^2 = limit-circle at inf).
  Right-end self-adjoint BC picks asymptotic phase theta: eigen-condition A sin(theta) - B cos(theta)=0.
    theta=0 => B=0 (the bare / Friedrichs build of §427).
  Left-end (x=sqrt2) knob beta: IC = cos(beta)*regular + sin(beta)*singular(log).
    beta=0 => regular (Friedrichs) = the bare build.
  (beta,theta) sweep the physically-relevant 2 real knobs of the 4-real-param U(2); the overall
  phase and the second diagonal phase do not change the real spectrum's fluctuation shape.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import pearsonr

import os
_ZP = 'zeros500.npy' if os.path.exists('zeros500.npy') else 'operators/zeros500.npy'
Z = np.load(_ZP)
FP = 8*np.pi**2; w = 2*np.pi*np.sqrt(2); sp = np.sqrt(2.0)
X = 60.0; GMAX = 95.0
eps = 1e-6

def rhs(x, y, alpha):
    xi, xip = y
    return [xip, -(2*x*xip + (FP*x*x - alpha)*xi)/(x*x - 2)]

def AB_for(gamma, beta):
    """Return (A,B): asymptotic sin/cos coeffs of u=x*xi at x=X for extension-knob beta at x=sqrt2."""
    alpha = -(gamma/2.0)**2
    dxi0_reg = -(2*FP - alpha)/(2*sp)
    IC_reg = np.array([1.0, dxi0_reg])
    IC_sing = np.array([np.log(eps), 1.0/eps])       # leading log solution at x=sqrt2+eps
    IC = np.cos(beta)*IC_reg + np.sin(beta)*IC_sing
    xs = np.linspace(X-6, X, 120)
    s = solve_ivp(rhs, [sp+eps, X], IC, args=(alpha,), method='DOP853',
                  rtol=1e-8, atol=1e-10, t_eval=xs)
    u = xs*s.y[0]
    ph = w*xs - (gamma/2.0)**2/(2*w*xs)          # WKB tail-corrected phase (accurate absolute spectrum)
    M = np.column_stack([np.sin(ph), np.cos(ph)])
    (A, B), *_ = np.linalg.lstsq(M, u, rcond=None)
    return A, B

def roots_from_h(gs, h):
    """Linear-interp roots of h(gs) (sign changes)."""
    r = []
    for i in range(len(gs)-1):
        if h[i] == 0.0:
            r.append(gs[i])
        elif h[i]*h[i+1] < 0:
            t = h[i]/(h[i]-h[i+1])
            r.append(gs[i] + t*(gs[i+1]-gs[i]))
    return np.array(r)

def Nsmooth(T):
    return (T/(2*np.pi))*np.log(T/(2*np.pi*np.e)) + 7.0/8.0

# S_true reference fluctuation
grid = np.linspace(14.0, GMAX, 300)
Ntrue = np.array([np.searchsorted(Z, t, side='right') for t in grid])
Ssm = Nsmooth(grid)
S_true = Ntrue - Ssm
nz = int(np.sum(Z < GMAX))
mean_spacing = 2*np.pi/np.log(GMAX/(2*np.pi))

# gamma grid for the shooting (counting-level resolution is enough)
gs = np.arange(3.0, GMAX, 0.20)

beta_grid = np.linspace(0.0, np.pi, 5)[:-1]  # left-end knob (0 = bare/Friedrichs); pi==0 so drop
theta_grid = np.linspace(0.0, np.pi, 37)    # right-end asymptotic-phase knob

print(f"U(2) x S(T) scan.  build lambda=sqrt2 (FP=8pi^2, w=2pi*sqrt2).  X={X}, GMAX={GMAX}")
print(f"n_true_zeros(<{GMAX})={nz}, mean_spacing~{mean_spacing:.3f}")
print(f"beta grid {len(beta_grid)} pts, theta grid {len(theta_grid)} pts.  Pre-reg: COINCIDE/PARTIAL/NULL\n")

results = []
for bi, beta in enumerate(beta_grid):
    # precompute A,B along gs for this beta (the heavy part; one pass)
    AB = np.array([AB_for(g, beta) for g in gs])
    A = AB[:,0]; B = AB[:,1]
    best_r_this = (0.0, None)
    for theta in theta_grid:
        h = A*np.sin(theta) - B*np.cos(theta)
        roots = roots_from_h(gs, h)
        roots = roots[(roots >= 3.0) & (roots < GMAX)]
        if len(roots) < 5:
            continue
        Nop = np.searchsorted(np.sort(roots), grid, side='right')
        Fop = Nop - Ssm
        if np.std(Fop) < 1e-9:
            continue
        r, _ = pearsonr(Fop, S_true)
        rms = np.sqrt(np.mean((Fop - S_true)**2))
        results.append((beta, theta, r, rms, len(roots)))
        if abs(r) > abs(best_r_this[0]):
            best_r_this = (r, theta)
    print(f"beta={beta:5.3f}: best |r|={abs(best_r_this[0]):.4f} at theta={best_r_this[1]}"
          if best_r_this[1] is not None else f"beta={beta:5.3f}: no valid theta", flush=True)

results.sort(key=lambda t: -abs(t[2]))
print("\n=== TOP 12 extensions by |pearsonr(F_theta, S_true)| ===")
print(f"{'beta':>7} {'theta':>7} {'r':>8} {'RMS':>7} {'nroot':>6}")
for beta, theta, r, rms, nr in results[:12]:
    print(f"{beta:7.3f} {theta:7.3f} {r:8.4f} {rms:7.3f} {nr:6d}")

# bare reference
bare = [x for x in results if abs(x[0])<1e-9 and abs(x[1])<1e-9]
maxr = max(abs(x[2]) for x in results)
print(f"\nmax|r| over ALL {len(results)} extensions = {maxr:.4f}")
print(f"bare (beta=0,theta=0) r = {bare[0][2]:.4f}" if bare else "bare not in grid")
verdict = "NULL (no extension carries arithmetic S(T) -> Hilbert-Polya gap)" if maxr < 0.35 \
    else ("COINCIDE" if maxr > 0.7 else "PARTIAL")
print(f"PRE-REGISTERED VERDICT: {verdict}")
