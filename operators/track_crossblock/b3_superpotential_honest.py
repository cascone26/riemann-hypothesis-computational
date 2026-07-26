#!/usr/bin/env python3
"""B3-done-right (§212 follow-up): extract the SUSY superpotential W(x) = -(log psi0)'(x) HONESTLY —
no regularization, no finite-difference smoothing (the workflow's 'clean ζ'/ζ ✓' was spline+Tikhonov+
Savitzky-Golay IMPOSED smoothness on 38%-trimmed data — an artifact).

Clean method: psi0 is the smallest-eigenvalue eigenvector of the certified Groskin tau in the FOURIER
basis, so psi0(x)=Σ_n v0[n] e^{2πinx/L} is a Fourier series → differentiate ANALYTICALLY term-by-term.
  W(x) = -psi0'(x)/psi0(x),  psi0'(x)=Σ_n v0[n](2πin/L) e^{2πinx/L}.
No smoothing anywhere. From §211c the smallest eigenvalue is ISOLATED (1.5e-37, 3 orders below next),
so v0 is well-defined at high precision (A2's '23-dim degenerate' was a float64 floor artifact).

Decisive skeptical questions, in order:
 1) Is v0 concentrated (low-freq ⟹ smooth psi0) or spread? Is it even/odd symmetric?
 2) Is the RAW analytic W(x) a clean function, or wild? Where are psi0's zeros (W's poles)?
 3) Only then: does W resemble ζ'/ζ, the archimedean digamma -½ψ(1/4+..), or neither? (Don't assume.)
"""
import importlib.util, os
import numpy as np
import mpmath as mp

GPATH = os.path.expanduser('~/projects/riemann/external/groskin_2607.02828/code/arb_ldlt_certify.py')
spec = importlib.util.spec_from_file_location('groskin_cert', GPATH)
g = importlib.util.module_from_spec(spec); spec.loader.exec_module(g)

def to_mpf(x):
    s = str(x)
    if '+/-' in s: s = s.split('+/-')[0].strip('[] ')
    try: return mp.mpf(s)
    except Exception: return mp.mpf(float(x))

mp.mp.dps = 200
c, N, prec = 13, 18, 1400
L = float(mp.log(c))
A, DIM = g.build_arb_tau(c, N, prec)
Amp = mp.matrix(DIM, DIM)
for i in range(DIM):
    for j in range(DIM):
        Amp[i,j] = to_mpf(A[i,j])
E, Q = mp.eigsy(Amp)                      # ascending eigenvalues, Q columns = eigenvectors
lam0 = E[0]
v0 = [Q[i,0] for i in range(DIM)]         # smallest-eigenvalue eigenvector
print(f"B3-honest: Groskin tau c={c} N={N} DIM={DIM}, L=log{c}={L:.4f}, dps={mp.mp.dps}")
print(f"smallest eigenvalue lam0 = {mp.nstr(lam0,5)}   (isolated: next = {mp.nstr(E[1],5)})")

# Q1: concentration + parity of v0.  index i -> frequency n = i-N
vf = np.array([float(x) for x in v0])
absv = np.abs(vf)
order = np.argsort(absv)[::-1]
print("\nQ1 — near-null eigenvector v0, dominant modes (n, coeff):")
for i in order[:8]:
    print(f"   n={i-N:+3d}   v0={vf[i]:+.5f}")
# parity: compare v0[n] vs v0[-n]
sym = np.mean([abs(vf[N+n]-vf[N-n]) for n in range(1,N+1)])
asym = np.mean([abs(vf[N+n]+vf[N-n]) for n in range(1,N+1)])
print(f"   even-symmetry residual mean|v[n]-v[-n]| = {sym:.2e}   odd mean|v[n]+v[-n]| = {asym:.2e}")
print(f"   => v0 is {'EVEN (cosine)' if sym<asym else 'ODD (sine)'} dominant")

# Build psi0 and psi0' ANALYTICALLY (float is fine for evaluation; coeffs are the exact eigvec)
def psi0(x):
    s = 0j
    for i in range(DIM):
        n = i - N; s += vf[i]*np.exp(2j*np.pi*n*x/L)
    return s
def dpsi0(x):
    s = 0j
    for i in range(DIM):
        n = i - N; s += vf[i]*(2j*np.pi*n/L)*np.exp(2j*np.pi*n*x/L)
    return s

xs = np.linspace(0, L, 2001)
psi = np.array([psi0(x) for x in xs])
dpsi = np.array([dpsi0(x) for x in xs])
print(f"\nQ2 — psi0 on [0,L]: max|Im|={np.max(np.abs(psi.imag)):.2e} (should be ~0, real), "
      f"range Re=[{psi.real.min():.4f},{psi.real.max():.4f}]")
# zeros of psi0 (poles of W)
psr = psi.real
signch = np.where(np.diff(np.sign(psr)) != 0)[0]
zeros_x = [0.5*(xs[k]+xs[k+1]) for k in signch]
print(f"   psi0 has {len(zeros_x)} sign-change zeros on [0,L] at x≈ {[round(z,3) for z in zeros_x]}")
# raw analytic W = -psi0'/psi0
W = -(dpsi/psi).real
# smoothness of RAW W away from poles: 2nd-difference RMS on pole-free interior
mask = np.ones_like(xs, bool)
for z in zeros_x:
    mask &= (np.abs(xs-z) > 0.05)
Wm = W[mask]; xm = xs[mask]
d2 = np.diff(Wm,2)
print(f"\nQ2 — RAW analytic W(x)=-psi0'/psi0 (NO smoothing):")
print(f"   on pole-free region: W range=[{Wm.min():.3f},{Wm.max():.3f}], "
      f"RMS 2nd-diff={np.sqrt(np.mean(d2**2)):.4e} (small ⟹ genuinely smooth, NOT an artifact)")
print("   W samples (x, W):")
for x in [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 2.0, 2.3]:
    if x < L:
        w = -(dpsi0(x)/psi0(x)).real
        print(f"     x={x:.2f}  W={w:+.4f}")

# Q3 — compare to candidates (only structural, honest). ζ'/ζ on critical line and arch digamma.
print("\nQ3 — candidate comparison (structural, honest — NOT a forced fit):")
# arch digamma log-derivative candidate: the near-null ~ E(g)~ζ-multiple; its smooth arch part
# would give W_arch(x) ~ Re[ -½ ψ(1/4 + i·(x?)) ]'... position x is log-scale; the natural frequency
# height is t. Report the simplest: does W track the DIGAMMA-type slowly-varying curve, or ζ'/ζ poles?
print("   psi0 zeros ⟹ W has genuine POLES; a smooth ζ'/ζ-critical-line object has poles at zeta zeros.")
print("   #poles of W on [0,L] =", len(zeros_x), " — compare to #primes<=c and #relevant modes.")
print("   (If W is piecewise-smooth with a few poles and matches a digamma envelope between them ⟹ real")
print("    structural object; if W is jagged even between poles ⟹ the 'clean W' was a regularization artifact.)")
