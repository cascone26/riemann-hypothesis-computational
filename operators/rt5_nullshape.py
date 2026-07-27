#!/usr/bin/env python3
"""§375: WHAT IS the arithmetic-critical direction? (§374 said it's arithmetic, not frame.)

At n=n_crit, v_full = bottom eigenvector of B_full = pole+arch-prime (zero-free). Form the
test function phi(r) = sum_j v_full[j] F_j(r) (symmetrized Gaussians). Question: does phi VANISH
at the zero ordinates gamma_k in the band, while being O(1) between them? If yes, the
arithmetic-critical direction realizes the Bondarenko-Radchenko-Seip picture (§336): the
collapse mode is exactly the function that interpolates ZERO on the zero set at critical density.
Compare against v0 = overlap-null (the frame mode, §374) which should NOT specially vanish at zeros.
Metric: |phi(gamma_k)| (at zeros) vs |phi| sampled at midpoints between consecutive zeros.
Zero-side check uses B_full's own near-null structure — B_full == B_zero (explicit formula) so
v_full is also the near-null of the ZERO-side Gram 2 sum F(gamma)F(gamma)^T, i.e. the combo
orthogonal to all zero-evaluation vectors => phi(gamma_k) ~ 0. Directly measurable.
"""
import numpy as np
import l2_edge_sweep as L
from scipy.special import digamma

def set_grid(hi):
    L.RGRID = np.arange(-(hi + 60.0), hi + 60.0, 0.002)
    L.OMEGA = np.real(digamma(0.25 + 0.5j * L.RGRID)) - np.log(np.pi)

def N_rvm(T):
    if T <= 2 * np.pi: return 0.0
    return T / (2 * np.pi) * np.log(T / (2 * np.pi)) - T / (2 * np.pi) + 7.0 / 8.0

def phi_at(cs, w, coef, r):
    val = 0.0
    for j, c in enumerate(cs):
        val += coef[j] * (np.exp(-((r - c) ** 2) / (2 * w * w)) + np.exp(-((r + c) ** 2) / (2 * w * w)))
    return val

w = 2.0
# use the low bands where the hardcoded ZEROS list is complete (gamma < 101)
BANDS = [(14.0, 44.0), (50.0, 80.0)]
print("Does the arithmetic-null mode phi VANISH at the zeros? (BRS interpolation picture, §336)")
for band in BANDS:
    b0, b1 = band
    set_grid(b1 + 30.0)
    ncrit = int(round(N_rvm(b1 + 2 * w) - N_rvm(b0 - 2 * w)))
    cs = np.linspace(b0, b1, ncrit)
    f = L.build_forms(cs, w)
    lf, Vf = np.linalg.eigh(f['full'])
    vfull = Vf[:, 0]                       # arithmetic-critical direction
    lo, Vo = np.linalg.eigh(f['overlap'])
    v0 = Vo[:, 0]                          # frame-null (control)
    zeros_in = [g for g in L.ZEROS if b0 - 2 * w <= g <= b1 + 2 * w]
    # normalize phi so max over the band ~ 1
    rgrid = np.linspace(b0 - 2 * w, b1 + 2 * w, 2000)
    def norm_phi(coef):
        vals = np.array([phi_at(cs, w, coef, r) for r in rgrid])
        return vals / (np.max(np.abs(vals)) + 1e-300)
    # at zeros vs at midpoints between zeros
    zs = np.array(zeros_in)
    mids = (zs[:-1] + zs[1:]) / 2 if len(zs) > 1 else zs
    for name, coef in [("v_full(arith)", vfull), ("v0(frame)", v0)]:
        vg = np.array([phi_at(cs, w, coef, r) for r in rgrid])
        scale = np.max(np.abs(vg)) + 1e-300
        at_zeros = np.array([abs(phi_at(cs, w, coef, g)) for g in zs]) / scale
        at_mids = np.array([abs(phi_at(cs, w, coef, m)) for m in mids]) / scale
        print(f"\n  band {band} n_crit={ncrit} lmin={lf[0]:+.2e}  mode={name}")
        print(f"    |phi(zero)| mean={np.mean(at_zeros):.3f} max={np.max(at_zeros):.3f}")
        print(f"    |phi(midpt)| mean={np.mean(at_mids):.3f} max={np.max(at_mids):.3f}")
        print(f"    ratio zeros/midpts = {np.mean(at_zeros)/(np.mean(at_mids)+1e-12):.3f} "
              f"(<<1 => vanishes at zeros = interpolation)")
print("\nreading: if v_full's phi is small at zeros / O(1) at midpoints (ratio << 1) while v0's is not,")
print("the arithmetic-critical mode IS the zero-interpolating function (BRS §336) — the collapse")
print("direction is 'the Gaussian combination that hits zero exactly on the zero set'. That names the")
print("arithmetic direction §374 could not fix a-priori: it's fixed by the ZEROS (why no zero-free choice).")
