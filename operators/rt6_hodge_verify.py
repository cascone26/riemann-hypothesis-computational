#!/usr/bin/env python3
"""§414c VERIFY the §414b claim: is the ONE negative mode of (arch−prime) aligned with the pole anchor,
and is full λmin≈0 a genuine anchor-balances-primitive criticality or just a low-probe (no-zeros-below-14)
degeneracy? Checks: (1) |⟨neg_evec(arch−prime), anchor⟩| (alignment; ~1 => exact 1-D Hodge structure).
(2) where full's near-zero eigvec lives (low probes=pole-region artifact, vs spread=real). (3) scalar
inequality pole_λ vs |neg_λ| along the anchor. (4) does the structure SURVIVE removing the empty low band
below the first zero (probes only where zeros exist, but still wide enough to see pole tail)?"""
import numpy as np, l2_edge_sweep as L
def whiten_full(M,G,rc=1e-12):
    s,V=np.linalg.eigh(G); keep=s>rc*s[-1]; Vs=V[:,keep]/np.sqrt(s[keep]); return Vs.T@M@Vs, Vs
def poleterm(cs,w):
    fp=np.exp(-((0.5j-cs)**2)/(2*w*w))+np.exp(-((0.5j+cs)**2)/(2*w*w)); return 2.0*np.real(np.outer(fp,fp))
w=3.0
for BAND in [(2.0,44.0),(8.0,44.0),(11.0,44.0)]:
    print(f"\n=== band {BAND} w={w} ===")
    print(f"{'s':>5} {'n':>3} {'|⟨neg,anchor⟩|':>13} {'pole_λ':>9} {'|neg_λ(a-p)|':>11} {'full_λmin':>10} "
          f"{'0-evec@lowprobes':>15}")
    for s in (6.0,5.0,4.0):
        cs=L.probes(BAND,s)
        if len(cs)<4: continue
        f=L.build_forms(cs,w); G=f['overlap']
        pole=poleterm(cs,w); R=(f['archpole']-pole)-f['prime']
        Rw,Vs=whiten_full(R,G); Pw,_=whiten_full(pole,G); Fw,_=whiten_full(f['full'],G)
        # anchor = top eigvec of whitened pole; neg mode = bottom eigvec of whitened arch-prime
        pev,pvec=np.linalg.eigh(Pw); anchor=pvec[:,-1]
        rev,rvec=np.linalg.eigh(Rw); negmode=rvec[:,0]
        align=abs(float(anchor@negmode))
        # full near-zero eigvec: is it concentrated on low (pole-region) probes? map back to probe space
        fev,fvec=np.linalg.eigh(Fw); z=fvec[:,0]              # whitened coords
        zc=Vs@z                                               # back to probe coefficients
        lowmask=cs<14.0                                       # probes below first zero
        low_weight=float(np.sum(zc[lowmask]**2)/np.sum(zc**2)) if lowmask.any() else 0.0
        print(f"{s:5.1f} {len(cs):3d} {align:13.4f} {pev[-1]:9.3f} {abs(rev[0]):11.3e} {fev[0]:+10.2e} "
              f"{low_weight:15.3f}",flush=True)
print("""
reading:
  |⟨neg,anchor⟩|~1 => the single negative primitive mode IS the anchor direction: exact 1-D Hodge structure,
    full≥0 ⟺ scalar pole_λ ≥ |neg_λ|. |⟨⟩|~0 => they're unrelated, §414b was projector coincidence.
  full_λmin≈0 with 0-evec low_weight~1 => it's the empty-low-band degeneracy (artifact); if it PERSISTS on
    band (11,44) (no empty probes, first zero 14.13 just above) with the SAME alignment => the structure is real.
""")
