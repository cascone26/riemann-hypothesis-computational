#!/usr/bin/env python3
"""§411 (ATTACK, positivity-bound family): the RH-content per §404/§410 = 'true primes keep B_analytic PSD
to the rank edge.' Test done RIGHT on the ANALYTIC side (primes-only), across the plunge, with quantified
margin. Discriminate: (i) TRUE Lambda, (ii) bent Lambda(2)x1.05, (iii) Poisson-random prime weights same
density, (iv) all-Lambda-scaled. Question: is 'edge margin >= 0' special to TRUE, and is the true edge margin
bounded below by an analytic formula (two-floor / gradient-floor §387)?"""
import numpy as np
import l2_edge_sweep as L
rng = np.random.default_rng(5)
w = 3.0
scale_bent = np.ones_like(L.LAMW)
for i,nn in enumerate(L.NS):
    m=nn
    while m%2==0: m//=2
    if m==1: scale_bent[i]=1.05
# poisson-random prime weights: keep positions, randomize magnitudes at matched mean/variance
scale_pois = np.ones_like(L.LAMW)
mask = L.LAMW>0
scale_pois[mask] = np.maximum(0, 1 + rng.normal(0,0.05,mask.sum()))   # 5% random jitter (same scale as bent)
print("edge margin (lmin at each fill) TRUE vs perturbed prime data. band[14,44] w=3.")
print(f"{'s':>5} {'n':>3} {'TRUE lmin':>11} {'bent(2)':>11} {'rand5%':>11} {'|| worst-neg ||':>14}")
for s in np.arange(6.0, 2.9, -0.3):
    cs = L.probes(L.BAND, s)
    lt = float(np.linalg.eigvalsh(L.build_forms(cs,w)['full'])[0])
    lb = float(np.linalg.eigvalsh(L.build_forms(cs,w,lam_scale=scale_bent)['full'])[0])
    lp = float(np.linalg.eigvalsh(L.build_forms(cs,w,lam_scale=scale_pois)['full'])[0])
    print(f"{s:5.1f} {len(cs):3d} {lt:+11.3e} {lb:+11.3e} {lp:+11.3e} "
          f"{'BENT NEG' if lb<-1e-12 else 'rand NEG' if lp<-1e-12 else 'all>=0':>14}", flush=True)
print("\nreading: if TRUE stays >=0 to the rank edge while BENT (and/or RAND) go NEGATIVE before it, the")
print("sign-survives-to-edge property is SPECIAL to true Lambda = the localized RH-content. The SCALE of")
print("perturbation that first breaks it (5% here) vs the two-floor bound (2C2-1 ~ 0.32 relative) is the")
print("quantitative margin. This is the analytic positivity-bound the second space would have to explain.")
