#!/usr/bin/env python3
"""Red's cheap probe (§791-addendum-3 → banks as §794, REPORT/float): is the near-kernel of the
ZEROS-ONLY sampling Gram on the height-localized class spanned by the in-band prime-node data?

CONJECTURE (three-way convergence: §476f corner / Gödel's t* undersampling / BRS joint nodes):
zeros alone undersample the band-limited class; BRS says the complementary data lives at log n.
So: near-kernel dim of S (zeros-only) should track the in-band prime-power count #{n: log n<4W},
and the prime-node functionals p_n[j] = tri(log n/(2W_time)) cos(c_j log n) should be NONdegenerate
on the near-kernel (the invisible-to-zeros directions carry exactly the prime data).

Setup: dense probe class (centers 14..44 step 1.0, 31 probes), zeros gamma<=800 (491), 4W in
{1.32, 1.45, 1.61} (prime nodes {2,3} -> {2,3,4} -> {2,3,4,5}).
S[j,k] = sum_gamma F_j(gamma)F_k(gamma); G = L2 overlap Gram; generalized spectrum of (S, G).
FLOAT = REPORT. If teeth: certify + BRS-restriction next.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import l2_fejer as fj

ZC = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'zeros500.npy'))
zeros = ZC[ZC <= 800.0]
CS = np.arange(14.0, 44.0001, 1.0)

for fourW in (1.32, 1.45, 1.61):
    fj.W = fourW/4.0
    fj.RGRID = np.linspace(-fj.BAND if hasattr(fj,'BAND') else -200, 200, 400001) if False else fj.RGRID
    n = len(CS)
    # sampling matrix at zeros and overlap Gram
    Z = np.array([fj.Frow(c, zeros) for c in CS])          # n x nz
    S = Z @ Z.T
    dr_grid = fj.RGRID
    P = np.array([fj.Frow(c, dr_grid) for c in CS])
    G = P @ P.T * (dr_grid[1]-dr_grid[0])
    # generalized eigenproblem S v = lam G v  (project out G's own near-null first for stability)
    gw, gv = np.linalg.eigh(G)
    keep = gw > gw.max()*1e-10
    B = gv[:, keep] / np.sqrt(gw[keep])                    # G^{-1/2} restricted
    lam = np.linalg.eigvalsh(B.T @ S @ B)
    lam = np.clip(lam, 0, None)
    lam_sorted = np.sort(lam)
    # in-band prime powers
    nodes = [m for m in range(2, 60) if abs(np.log(m)) < fourW and
             any(m == p**k for p in (2,3,5,7,11,13) for k in (1,2,3,4,5))]
    # gap detection: ratios of consecutive eigenvalues (smallest 8)
    small = lam_sorted[:8]
    ratios = small[1:]/np.maximum(small[:-1], 1e-300)
    kgap = int(np.argmax(ratios)) + 1 if small[0] < 1e-6*lam_sorted[-1] else 0
    print(f"\n===== 4W={fourW}  (in-band prime-power nodes: {nodes}, count {len(nodes)}) =====")
    print(f"  class dim kept: {keep.sum()}/{n} | zeros used: {len(zeros)}")
    print(f"  smallest 8 gen-eigs: {['%.3e'%x for x in small]}")
    print(f"  largest ratio-gap after index: {kgap} (near-kernel dim by gap heuristic)")
    # prime functionals on the near-kernel: p_n[j] = tri(logn/(2W_t)) cos(c_j logn); time half-width 2W
    Wt = 2*fj.W
    ew, ev = np.linalg.eigh(B.T @ S @ B)
    kk = max(kgap, len(nodes))
    K = B @ ev[:, :kk]                                     # near-kernel directions in probe coords
    Pm = []
    for m in nodes:
        u = np.log(m); tri = max(0.0, 1 - u/(2*Wt))
        Pm.append(tri*np.cos(CS*u))
    if Pm:
        A = np.array(Pm) @ K                               # nodes x kk : functionals on near-kernel
        sv = np.linalg.svd(A, compute_uv=False)
        print(f"  prime functionals on near-kernel (top-{kk}): singular values {['%.3e'%x for x in sv]}")
        print(f"  nondegenerate rank (sv > 1e-3*max): {(sv > 1e-3*sv.max()).sum()} vs node count {len(nodes)}")
