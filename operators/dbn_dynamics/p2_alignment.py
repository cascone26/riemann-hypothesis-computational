#!/usr/bin/env python3
"""§188: the explicit alignment condition. full = arch + w·J, J involution (eigenvalues ±w, w=Λ(2)/√2).
Split J's ±w eigenspaces V±. Diagonal prediction: an arch mode u with arch-value α and +w-overlap
p=‖P₊u‖² gets lifted to α + w(2p−1). λ_min(full) is the near-cancellation of arch's 2 negative modes
(α₁,α₂) against their +w overlaps. Compute, across a∈(½log2,½log3):
  - the 2 arch negative eigenvalues α₁,α₂ and their +w-eigenspace overlaps p₁,p₂
  - the diagonal-predicted lift α_i + w(2p_i−1) vs the ACTUAL full λ_min (how much off-diagonal mixing?)
  - the deficit D(a) = w − |α₁| (perturbation-norm minus arch-negativity-depth) — the room available.
If p_i and D(a) have clean a-dependence, the p=2-only inequality is pinned to an explicit scalar form.
Reuses weil_matrices (mpmath). NO zeros involved — a finite explicit inequality, not RH-equivalent."""
import sys, math
import mpmath as mp
sys.path.insert(0, '/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import weil_matrices, gen_eig_sym

def massinner(x, y, Mass, K):
    return mp.fsum(x[i] * Mass[i, j] * y[j] for i in range(K) for j in range(K))

if __name__ == "__main__":
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    K   = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    M   = int(sys.argv[3]) if len(sys.argv) > 3 else 600
    avals = [mp.mpf(x) for x in (sys.argv[4].split(',') if len(sys.argv) > 4 else ['0.45','0.50','0.549'])]
    mp.mp.dps = dps
    w = mp.log(2) / mp.sqrt(2)
    print(f"§188 p=2 alignment | dps={dps} K={K} M={M}  w=Λ(2)/√2={mp.nstr(w,8)}")
    for a in avals:
        Qf, Qa, Mass, Vv, xs, ww = weil_matrices(a, K, M)
        Ea, Va = gen_eig_sym(Qa, Mass)                 # arch eigenpairs (Mass-orthonormal columns)
        Qp = mp.matrix(K, K)
        for i in range(K):
            for j in range(K): Qp[i, j] = Qf[i, j] - Qa[i, j]
        Ep, Vp = gen_eig_sym(Qp, Mass)                 # p2-part eigenpairs; +w eigenspace = columns with Ep>0
        Ef, _ = gen_eig_sym(Qf, Mass)
        plus = [k for k in range(K) if Ep[k] > 0]      # V₊ index set
        def col(V, k): return [V[r, k] for r in range(K)]
        def overlap_plus(u):                            # ‖P₊u‖² in Mass metric (u Mass-normalized)
            return mp.fsum(massinner(u, col(Vp, k), Mass, K) ** 2 for k in plus)
        print(f"\n=== a={mp.nstr(a,4)}  full λ_min={mp.nstr(Ef[0],6)}  D=w−|α₁|={mp.nstr(w-abs(Ea[0]),6)} ===")
        for i in (0, 1):                                # the two arch negative modes
            u = col(Va, i); alpha = Ea[i]
            p = overlap_plus(u)
            pred = alpha + w * (2 * p - 1)
            print(f"  arch mode {i}: α={mp.nstr(alpha,6)}  +w-overlap p={mp.nstr(p,6)}  diag-lift α+w(2p−1)={mp.nstr(pred,6)}")
    print("\nREAD: if the diag-lift ≈ full λ_min, positivity ≈ (α_i + w(2p_i−1) ≥ 0) — an explicit scalar per mode.")
    print("If diag-lift ≫ full λ_min, the 2 negative modes MIX through J (off-diagonal) — the near-cancellation is 2D.")
    print("Either way it is a FINITE inequality in {α_i(a), p_i(a), w} — no zeros; the a-dependence is the proof target.")
