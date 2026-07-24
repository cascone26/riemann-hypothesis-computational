#!/usr/bin/env python3
"""§199 tail-bound analysis: compute even-sector λmin vs truncation size d=1,2,3,4,6,8,10,12,14,16,20
to quantify convergence and establish rigorous tail bounds for d>trunc. This determines the minimum
truncation size for which a rigorous tail bound closes the proof."""
import sys, math, json
import mpmath as mp
sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import weil_matrices, gen_eig_sym
mp.mp.dps=40; K=24; M=800

ev=[c for c in range(K) if (c+1)%2==1]   # odd-k modes: k=1,3,5,...

def block(Q,idx):
    """Extract a submatrix by indices."""
    B=mp.matrix(len(idx),len(idx))
    for i in range(len(idx)):
        for j in range(len(idx)): B[i,j]=Q[idx[i],idx[j]]
    return B

def eigenvalues_by_dim(a_val, K, M, max_d=20):
    """Compute even-sector λmin vs dimension d."""
    am=mp.mpf(str(a_val))
    Qf,Qa,Mass,Vv,xs,w=weil_matrices(am,K,M)

    results={}
    # Test dimensions
    test_dims=[1,2,3,4,6,8,10,12,14,16,20]
    for d in test_dims:
        if d > len(ev):
            break
        idx=ev[:d]
        E,_=gen_eig_sym(block(Qf,idx),block(Mass,idx))
        results[d]=E[0]

    # Also get the full even-sector minimum
    idx_full=ev  # all available even modes
    E_full,_=gen_eig_sym(block(Qf,idx_full),block(Mass,idx_full))
    results['full']=E_full[0]

    return results

print("§199 tail-bound convergence: even-sector λmin vs truncation dimension d")
print("="*80)
print("\nComputing for K={}, M={}, dps=40...".format(K,M))
print()

# Test across the critical a-range
test_a=[0.40, 0.44, 0.46, 0.48, 0.50, 0.52, 0.549]
all_results={}

for a in test_a:
    print(f"a = {a}: ", end='', flush=True)
    results=eigenvalues_by_dim(a, K, M)
    all_results[a]=results

    # Display as: d=1, d=2, ..., d=full
    dims=[d for d in sorted([k for k in results.keys() if k != 'full'])]

    row=[mp.nstr(results[d],4) for d in dims]
    full_val=results['full']
    row.append(mp.nstr(full_val,4))

    print("  ".join(row[:5]))
    if len(row) > 5:
        print(" "*18 + "  ".join(row[5:]))
    print()

print("\n" + "="*80)
print("CONVERGENCE TABLE (detailed): λmin(d) and tail decrement")
print("="*80)
print(f"{'a':>7} ", end='')
for d in [1,2,3,4,6,8,10,12,14,16]:
    print(f"d={d:>2}  ", end='')
print("full  Tail(d=8)  Tail(d=10) Tail(d=12)")
print("-"*100)

for a in test_a:
    if a in all_results:
        print(f"{a:>7.3f} ", end='')
        results=all_results[a]
        for d in [1,2,3,4,6,8,10,12,14,16]:
            if d in results:
                print(f"{mp.nstr(results[d],3):>7} ", end='')
        print(f"{mp.nstr(results['full'],3):>7} ", end='')

        # Estimate tail bounds: difference between d and full
        if 8 in results and 'full' in results:
            tail_8=results[8]-results['full']
            print(f"{mp.nstr(tail_8,3):>10} ", end='')
        else:
            print("   -      ", end='')

        if 10 in results and 'full' in results:
            tail_10=results[10]-results['full']
            print(f"{mp.nstr(tail_10,3):>10} ", end='')
        else:
            print("   -      ", end='')

        if 12 in results and 'full' in results:
            tail_12=results[12]-results['full']
            print(f"{mp.nstr(tail_12,3):>10}", end='')
        else:
            print("   -     ", end='')

        print()

print("\n" + "="*80)
print("TAIL BOUND ANALYSIS")
print("="*80)

for a in test_a:
    if a in all_results:
        results=all_results[a]
        print(f"\na = {a}:")
        print(f"  λmin(full) = {mp.nstr(results['full'],6)}")

        for d in [4,6,8,10,12]:
            if d in results:
                lam_d=results[d]
                tail=lam_d-results['full']
                ratio_to_full=tail/results['full'] if results['full']>0 else mp.mpf('inf')
                print(f"  λmin(d={d:2}) = {mp.nstr(lam_d,6)}  |  tail = {mp.nstr(tail,6)}  ({mp.nstr(ratio_to_full*100,5)}% of full)")

print("\n" + "="*80)
print("KEY FINDINGS:")
print("="*80)

# Find where tail becomes small relative to full margin
print("\nMin tail-size by dimension (over all tested a):")
for d in [4,6,8,10,12]:
    min_tail=None
    a_min=None
    for a in test_a:
        if a in all_results:
            results=all_results[a]
            if d in results and 'full' in results:
                tail=results[d]-results['full']
                if min_tail is None or tail < min_tail:
                    min_tail=tail
                    a_min=a
    if min_tail is not None:
        print(f"  d={d:2}: tail_min = {mp.nstr(min_tail,6)} (at a={a_min})")

# Estimate slope of convergence
print("\nConvergence slope (tail reduction per +2 modes):")
for a in [0.48, 0.50, 0.549]:
    if a in all_results:
        results=all_results[a]
        print(f"  a={a}:")
        if 8 in results and 6 in results:
            slope_68=(results[6]-results[8])/(2.0)
            print(f"    (d=6→8): {mp.nstr(slope_68,6)} reduction per mode")
        if 10 in results and 8 in results:
            slope_8_10=(results[8]-results[10])/(2.0)
            print(f"    (d=8→10): {mp.nstr(slope_8_10,6)} reduction per mode")
        if 12 in results and 10 in results:
            slope_10_12=(results[10]-results[12])/(2.0)
            print(f"    (d=10→12): {mp.nstr(slope_10_12,6)} reduction per mode")

print("\n" + "="*80)
print("RIGOROUS TAIL BOUND RECOMMENDATION")
print("="*80)
print("\nFor d>trunc, a crude rigorous bound can be established via diagonal dominance:")
print("The d>trunc tail contributes at most the sum of |off-diagonal| entries in rows d,d+2,d+4,...")
print("\nFor a rigorous proof:")
print("  • If λmin(d=4) tail ≈ O(margin):  need d=6-8 + rigorous tail bound")
print("  • If λmin(d=8) tail ≈ O(1e-7):    d=8 likely sufficient with small tail correction")
print("  • If λmin(d=12) tail ≈ o(1e-8):   d=12 strong; even simpler rigorous tail")
