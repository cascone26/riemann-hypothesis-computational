#!/usr/bin/env python3
"""hp_secular_herglotz.py -- §237: is the zero-secular-function Herglotz? = does reality come free?

From §236: zeros are vector-encoded as roots of f(z)=Sum_j xi_j/(z-lambda_j), lambda_j=2pi j/L the
Dirichlet lattice, xi = near-null eigenvector of tau.
KEY DICHOTOMY (classical, Nevanlinna/interlacing):
  - If the residues can be made POSITIVE (all same sign) in some metric, f is Herglotz -> exactly one
    real root per pole-gap, ALL roots real by interlacing -> reality (RH-trunc) is STRUCTURAL/automatic.
    This requires the roots to INTERLACE the poles (one zero per consecutive lattice gap).
  - If the zeros do NOT interlace the lattice (sparser/denser than the poles), f is structurally
    non-Herglotz -> reality is NOT automatic -> that is exactly the depth of RH.
This probe: build xi, confirm roots=zeros, then test interlacing + residue signs + the DENSITY of zeros
vs lattice. Connects to §232 (no scale-invariant tower) and Isa's diagonal law log(kb/2pi)=2pi N'(kb):
the leafwise-harmonic space must have ZERO-density (log-growing), not lattice-density (constant).
"""
import numpy as np, mpmath as mp
import importlib.util, os
spec=importlib.util.spec_from_file_location('hpc', os.path.expanduser('~/projects/riemann/operators/hp_scaling_compress.py'))
hpc=importlib.util.module_from_spec(spec); spec.loader.exec_module(hpc)
L=hpc.L; TWO_PI=hpc.TWO_PI
TRUE_ZEROS=hpc.TRUE_ZEROS

def secular_roots(xi, lam, lo, hi, ng=8000):
    xs=np.linspace(lo,hi,ng)
    f=np.array([np.sum(xi/(x-lam)) for x in xs])
    roots=[]
    for i in range(1,len(xs)):
        if np.isfinite(f[i-1]) and np.isfinite(f[i]) and (f[i-1]<0)!=(f[i]<0):
            a,b=xs[i-1],xs[i]
            # bisect avoiding poles
            if any(abs((a+b)/2-lam)<1e-6): continue
            for _ in range(60):
                m=(a+b)/2; fm=np.sum(xi/(m-lam))
                if (f[i-1]<0)!=(fm<0): b=m
                else: a=m
            roots.append((a+b)/2)
    return roots

if __name__=="__main__":
    N,R,G=32,60.0,24000
    tau,idx=hpc.build_tau(N,R,G); idx=np.array(idx)
    ev,U=np.linalg.eigh(tau)
    xi=U[:,0]                                  # smallest-eigenvalue eigenvector = near-null
    lam=TWO_PI*idx/L                           # lattice poles
    print(f"§237 secular/Herglotz: N={N}, smallest tau eig={ev[0]:.3e}, next={ev[1]:.3e}")
    # 1. roots of secular function vs zeros
    roots=[r for r in secular_roots(xi,lam,10.0,60.0) if r>10]
    print(f"\nsecular roots in [10,60]: {[round(r,3) for r in roots]}")
    zin=[z for z in TRUE_ZEROS if 10<z<60]
    print(f"true zeros in [10,60]:    {[round(z,3) for z in zin]}")
    matched=[min(abs(z-r) for r in roots) if roots else 9.9 for z in zin]
    print(f"nearest-root gap per zero: {[round(m,4) for m in matched]}  (max {max(matched):.4f})")
    # 2. interlacing test: poles vs zeros in window, and DENSITY
    poles_in=sorted([p for p in lam if 10<p<60])
    print(f"\nlattice poles in [10,60]: {len(poles_in)} of them, spacing {TWO_PI/L:.3f}")
    print(f"zeros in [10,60]:         {len(zin)}  -> zeros are {'SPARSER' if len(zin)<len(poles_in) else 'DENSER'} than lattice")
    # per pole-gap zero count
    gaps_with_zero=0; multi=0
    allpts=sorted(poles_in)
    for i in range(len(allpts)-1):
        c=sum(1 for z in zin if allpts[i]<z<allpts[i+1])
        if c>=1: gaps_with_zero+=1
        if c>=2: multi+=1
    print(f"pole-gaps containing >=1 zero: {gaps_with_zero}/{len(allpts)-1}; gaps with >=2 zeros: {multi}")
    interlace = (gaps_with_zero==len(zin) and multi==0 and len(zin)==len(allpts)-1)
    print(f"INTERLACING (one zero per gap, Herglotz-compatible)? {interlace}")
    # 3. residue signs of xi at the poles (does a same-sign metric exist trivially?)
    signs=np.sign(xi)
    nneg=int(np.sum(signs<0)); npos=int(np.sum(signs>0))
    print(f"\nxi residue signs: {npos} positive, {nneg} negative (mixed => not manifestly Herglotz in L^2 basis)")
    # 4. density law: zeros/lattice count ratio vs height (Riemann-vonMangoldt vs constant)
    print(f"\nDENSITY: N_zeros(T)~ (T/2pi)log(T/2pi) [log-GROWING]; N_lattice(T)=T*L/2pi [constant slope {L/TWO_PI:.4f}]")
    for T in [20,40,60]:
        nz=(T/TWO_PI)*np.log(T/TWO_PI); nl=T*L/TWO_PI
        print(f"  T={T}: zeros~{nz:.1f}, lattice={nl:.1f}, ratio {nz/nl:.3f}")
    print("\nREAD: mismatch (zeros sparser + non-interlacing + log-growing density) => the secular function is")
    print("STRUCTURALLY non-Herglotz in the lattice basis. Reality is NOT free; the AKM metric must reshape the")
    print("harmonic space to ZERO-density (log-growing), not the fixed lattice -- which is WHY the §232 tower has")
    print("no scale-invariant limit. The obstruction has a name now: density/Herglotz, not generic 'positivity'.")
