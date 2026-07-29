#!/usr/bin/env python3
"""§435: the BSLR equality-gap. d_N²logN ≥ Σ1/|ρ|²=0.046191 (proven); equality conjectural. Push N, characterize
gap(N)=d_N²logN−const: does it →0, and how? Monitor Gram conditioning (γ_n become near-dependent at large N)."""
import numpy as np
EULER=0.5772156649015329; const=2+EULER-np.log(4*np.pi)
J=1500000; j=np.arange(1,J+1,dtype=np.float64); w=1.0/(j*(j+1))
print(f"§435 BSLR gap. const={const:.6f}. J={J}, tail~{1/(4*J):.1e}")
Nmax=160
# build V incrementally to save memory: process in chunks
V=np.empty((Nmax-1,J),dtype=np.float64)
for i,n in enumerate(range(2,Nmax+1)):
    V[i]=(j%n)/n
Vw=V*w[None,:]; Gf=Vw@V.T; bf=-(Vw@np.ones(J)); one=w.sum()
print(f"{'N':>4} {'d_N^2':>11} {'d_N^2logN':>11} {'gap':>11} {'gap/const':>9} {'cond(G)':>9}")
for N in [20,40,60,80,100,120,140,160]:
    k=N-1; G=Gf[:k,:k]; b=bf[:k]
    cnd=np.linalg.cond(G)
    c=np.linalg.solve(G,b); d2=max(one-b@c,0.0)
    dl=d2*np.log(N); gap=dl-const
    print(f"{N:>4} {d2:11.6e} {dl:11.6f} {gap:+11.6f} {gap/const:9.3f} {cnd:9.1e}",flush=True)
print("""
reading: gap decreasing toward 0 with clean rate => supports the BSLR equality conjecture (d_N²logN→exactly the
const). gap flattening at >0 or cond(G) exploding => numerical limit reached / or genuine positive excess. The
gap's DECAY RATE (if clean) is the refinement content. (Note: cond(G) blows up as γ_n become near-dependent — the
numerics cap N.)
""")
