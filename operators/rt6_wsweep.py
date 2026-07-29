import numpy as np, l2_edge_sweep as L
BAND=(14.0,44.0); PRIMES=(2,3,5,7,11,13,17); h=1e-4
def scale_for(p,eps):
    s=np.ones_like(L.LAMW)
    for i,n in enumerate(L.NS):
        m=n
        while m%p==0: m//=p
        if m==1: s[i]=eps
    return s
def lmin(cs,w,sc): return float(np.linalg.eigvalsh(L.build_forms(cs,w,lam_scale=sc)['full'])[0])
print("gradient dλmin/dε_p vs probe width w (deep s=6.0). does high-p couple as w sharpens?")
print(f"{'w':>4}  " + "  ".join(f"p={p}" for p in PRIMES))
for w in (4.0,3.0,2.0,1.5,1.0,0.7):
    cs=L.probes(BAND,6.0)
    row=[]
    for p in PRIMES:
        g=(lmin(cs,w,scale_for(p,1+h))-lmin(cs,w,scale_for(p,1-h)))/(2*h)
        row.append(g)
    print(f"{w:4.1f}  " + "  ".join(f"{g:+.3f}" for g in row), flush=True)
print("\n1/w = probe bandwidth in r; couples to primes with log p <~ 1/w * O(1).")
print("log p:", {p: round(np.log(p),2) for p in PRIMES})
