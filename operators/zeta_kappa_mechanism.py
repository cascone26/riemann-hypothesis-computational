import numpy as np
# zeta zeros (heights gamma), on the line rho=1/2+i*gamma
G = np.array([float(x) for x in open("operators/_zeros150.txt").read().split()][:40])

def kappa_of(u, res, zpts, tol=1e-9):
    def N(z): return sum(r/(ui-z) for ui,r in zip(u,res))
    Z=np.array(zpts); mv=np.array([N(z) for z in Z]); n=len(Z); K=np.empty((n,n),complex)
    for i in range(n):
        for j in range(n):
            d=Z[i]-np.conj(Z[j]); K[i,j]=np.imag(mv[i])/np.imag(Z[i]) if abs(d)<1e-12 else (mv[i]-np.conj(mv[j]))/d
    K=0.5*(K+K.conj().T); ev=np.linalg.eigvalsh(K); s=max(1.0,abs(ev).max())
    return int(np.sum(ev<-tol*s))

def build_zeta(off_pairs=0, off_beta=0.30, jitter=0.0, seed=0):
    """u-heights: on-line zeros at real +/-gamma; pole s=1,s=0 -> conjugate pair at -/+ i/2;
    off_pairs off-line zeros at beta != 1/2 -> conjugate pair at gamma +/- i(beta-1/2)."""
    r=np.random.default_rng(seed); u=[]; res=[]
    for g in G:
        gg = g*(1+jitter*r.standard_normal())
        u += [ gg, -gg ]; res += [1.0,1.0]        # +/- gamma both real (on line)
    # ζ pole at s=1 -> z=-i/2, FE image s=0 -> z=+i/2 : a conjugate pair OFF the real axis
    u += [ 0.5j, -0.5j ]; res += [1.0,1.0]
    # off-line zeros (RH violations): beta != 1/2 -> height leaves the real axis, conjugate pair
    for k in range(off_pairs):
        g = G[k]; d = off_beta-0.5
        u += [ g + 1j*d, g - 1j*d ]; res += [1.0,1.0]
    return u,res

# test points near real axis (resolution) in UHP
r=np.random.default_rng(7)
zp=[complex(x,y) for x,y in zip(r.uniform(-60,60,90), r.uniform(0.05,4.0,90))]

print("="*80)
print("§983 ζ-SIDE mechanism test (corrected Herglotz + conjugate-pair counter, real ζ zeros)")
print("="*80)
print(f"  on-line only (all zeros real + ζ pole)        -> κ = {kappa_of(*build_zeta(0),zp)}   (expect 1 = the pole)")
for k in [1,2,3]:
    print(f"  + {k} off-line zero-pair(s) at β=0.30           -> κ = {kappa_of(*build_zeta(k),zp)}   (expect {1+k})")
print("-"*80)
print("  ROBUSTNESS of on-line κ=1 under ±0.5% zero-height jitter (30 trials):")
vals=[kappa_of(*build_zeta(0,jitter=0.005,seed=s),zp) for s in range(1,31)]
print(f"    min={min(vals)} max={max(vals)} stable={len(set(vals))==1 and vals[0]==1}")
print("  ROBUSTNESS of κ=2 (1 off-line pair) under ±0.5% jitter (30 trials):")
vals2=[kappa_of(*build_zeta(1,jitter=0.005,seed=s),zp) for s in range(1,31)]
print(f"    min={min(vals2)} max={max(vals2)} stable={len(set(vals2))==1 and vals2[0]==2}")
print("="*80)
