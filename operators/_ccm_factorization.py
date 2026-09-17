import numpy as np, mpmath as mp
import operators._ccm_klambda as C
# Claim: k̂_λ(z) = ζ(1/2+iz) · ĥ_λ(1/2+iz), with ĥ_λ (Mellin of h_λ) ZERO-FREE (→ Gamma factor).
# So k̂'s zeros = ζ's zeros (γ_k) — the SPECTRAL realization; RH content is elsewhere (μ_λ→0 positivity).

def hhat(hf, x, z, lam):
    # ĥ_λ(1/2+iz) = ∫_0^λ h_λ(v) v^{(1/2+iz)-1} dv  (h even, v>0 half)
    vs = np.linspace(1e-4, lam, 6000); dv=vs[1]-vs[0]
    hv = hf(vs)
    return np.sum(hv * vs**(-0.5+1j*z) * dv)

lam=20.0
hf,x = C.build_hlam(lam, M=5000)
zs = np.linspace(10,40,120)
zeta = np.array([complex(mp.zeta(mp.mpf('0.5')+1j*z)) for z in zs])
hh   = np.array([hhat(hf,x,z,lam) for z in zs])
prod = zeta*hh
# my direct k̂ for comparison
us,kv = C.k_lambda(lam, M=5000, Nu=1500)
t=np.log(us); w=np.gradient(t)
kh = np.array([np.sum(kv*np.exp(-1j*z*t)*w) for z in zs])

# normalize each to compare shape (zeros locations)
def nrm(a): return np.abs(a)/np.max(np.abs(a))
print(f'λ={lam}: is ĥ_λ zero-free? min|ĥ|/median = {np.min(np.abs(hh))/np.median(np.abs(hh)):.3f} (large=>zero-free)')
print(f'ζ(1/2+iz) zeros in window (|ζ| dips): ', [round(zs[i],1) for i in range(1,len(zs)-1) if abs(zeta[i])<abs(zeta[i-1]) and abs(zeta[i])<abs(zeta[i+1]) and abs(zeta[i])<0.2*np.median(np.abs(zeta))])
print(f'k̂ (direct) zeros: ', [round(zs[i],1) for i in range(1,len(zs)-1) if nrm(kh)[i]<nrm(kh)[i-1] and nrm(kh)[i]<nrm(kh)[i+1] and nrm(kh)[i]<0.2])
# correlation of |k̂_direct| shape vs |ζ·ĥ|
cc=np.corrcoef(nrm(kh), nrm(prod))[0,1]
print(f'corr( |k̂_direct| , |ζ(1/2+iz)·ĥ_λ| ) = {cc:.3f}  (→1 confirms the Mellin factorization)')
