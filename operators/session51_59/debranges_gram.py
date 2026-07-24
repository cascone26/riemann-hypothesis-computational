# DE BRANGES door, step 1: the HB structure = interlacing of Gram points (zeros of phase theta) with
# zeros of Z (the gammas). Gram's law fails infinitely often yet RH holds => the de Branges reproducing
# kernel must stay POSITIVE-DEFINITE across Gram failures. Test that directly.
import mpmath as mp
mp.mp.dps=30
# Gram points g_n: theta(g_n) = n*pi.  zeros gamma_n of Z (=zeta zeros imag parts).
N=250
theta=mp.siegeltheta
gram=[mp.grampoint(n) for n in range(-1,N)]      # g_{-1}, g_0, g_1, ...
gam=[mp.im(mp.zetazero(n)) for n in range(1,N+1)]
# Gram's law: g_{n-1} < gamma_n < g_n  (one zero per Gram interval). Count/locate failures.
fails=[]
for n in range(1,N):
    # Gram interval (g_{n-1}, g_n) should contain exactly gamma_n (Rosser-style). Simple test:
    lo,hi=gram[n], gram[n+1]   # gram index offset (gram[0]=g_{-1})
    cnt=sum(1 for g in gam if lo< g <hi)
    if cnt!=1: fails.append((n,float(lo),float(hi),cnt))
print(f"Gram intervals checked: {N-1}.  Gram's-law failures (interval not containing exactly 1 zero): {len(fails)}")
for f in fails[:6]: print(f"   n={f[0]}: interval ({f[1]:.3f},{f[2]:.3f}) contains {f[3]} zeros")
# de Branges reproducing kernel for HB function E with phase theta:
#   K(w,z) = [ sin(theta(z)-theta(w)) ... ]  diagonal K(t,t) = theta'(t)/pi  (phase density) > 0 always.
# The real HB test: 2x2 Gram matrix [[K(a,a),K(a,b)],[K(b,a),K(b,b)]] pos-def across points spanning a
# failure.  K(a,b) = |E(a)||E(b)| sin(theta(a)-theta(b))/(pi(a-b)) form; use phase-only model K=sin(th_a-th_b)/(a-b).
def thp(t): return float(mp.diff(theta,t))
def K(a,b):
    ta,tb=float(theta(a)),float(theta(b))
    if abs(a-b)<1e-9: return thp(a)/mp.pi
    return mp.sin(ta-tb)/(mp.pi*(a-b))
# test kernel positivity on a cluster of points straddling the first failure
if fails:
    n0=fails[0][0]; center=(fails[0][1]+fails[0][2])/2
else:
    n0=10; center=float(gram[11])
pts=[center-3,center-1.5,center,center+1.5,center+3]
G=mp.matrix(5,5)
for i in range(5):
    for j in range(5): G[i,j]=K(pts[i],pts[j])
eigs=sorted(float(e) for e in mp.eigsy(G,eigvals_only=True))
print(f"\nde Branges phase-kernel Gram matrix at 5 pts straddling {'first failure n='+str(n0) if fails else 'a Gram interval'}:")
print(f"   eigenvalues: {['%.3e'%e for e in eigs]}   (all >0 => reproducing kernel POS-DEF here)")
print(f"   min eig = {eigs[0]:.3e}  => {'POSITIVE (HB structure intact across failure)' if eigs[0]>-1e-12 else 'INDEFINITE'}")
