"""T3 exp3: BUILD the Yakaboylu Hilbert-Pólya operator (arXiv:2408.15135, adversarially-verified recipe) and
test its spectrum vs real ζ zeros. R = -D_hat - i*mu(T), D_hat=½(px+xp) Berry-Keating, T=½(x²+p²), mu(T)=T tanh(T/2)-I.
On L²([0,L]) with Dirichlet ψ(0)=0. Point spectrum should be {i(½-s): Λ(s)=0} = {γ_k} (REAL) if it works.
HONEST: matching = the construction works (Handle); it does NOT prove RH (RH hides in W≥0/simplicity)."""
import numpy as np
from scipy.linalg import eigh, eig
G = np.array([float(x) for x in open("operators/_zeros150.txt").read().split()][:40])

def build_R(L, N):
    t = np.linspace(L/N, L, N)          # interior grid, Dirichlet ψ(0)=0 (exclude t=0)
    dt = t[1]-t[0]
    # d/dt central difference (antisymmetric) -> p = -i d/dt is hermitian
    D = np.zeros((N,N))
    for j in range(N):
        if j+1<N: D[j,j+1]=1/(2*dt)
        if j-1>=0: D[j,j-1]=-1/(2*dt)
    D = (D - D.T)/2                      # enforce antisymmetric
    p = -1j*D
    x = np.diag(t.astype(complex))
    Dhat = 0.5*(p@x + x@p)               # Berry-Keating ½(px+xp), hermitian
    Ttil = 0.5*(x@x + p@p)               # ½(x²+p²), hermitian pos
    Ttil = (Ttil + Ttil.conj().T)/2
    w,V = eigh(Ttil)                     # matrix function mu(T)=T tanh(T/2)-1
    w = np.clip(w,1e-9,None)
    mu_w = w*np.tanh(w/2) - 1.0
    muT = (V*mu_w)@V.conj().T
    R = -Dhat - 1j*muT
    return R

for (L,N) in [(30,400),(40,600)]:
    R = build_R(L,N)
    ev = eig(R, right=False)
    # eigenvalues should be ~real (imag small) and = γ_k; filter real-ish, positive, in-window
    realish = ev[np.abs(ev.imag) < 0.5*np.abs(ev.real).clip(0.1,None)]
    cand = np.sort(realish.real[(realish.real>1)&(realish.real<L)])
    # match to real γ
    matched=[]
    for g in G[G<L-3]:
        if len(cand)==0: break
        k=np.argmin(np.abs(cand-g)); matched.append(abs(cand[k]-g))
    print(f"L={L} N={N}: {len(cand)} real-ish eigenvalues in ({1},{L}); "
          f"first few: {np.round(cand[:8],2)}")
    print(f"    real γ:       {np.round(G[:8],2)}")
    if matched: print(f"    mean |eig-γ| over first {len(matched)} zeros: {np.mean(matched):.3f}   median: {np.median(matched):.3f}")
    print(f"    max |Im(eig)| among candidates: {np.max(np.abs(realish.imag)) if len(realish) else float('nan'):.3f}")
