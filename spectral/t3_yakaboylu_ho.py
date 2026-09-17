"""Yakaboylu R = -D_hat - i*mu(T) in the HARMONIC-OSCILLATOR basis (natural: T=½(x²+p²)=HO diagonal).
D_hat=½(px+xp) = -i/2(a²-a†²). mu(T)=T tanh(T/2)-I. Eigenvalues of R should be the γ_k (real).
Half-line Dirichlet ψ(0)=0 => odd HO states vanish at 0; also try full basis. Diagnose raw spectrum."""
import numpy as np
G = np.array([float(x) for x in open("operators/_zeros150.txt").read().split()][:40])

def build_R_ho(N):
    n = np.arange(N)
    T = n + 0.5
    a2 = np.zeros((N,N)); ad2 = np.zeros((N,N))
    for k in range(N):
        if k>=2: a2[k-2,k] = np.sqrt(k*(k-1))       # a² lowers by 2
        if k+2<N: ad2[k+2,k] = np.sqrt((k+1)*(k+2))  # a†² raises by 2
    Dhat = -0.5j*(a2 - ad2)
    mu = T*np.tanh(T/2) - 1.0
    R = -Dhat - 1j*np.diag(mu)
    return R

for N in [200, 400]:
    R = build_R_ho(N)
    ev = np.linalg.eigvals(R)
    # trim edge-of-truncation junk: keep |Re|<N/4
    ev = ev[np.abs(ev.real) < N/4]
    # sort by imaginary part to see structure
    order = np.argsort(ev.imag)
    print(f"=== N={N} : {len(ev)} eigenvalues (edge-trimmed) ===")
    print("  Im range:", round(ev.imag.min(),2), "to", round(ev.imag.max(),2),
          " | Re range:", round(ev.real.min(),2), "to", round(ev.real.max(),2))
    realish = ev[np.abs(ev.imag) < 0.3]
    rr = np.sort(realish.real[realish.real>1])
    print(f"  real-ish (|Im|<0.3), Re>1, sorted first 10: {np.round(rr[:10],3)}")
    print(f"  real γ first 10:                            {np.round(G[:10],3)}")
