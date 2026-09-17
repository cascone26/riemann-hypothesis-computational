"""NOVEL probe on the de Bruijn-Newman reframe: does 'heating' the Weil form convert the saturation
(λ_min=0) into a strict positive margin? The heat flow damps prime contributions at log n. Model:
Q(δ) = pole + core + tail - Σ Λ(n)/√n · e^{-δ(log n)²} · S(log n)   [primes Gaussian-damped by δ≥0].
If λ_min(Q(δ))>0 for δ>0 and →0 as δ→0, the saturation IS the δ→0 (Λ→0) edge — concrete flow connection."""
import numpy as np, sys
sys.path.insert(0,'operators')
from suzuki_exact_Q import build_Q, S_mat, E_p, C_CORE, Q_tail, vonmangoldt_upto

def build_Q_heated(a,K,delta,NG=400):
    kk=np.arange(1,K+1)
    Ehp=E_p(0.5,a,kk); Ehm=E_p(-0.5,a,kk)
    Qpole=np.outer(Ehp,Ehm)+np.outer(Ehm,Ehp)
    Qcore=-2*C_CORE*a*np.eye(K)
    Qtail=Q_tail(a,K,NG)
    L=2*a; Nmax=int(np.floor(np.exp(L))); Lam=vonmangoldt_upto(max(Nmax,2))
    Qprime=np.zeros((K,K))
    for n in range(2,Nmax+1):
        if Lam[n]>0 and np.log(n)<L:
            ln=np.log(n)
            damp=np.exp(-delta*ln*ln)
            Qprime+= -(Lam[n]/np.sqrt(n))*damp*S_mat(ln,a,K)
    Q=Qpole+Qcore+Qtail+Qprime
    return (Q+Q.T)/2

for a in [1.0,1.5,2.0]:
    K=max(70,int(30*a)+40)
    print(f"a={a}:")
    for delta in [0.0, 0.001, 0.005, 0.02, 0.05, 0.1]:
        Q=build_Q_heated(a,K,delta)
        w=np.linalg.eigvalsh(Q)
        print(f"   δ={delta:5.3f}:  λ_min={w[0]/a:+.5e}   κ=#neg={int(np.sum(w<-1e-9))}   λ_2={w[1]/a:+.4e}")
print("\nRead: δ=0 is our saturated Weil form. If λ_min DECREASES below 0 as δ grows (over-damping) or")
print("INCREASES to +margin — the sign of dλ_min/dδ at δ=0 tells whether heating helps (positive margin)")
print("or hurts. Rodgers-Tao Λ≥0 predicts you CAN'T get a free positive margin — test if that shows here.")
