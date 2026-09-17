# Degree-2 Weil criterion BACKWARD direction (cycle 253): the genuine Ramanujan-Delta Weil form
# DETECTS off-line zeros. On-line zero => positive rank-1 phihat phihat^T; off-line zero at
# 1/2+delta+i*gamma => 2 Re[phihat_a(g-i d) phihat_b(g-i d)] = 2[Re(w)Re(w)^T - Im(w)Im(w)^T],
# a rank-2 INDEFINITE term. With a dense Gaussian basis the on-line form is marginally PSD (the
# RH razor); injecting any off-line zero (delta>0) drives min-eig negative, monotone in delta.
import numpy as np
Z=np.load('Delta_zeros_complete.npy')   # 228 on-line ordinates of L(s,Delta)
SQRT2PI=np.sqrt(2*np.pi)
def phihat(a,z): return SQRT2PI/2*(np.exp(-(z-a)**2/2)+np.exp(-(z+a)**2/2))
A=np.arange(6.0,27.0,1.0); n=len(A)      # dense basis (21 test functions)
def Q0():
    M=np.zeros((n,n))
    for g in Z: v=phihat(A,g); M+=2*np.outer(v,v)
    return M
Qb=Q0()
def Qghost(kidx,delta):
    gk=Z[kidx]; M=Qb.astype(complex); v=phihat(A,gk); M-=2*np.outer(v,v)
    w=phihat(A,gk-1j*delta); M+=2*np.real(np.outer(w,w)); return np.real(M)
if __name__=="__main__":
    print(f"baseline (all on-line, RH razor): min-eig={np.linalg.eigvalsh(Qb).min():+.3e} n_neg={(np.linalg.eigvalsh(Qb)<-1e-9).sum()}")
    for kidx in [0,1,2]:
        print(f"\n zero #{kidx} (g={Z[kidx]:.2f}) off-line:")
        for d in [0.0,0.1,0.2,0.3,0.4,0.5]:
            ev=np.linalg.eigvalsh(Qghost(kidx,d))
            print(f"   delta={d:.2f} (Re s={0.5+d:.2f}): min-eig={ev.min():+.4e} n_neg={(ev<-1e-9).sum()}")
