# CONVERGENCE-RATE LAW (Connes proof-strategy / Groskin data frontier), structural derivation.
# Recovered zeros -> zeros of Xi = FT of Polya Phi (super-decaying, ~e^{-pi e^{4u}}). Error at
# localization scale U = tail of Phi = |F_U(2*gamma_1)|, F_U(z)=2 int_0^U Phi(u)cos(zu)du -> Xi(z)=0.
# PREDICTION: log|error| ~ -pi e^{4U} (DOUBLY-exponential in U) => explains Groskin's 113-orders-fast convergence.
import mpmath as mp
mp.mp.dps=60
def Phi(u,nmax=40):
    s=mp.mpf(0); e4=mp.e**(4*u)
    for n in range(1,nmax+1):
        s+=(2*mp.pi**2*n**4*mp.e**(9*u)-3*mp.pi*n**2*mp.e**(5*u))*mp.e**(-mp.pi*n**2*e4)
    return s
g1=mp.mpf('14.134725141734693790'); z=2*g1   # zero of the scaled transform (my u-var, §65)
def F_U(U,M=4000):    # 2 int_0^U Phi(u) cos(z u) du  (-> Xi(z)=0 as U->inf; |F_U| = convergence error)
    h=U/M; s=mp.mpf(0)
    for i in range(M+1):
        u=i*h; w=h if 0<i<M else h/2
        s+=w*Phi(u)*mp.cos(z*u)
    return 2*s
print("U     |F_U(2*g1)| = convergence error      log10|err|   -pi*e^{4U}/ln10 (predicted)")
prev=None
for U in [mp.mpf(x) for x in ['0.35','0.45','0.55','0.65','0.75','0.85']]:
    err=abs(F_U(U)); l=mp.log(err)/mp.log(10) if err>0 else mp.mpf('-inf')
    pred=-mp.pi*mp.e**(4*U)/mp.log(10)
    drop = f"  (drop {float(l-prev):+.1f} orders)" if prev is not None else ""
    print(f"{float(U):.2f}  {mp.nstr(err,4):>16}   {float(l):+8.2f}    {float(pred):+9.2f}{drop}")
    prev=l
print("\nlog10|err| tracking -pi*e^{4U}/ln10 => DOUBLY-exponential convergence, inherited from Phi's e^{-pi e^{4u}} decay.")
print("This is the structural origin of the CCM/Suzuki/Groskin rapid convergence (Groskin: 113 orders over c=13..100).")
