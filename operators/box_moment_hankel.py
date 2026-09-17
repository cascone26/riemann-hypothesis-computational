import numpy as np
# box-value moment / Hankel test (cycle 329): box(rho)=rho(1-rho)=1/4+gamma^2 (the SS1161 correct positive object,
# unlike the dead lambda_n Li coeffs). m_k=sum_rho box^{-k}. Hankel-PSD => a positive measure with support {box}
# exists => a self-adjoint Jacobi (HP) operator with that spectrum exists. CIRCULAR if m_k computed from on-line
# zeros (box trivially real+ under RH). The REAL test: compute m_k UNCONDITIONALLY from Xi's log-derivative.
G=np.load('zeros2000.npy').astype(float); G=G[G>0]; box=0.25+G**2
def m(k): return 2.0*np.sum(box**(-k))
for n in range(2,8):
    H=np.array([[m(i+j) for j in range(n)] for i in range(1,n+1)]); ev=np.linalg.eigvalsh(H)
    print(f'n={n}: min-eig={ev.min():+.3e}  {"PSD" if ev.min()>-1e-9*abs(ev.max()) else "NOT PSD"}')
