import numpy as np, suzuki_exact_Q as EQ
from scipy.linalg import eigh
print("Birman-Schwinger margin 1-rho(a) to high digits (a where A is PD):",flush=True)
print(f"{'a':>5} {'rho(a)':>18} {'1-rho (margin)':>16} {'lambda_a=mineig(A-P)':>20} {'sign check'}",flush=True)
K=100
rows=[]
for a in [0.5,0.8,1.1,1.4,1.7,2.0]:
    Q,pc=EQ.build_Q(a,K,return_pieces=True)
    sym=lambda M:(M+M.T)/2
    A=sym(pc["tail"]+pc["pole"]); P=sym(-(pc["core"]+pc["prime"]))
    ea=eigh(A,eigvals_only=True)
    if ea[0]<=0: print(f"{a:5.2f}  A not PD (min={ea[0]:.2e}) skip",flush=True); continue
    # rho = max generalized eig P v = rho A v
    ev=eigh(P,A,eigvals_only=True); rho=ev[-1]
    lam=eigh(A-P,eigvals_only=True)[0]
    rows.append((a,1-rho,lam))
    print(f"{a:5.2f}  {rho:.12f}  {1-rho:16.3e}  {lam:20.3e}  {'OK λ>0,ρ<1' if (lam>0 and rho<1) else 'CHECK'}",flush=True)
print("\nmargin 1-rho and lambda_a both >0 (marginal), shrinking with a = the knife-edge.",flush=True)
if len(rows)>=3:
    a_=np.array([r[0] for r in rows]); m_=np.array([r[1] for r in rows])
    print(f"log(1-rho) vs a slope (last 4): {np.polyfit(a_[-4:],np.log(m_[-4:]),1)[0]:.3f}  (negative=shrinking)",flush=True)
