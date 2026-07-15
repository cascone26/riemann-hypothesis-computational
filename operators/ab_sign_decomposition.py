import numpy as np
from scipy.special import digamma as _digamma
import mpmath as mp
mp.mp.dps=30
TWO_PI=2*np.pi; LOGPI=float(np.log(np.pi)); LAMBDA_SQ=14
L=float(mp.log(LAMBDA_SQ)); sh=float(np.sinh(L/4)**2)
def vm_table(kmax):
    tab=[]
    for k in range(2,kmax+1):
        n,p,f,ok=k,0,2,True
        while f*f<=n:
            if n%f==0:
                if p==0:p=f
                elif p!=f:ok=False;break
                while n%f==0:n//=f
            f+=1
        if n>1:
            if p==0:p=n
            elif p!=n:ok=False
        if ok and p>0: tab.append((k,float(mp.log(p))))
    return tab
def build(N):
    VM=vm_table(LAMBDA_SQ); idx=list(range(-N,N+1)); D=len(idx)
    R=TWO_PI*N/L+80.0; G=60000; grid=np.linspace(-R,R,G); h=grid[1]-grid[0]
    w=np.full(G,h); w[0]=w[-1]=0.5*h
    Om=np.real(_digamma(0.25+0.5j*grid))-LOGPI
    E=np.empty((D,G),dtype=complex)
    for k,n in enumerate(idx):
        s=TWO_PI*n/L; d=s-grid; sm=np.abs(d)<1e-12
        E[k,sm]=np.sqrt(L); dd=d[~sm]
        E[k,~sm]=(1/np.sqrt(L))*np.exp(0.5j*dd*L)*2*np.sin(0.5*dd*L)/dd
    arch=((E*(w*Om))@E.conj().T).real/(2*np.pi)
    ns=np.array(idx,dtype=float)
    B=32*L*sh*(L**2-16*np.pi**2*np.outer(ns,ns))/np.outer(L**2+16*np.pi**2*ns**2,L**2+16*np.pi**2*ns**2)
    log_k=np.array([np.log(k) for k,_ in VM]); wts=np.array([lk*k**(-0.5) for k,lk in VM])
    P=np.zeros((D,D))
    for ii,n in enumerate(idx):
        for jj in range(ii,D):
            m=idx[jj]
            q=np.where(n==m,2*(1-np.abs(log_k)/L)*np.cos(TWO_PI*n*log_k/L),
                (np.sin(TWO_PI*m*log_k/L)-np.sin(TWO_PI*n*log_k/L))/(np.pi*(n-m)+1e-300))
            P[ii,jj]=P[jj,ii]=float(wts@q)
    return arch,B,P,idx

N=64
arch,B,P,idx=build(N)
AB=(arch+B+(arch+B).T)/2
ev,evec=np.linalg.eigh(AB)
print(f"N={N}: A+B eigenvalues (smallest 8):")
for j in range(8):
    v=evec[:,j]
    # dominant modes
    pw=v**2; top=np.argsort(pw)[-3:][::-1]
    modes=", ".join(f"n={idx[i]}({pw[i]:.2f})" for i in top)
    Pvv=float(v@P@v); ABvv=float(v@AB@v)
    print(f"  [{j}] eig(A+B)={ev[j]:+.4f}  P(v)={Pvv:+.4f}  (A+B-P)(v)={ABvv-Pvv:+.4f}  modes: {modes}")

# The 4-dim negative subspace of A+B: is P large-positive there (so -P... wait tau=A+B-P)
neg = evec[:, ev<0]  # D x 4
print(f"\n4-dim A+B-negative subspace: projected forms")
AB_sub = neg.T@AB@neg
P_sub = neg.T@P@neg
tau_sub = neg.T@(AB-P)@neg
print(f"  eig(A+B | sub) = {np.linalg.eigvalsh(AB_sub)}")
print(f"  eig(P    | sub) = {np.linalg.eigvalsh(P_sub)}")
print(f"  eig(A+B-P| sub) = {np.linalg.eigvalsh(tau_sub)}")

# On the complement (A+B >= 0 subspace), is A+B >= P?  i.e. is tau >= 0 there trivially?
pos = evec[:, ev>1e-9]
tau_pos = pos.T@(AB-P)@pos
print(f"\nOn A+B-positive subspace ({pos.shape[1]}-dim): min eig(A+B-P)={np.linalg.eigvalsh(tau_pos)[0]:+.4f}")
print(f"  min eig(A+B|pos)={np.linalg.eigvalsh(pos.T@AB@pos)[0]:+.4f}, max eig(P|pos)={np.linalg.eigvalsh(pos.T@P@pos)[-1]:+.4f}")

print("\n\n=== P tracks A+B: full-spectrum comparison (N=64) ===")
# For each A+B eigvec v_j (eigenvalue a_j), compute P(v_j). Is P(v_j) <= a_j always, gap small?
gaps=[]
Pvals=[]
for j in range(len(ev)):
    v=evec[:,j]
    Pvv=float(v@P@v)
    Pvals.append(Pvv)
    gaps.append(ev[j]-Pvv)  # = (A+B-P)(v_j), the diagonal of tau in A+B basis
gaps=np.array(gaps); Pvals=np.array(Pvals)
print(f"  min gap (A+B-P)(v_j) over all j: {gaps.min():+.5f} at j={gaps.argmin()} (eig(A+B)={ev[gaps.argmin()]:+.3f})")
print(f"  max gap: {gaps.max():+.4f}")
print(f"  ALL gaps > 0? {np.all(gaps>0)}")
print(f"  Correlation P(v_j) vs eig(A+B)_j: {np.corrcoef(ev,Pvals)[0,1]:.6f}")
# relative tracking: P(v_j)/a_j for |a_j| not tiny
mask=np.abs(ev)>0.3
print(f"  P(v_j)/(A+B)(v_j) ratio for |a_j|>0.3: min={np.min(Pvals[mask]/ev[mask]):.4f}, max={np.max(Pvals[mask]/ev[mask]):.4f}")
# The diagonal-in-AB-basis gaps sum = trace(tau) = trace(A+B) - trace(P)
print(f"  sum of gaps = trace(A+B-P) = {gaps.sum():.4f}  (check trace tau={np.trace(arch+B-P):.4f})")
# BUT: min eig(tau) is NOT min gap (off-diagonal in AB basis). Compare:
tauev=np.linalg.eigvalsh((arch+B-P+(arch+B-P).T)/2)
print(f"  min eig(tau)={tauev[0]:.2e} vs min diagonal-gap={gaps.min():.4f}")
print(f"  -> off-diagonal (in A+B eigenbasis) of P pushes min eig from {gaps.min():.4f} down to ~0")
