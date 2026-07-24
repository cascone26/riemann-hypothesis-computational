import numpy as np, os
# de Bruijn-Newman zero dynamics: dx_j/dt = 2 sum_{k!=j} 1/(x_j-x_k). Flow the actual zeta zeros.
# FORWARD (t>0): spread/stay real. BACKWARD (t<0): compress -> collide at Lambda. Look for a MONOTONE (Lyapunov) functional.
G=[]
for fn in ['data/zeros_100k.txt']:
    p=os.path.join('/Users/scones/projects/riemann',fn)
    if os.path.exists(p):
        with open(p) as f:
            for line in f:
                line=line.strip()
                if line:
                    try: G.append(float(line.split()[0]))
                    except: pass
                if len(G)>=350: break
G=np.array(sorted(set(G))[:250]); N=len(G)
print(f"loaded {N} zeros, first={G[0]:.4f}, last={G[-1]:.2f}")
X0=np.concatenate([-G[::-1],G])   # symmetric +/- system
def deriv(X):
    D=X[:,None]-X[None,:]; np.fill_diagonal(D,np.inf); return 2*np.sum(1.0/D,axis=1)
def flow(X,dt,ns):
    X=X.copy()
    for _ in range(ns):
        k1=deriv(X);k2=deriv(X+dt/2*k1);k3=deriv(X+dt/2*k2);k4=deriv(X+dt*k3)
        X=X+dt/6*(k1+2*k2+2*k3+k4)
    return X
def diag(X):
    pos=np.sort(X[X>0]); gap=np.min(np.diff(pos))
    wmean=np.sum(2/pos**2)
    D=np.abs(X[:,None]-X[None,:]); np.fill_diagonal(D,1.0); inter=0.5*np.sum(np.log(D))
    # Rodgers-Tao-style energy H = -sum_{j<k} log|x_j-x_k| + (1/2)sum x_j^2 -type; report -inter and second moment
    m2=np.sum(X**2)
    return gap,wmean,inter,m2
print("\n t        min-gap    Weil-mean(sum 2/g^2)   log-interaction   2nd-moment(sum x^2)")
for t in [-0.03,-0.02,-0.01,-0.005,0.0,0.005,0.01,0.02,0.03,0.05,0.08]:
    ns=max(1,int(abs(t)/0.0005)); X=flow(X0,np.sign(t)*0.0005,ns) if t!=0 else X0
    g,wm,it,m2=diag(X)
    print(f"{t:+.3f}   {g:.5f}   {wm:.6f}         {it:12.3f}   {m2:.1f}", flush=True)
print("\nwhich functional is MONOTONE in t across 0? (Lyapunov candidate for controlling Lambda). min-gap->0 backward = Lambda.")
