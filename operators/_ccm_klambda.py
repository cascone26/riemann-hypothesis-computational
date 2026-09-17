import numpy as np
# Real CCM k_λ (arXiv:2511.22755 Eq 7.2,7.5,7.6, Lemma 7.1). Build h_{0,λ},h_{4,λ} from PW_λ,
# form h_λ, apply E-map k_λ(u)=√u Σ_{n≥1} h_λ(nu), transform, check zeros vs γ_k.
GAM = [14.1347,21.0220,25.0109,30.4249,32.9351,37.5862]

def prolate_efns(lam, M=4000):
    x = np.linspace(-lam, lam, M+1); h = x[1]-x[0]
    xm = 0.5*(x[:-1]+x[1:]); p = lam*lam - xm*xm
    N=M+1; A=np.zeros((N,N))
    for i in range(1,N-1):
        A[i,i-1]=-p[i-1]/h**2; A[i,i+1]=-p[i]/h**2
        A[i,i]=(p[i-1]+p[i])/h**2 + (2*np.pi*lam*x[i])**2
    A[0,0]=A[-1,-1]=1e12
    ev,evec=np.linalg.eigh(A)
    def norm(v): return v/np.sqrt(np.trapezoid(v*v,x))
    return x, norm(evec[:,0]), norm(evec[:,4])   # h_0, h_4 (even)

def build_hlam(lam, M=4000):
    x,h0,h4 = prolate_efns(lam,M)
    c4=np.sqrt(3)/2**(11/4); c0=3/2**(17/4)
    hl = c4*h4 - c0*h0
    # ensure even-sign convention h0>0 center
    if hl[len(hl)//2] < 0: hl=-hl
    from scipy.interpolate import interp1d
    hf = interp1d(x, hl, bounds_error=False, fill_value=0.0)
    return hf, x

def k_lambda(lam, M=4000, Nu=400):
    hf,x = build_hlam(lam,M)
    us = np.exp(np.linspace(-np.log(lam), np.log(lam), Nu))   # u in [1/λ, λ], log-spaced
    kv = np.zeros(Nu)
    for i,u in enumerate(us):
        nmax = int(np.floor(lam/u))
        if nmax<1: continue
        ns = np.arange(1,nmax+1)
        kv[i] = np.sqrt(u)*np.sum(hf(ns*u))
    return us, kv

def khat(us, kv, zs):
    # scaling-operator transform: k̂(z)=∫ k(u) u^{-iz} du/u ; log var t=log u
    t = np.log(us); w = np.gradient(t)
    out=[]
    for z in zs:
        out.append(np.abs(np.sum(kv*np.exp(-1j*z*t)*w)))
    return np.array(out)

for lam in [3.6, 6.0, 10.0]:
    us,kv = k_lambda(lam, M=4000, Nu=600)
    zs = np.linspace(2,45,2000)
    ah = khat(us,kv,zs)
    # find minima (zeros of k̂)
    mins=[zs[i] for i in range(1,len(zs)-1) if ah[i]<ah[i-1] and ah[i]<ah[i+1] and ah[i]<0.15*np.median(ah)]
    print(f"λ={lam}: k̂ zeros (first 6): {[round(m,2) for m in mins[:6]]}")
    print(f"       γ_k target        : {GAM}")
    if mins:
        errs=[min(abs(np.array(mins)-g)) for g in GAM[:len(mins) if len(mins)<6 else 6]]
        print(f"       nearest-γ errs     : {[round(e,2) for e in errs]}", flush=True)
