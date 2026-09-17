"""
Extract the Weil-form minimizer w*(x) on [-a,a] (lowest eigenvector of Q_full) and test the
CCM/VJ claim that it's a PROLATE SPHEROIDAL wave function. Compare w* to prolate PSWF ps_0.
"""
import numpy as np
from numpy.polynomial import legendre
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def minimizer(a, M=4000, K=90):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gf=np.array([screw_g(d,Lam,Nmax) for d in diffs])
    idx=np.arange(M+1); shift=(idx[:,None]-idx[None,:])+M
    Gf=gf[shift]
    kk=np.arange(1,K+1)
    arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    Wp=Vp*w[:,None]
    Qf=Wp.T@Gf@Wp; Qf=(Qf+Qf.T)/2
    ev,V=np.linalg.eigh(Qf/a); c=V[:,0]
    # reconstruct w*(x)=sum c_k sin(k pi (x+a)/2a)
    Vb=np.sin(np.pi*np.outer((x+a),kk)/(2*a))
    wf=Vb@c
    if wf[np.argmax(np.abs(wf))]<0: wf=-wf   # sign convention
    return x, wf, ev[0]

def prolate_pswf_ground(a, c_bw, M):
    """Ground prolate spheroidal wave fn on [-a,a] with bandwidth c_bw, as top eigenvector of the
    concentration (sinc) kernel restricted to [-a,a]."""
    x=np.linspace(-a,a,M); h=(2*a)/(M-1)
    D=x[:,None]-x[None,:]
    Kk=np.where(np.abs(D)<1e-12, c_bw/np.pi, np.sin(c_bw*D)/(np.pi*D))*h
    ev,V=np.linalg.eigh(Kk)
    g=V[:,-1]
    if g[np.argmax(np.abs(g))]<0: g=-g
    return x, g/np.sqrt(np.sum(g*g)*h)

if __name__=="__main__":
    a=1.0
    x,wf,lam=minimizer(a, M=3000, K=90)
    wf=wf/np.sqrt(np.trapezoid(wf*wf,x))
    print(f"a={a}, lambda_a={lam:.3e}")
    print(f"w* even? symmetry err = {np.max(np.abs(wf-wf[::-1])):.3e} (small=even)")
    print(f"w* nodes (sign changes interior): {np.sum(np.diff(np.sign(wf))!=0)}")
    # best-fit prolate bandwidth by overlap
    xg=x
    best=(-1,None,None)
    for c_bw in np.linspace(0.5,12,40):
        xp,g=prolate_pswf_ground(a,c_bw,len(x))
        gi=np.interp(xg,xp,g)
        gi=gi/np.sqrt(np.trapezoid(gi*gi,xg))
        ov=abs(np.trapezoid(wf*gi,xg))
        if ov>best[0]: best=(ov,c_bw,gi)
    ov,c_bw,gi=best
    print(f"best prolate overlap |<w*,ps_0>| = {ov:.6f} at bandwidth c={c_bw:.3f}")
    print(f"(overlap ~1 => w* IS the ground prolate spheroidal wave function, per CCM/VJ)")
    # residual shape
    print(f"max|w*-ps_0| after norm = {np.max(np.abs(wf-gi*np.sign(np.trapezoid(wf*gi,xg)))):.4f}")

def sweep():
    print("\n--- c(a) scaling: does bandwidth scale cleanly with a? ---")
    print(f"{'a':>5} {'lambda_a':>12} {'best c':>8} {'overlap':>9} {'c/a':>7} {'c/a^2':>8}")
    for a in [0.5,0.75,1.0,1.5,2.0,3.0]:
        M=min(6000,max(2500,int(500/a)))
        x,wf,lam=minimizer(a,M=M,K=90)
        wf=wf/np.sqrt(np.trapezoid(wf*wf,x))
        best=(-1,None)
        for c_bw in np.linspace(0.5,20,60):
            xp,g=prolate_pswf_ground(a,c_bw,len(x))
            gi=np.interp(x,xp,g); gi=gi/np.sqrt(np.trapezoid(gi*gi,x))
            ov=abs(np.trapezoid(wf*gi,x))
            if ov>best[0]: best=(ov,c_bw)
        ov,c_bw=best
        print(f"{a:5.2f} {lam:12.3e} {c_bw:8.3f} {ov:9.5f} {c_bw/a:7.3f} {c_bw/a**2:8.3f}",flush=True)
