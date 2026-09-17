"""
Resolve cycle-85c tension: compute L(w*)/||w*||^2 (archimedean Dirichlet form, eq 4.4) on the SCALED
minimizer w on [-1,1], and check it reconciles with arch_part = L/||w||^2 - log a - (2A+1) [=g_no_primes
derivative Rayleigh, cycle 83]. Determine the growth law of L/||w||^2 vs a.
L(w)=(1/4) int int |w(x)-w(y)|^2/|x-y| dxdy - (1/2) int |w(x)|^2 log(1-x^2) dx.
"""
import numpy as np
from suzuki_screw_lambda_a import screw_g, screw_g_no_primes, vonmangoldt_upto

TWO_A1 = np.log(2*np.pi) + 0.5772156649015329   # 2A+1 = log2pi + C0

def minimizer_scaled(a, M, K):
    """lowest eigenvector of full Q on [-a,a]; return t in [-1,1], w(t)=v(a t), and arch_part."""
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gf=np.array([screw_g(d,Lam,Nmax) for d in diffs])
    ga=np.array([screw_g_no_primes(d,Lam,Nmax) for d in diffs])
    idx=np.arange(M+1); shift=(idx[:,None]-idx[None,:])+M
    Gf=gf[shift]; Ga=ga[shift]
    kk=np.arange(1,K+1)
    arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg)
    ww=np.full(M+1,h); ww[0]=ww[-1]=0.5*h
    Wp=Vp*ww[:,None]
    Qf=Wp.T@Gf@Wp; Qf=(Qf+Qf.T)/2
    Qa=Wp.T@Ga@Wp; Qa=(Qa+Qa.T)/2
    ev,V=np.linalg.eigh(Qf/a); c=V[:,0]
    mass=a*np.dot(c,c)
    arch_part=(c@Qa@c)/mass
    lam=ev[0]
    Vb=np.sin(np.pi*np.outer((x+a),kk)/(2*a))
    wf=Vb@c
    t=x/a
    return t, wf, arch_part, lam

def L_form(t, w):
    """eq 4.4 on [-1,1]. Diagonal integrand ~|x-y| -> 0, safe."""
    ht=t[1]-t[0]
    # jumping part: (1/4) int int |w(x)-w(y)|^2/|x-y|
    D=t[:,None]-t[None,:]
    Wd=w[:,None]-w[None,:]
    absD=np.abs(D); np.fill_diagonal(absD,1.0)
    integ=(Wd*Wd)/absD; np.fill_diagonal(integ,0.0)
    jump=0.25*integ.sum()*ht*ht
    # killing part: -(1/2) int |w|^2 log(1-t^2)
    kern=np.log(np.clip(1-t*t,1e-15,None))
    kill=-0.5*np.sum(w*w*kern)*ht
    norm=np.sum(w*w)*ht
    return jump+kill, norm

if __name__=="__main__":
    print("Reconcile arch_part vs L/||w||^2 - (log a + 2A+1); growth law of L/||w||^2.")
    print(f"2A+1 = {TWO_A1:.5f}")
    print(f"{'a':>5} {'arch_part':>11} {'L/||w||^2':>11} {'L/n - log a -(2A+1)':>20} {'match?':>8} {'L/n':>8}")
    for a in [0.35,0.5,0.75,1.0,1.5,2.0,3.0]:
        M=min(2400,max(1200,int(300/a)))
        t,w,arch,lam=minimizer_scaled(a,M=M,K=70)
        L,norm=L_form(t,w)
        Ln=L/norm
        pred=Ln-np.log(a)-TWO_A1
        ok="OK" if abs(pred-arch)<max(0.05,0.05*abs(arch)) else "MISMATCH"
        print(f"{a:5.2f} {arch:11.5f} {Ln:11.5f} {pred:20.5f} {ok:>8} {Ln:8.4f}",flush=True)

def boundary_kernel(t):
    """pole/boundary term of g: -4(e^{t/2}+e^{-t/2}-2)."""
    return -4.0*(np.exp(t/2)+np.exp(-t/2)-2.0)

def check_boundary():
    from suzuki_screw_lambda_a import screw_g, screw_g_no_primes, vonmangoldt_upto
    print("\nIsolate: does the pole/boundary term carry the a-growth of cycle-83 arch_part?")
    print(f"{'a':>5} {'arch_part(g_np)':>15} {'L/||w||^2':>11} {'boundary-Rayleigh':>17} {'L-log a-(2A+1)+bdry':>20}")
    for a in [0.5,1.0,2.0,3.0]:
        M=min(2400,max(1200,int(300/a)))
        t,w,arch,lam=minimizer_scaled(a,M=M,K=70)
        L,norm=L_form(t,w)
        # boundary derivative-Rayleigh on same w* (on [-a,a]): need w' ; reconstruct on [-a,a]
        x=t*a; h=x[1]-x[0]
        wp=np.gradient(w,x)
        diffs=np.arange(-(len(x)-1),len(x))*h
        gb=boundary_kernel(diffs)
        idx=np.arange(len(x)); shift=(idx[:,None]-idx[None,:])+(len(x)-1)
        Gb=gb[shift]
        ww=np.full(len(x),h); ww[0]=ww[-1]=0.5*h
        wpw=wp*ww
        bdry=(wpw@Gb@wpw)/(norm*a)   # ||w||^2 on [-a,a] = norm*a
        Ln=L/norm
        recon=Ln-np.log(a)-TWO_A1+bdry
        print(f"{a:5.2f} {arch:15.5f} {Ln:11.4f} {bdry:17.5f} {recon:20.5f}",flush=True)
    print("if recon ~ arch_part, then arch_part = L/||w||^2 - log a - (2A+1) + boundary (bookkeeping closed).")

check_boundary()
