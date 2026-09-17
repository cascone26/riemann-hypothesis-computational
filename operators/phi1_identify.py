"""
Identify phi_1 = ground state of Q_0 (no-pole Weil form) = the single negative square.
Compare to: w* (ground state of Q_full, = prolate ps_0 at 99.6%, cycle 85) and cosh pole direction.
Reconstruct all as functions on [-a,a], compute overlaps + node structure.
"""
import numpy as np
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def boundary_kernel(t): return -4.0*(np.exp(t/2)+np.exp(-t/2)-2.0)

def get(a, M, K):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gf=np.array([screw_g(d,Lam,Nmax) for d in diffs]); gb=boundary_kernel(diffs)
    idx=np.arange(M+1); shift=(idx[:,None]-idx[None,:])+M
    Gf=gf[shift]; Gb=gb[shift]
    kk=np.arange(1,K+1)
    arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg)
    ww=np.full(M+1,h); ww[0]=ww[-1]=0.5*h
    Wp=Vp*ww[:,None]
    Qf=(Wp.T@Gf@Wp)/a; Qf=(Qf+Qf.T)/2
    Pp=(Wp.T@Gb@Wp)/a; Pp=(Pp+Pp.T)/2
    Q0=Qf-Pp
    Vb=np.sin(np.pi*np.outer((x+a),kk)/(2*a))          # basis fns (value, not deriv)
    # phi_1 = Q0 ground state; wstar = Qf ground state (coeffs in derivative-basis eigenproblem)
    _,V0=np.linalg.eigh(Q0); c_phi=V0[:,0]
    _,Vf=np.linalg.eigh(Qf); c_w=Vf[:,0]
    phi=Vb@c_phi; wst=Vb@c_w
    # normalize (L2 on [-a,a]); fix sign
    def norm(f):
        f=f/np.sqrt(np.trapezoid(f*f,x));  return f*np.sign(f[np.argmax(np.abs(f))])
    phi=norm(phi); wst=norm(wst)
    cosh=np.cosh(x/2); cosh=cosh/np.sqrt(np.trapezoid(cosh*cosh,x))
    def ov(f,g): return abs(np.trapezoid(f*g,x))
    def nodes(f): return int(np.sum(np.diff(np.sign(f))!=0))
    return x, phi, wst, cosh, ov, nodes

if __name__=="__main__":
    print("phi_1 (neg square of Q_0) structure vs w*(=ps_0) and cosh pole dir:")
    print(f"{'a':>5} {'<phi1,w*>':>10} {'<phi1,cosh>':>12} {'<w*,cosh>':>10} {'phi1 nodes':>11} {'w* nodes':>9}")
    for a in [0.5,1.0,2.0,3.0]:
        M=min(2600,max(1500,int(400/a)))
        x,phi,wst,cosh,ov,nodes=get(a,M=M,K=80)
        print(f"{a:5.2f} {ov(phi,wst):10.5f} {ov(phi,cosh):12.5f} {ov(wst,cosh):10.5f} {nodes(phi):11d} {nodes(wst):9d}",flush=True)
    print("\nIf <phi1,w*>~1: the neg square = the full minimizer = ps_0 (all the same prolate ground state).")
