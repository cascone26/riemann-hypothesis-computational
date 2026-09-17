#!/usr/bin/env python3
"""Cycle 10 (redo, correct frame): Delta_prime lift in the screw g-kernel frame (suzuki_screw_lambda_a).
G_arch from screw_g_no_primes, G_prime = G_full - G_arch. v_arch = ground eigvec of Q_arch=Wp G_arch Wp.
(1) check lambda_arch(screw) =? mu1 - log a  (reconcile with the (2.4) Dirichlet build).
(2) pairing = <c_arch|Q_prime|c_arch>/a  (mass=a*I, c normalized) = 1st-order prime lift.
(3) actual Delta = lambda_full - lambda_arch. All in ONE consistent discretization."""
import numpy as np, suzuki_screw_lambda_a as S
mu1=0.23676
def pieces(a,M=3000,K=60):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=S.vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; diffs=np.arange(-M,M+1)*h; idx=np.arange(M+1)
    gf=np.array([S.screw_g(d,Lam,Nmax) for d in diffs]); ga=np.array([S.screw_g_no_primes(d,Lam,Nmax) for d in diffs])
    Gf=gf[(idx[:,None]-idx[None,:])+M]; Ga=ga[(idx[:,None]-idx[None,:])+M]; Gp=Gf-Ga
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a); Vp=(kk*np.pi/(2*a))*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h; Wp=Vp*w[:,None]
    Qf=Wp.T@Gf@Wp; Qa=Wp.T@Ga@Wp; Qp=Wp.T@Gp@Wp
    for Q in (Qf,Qa,Qp): Q[:]=(Q+Q.T)/2
    return Qf,Qa,Qp
print(f"{'a':>5} {'lam_arch':>9} {'mu1-log a':>10} {'pair(1st)':>10} {'log a-mu1':>10} {'lam_full':>9} {'Dact':>8}")
for a in [0.4,0.7,1.0,1.5,2.0,2.5,3.0]:
    Qf,Qa,Qp=pieces(a)
    ea,Va=np.linalg.eigh(Qa); c=Va[:,0]; lam_arch=ea[0]/a
    pair=float(c@Qp@c)/float(c@c)/a
    lam_full=np.linalg.eigvalsh(Qf)[0]/a
    print(f"{a:5.1f} {lam_arch:+9.4f} {mu1-np.log(a):+10.4f} {pair:+10.4f} {np.log(a)-mu1:+10.4f} {lam_full:+9.4f} {lam_full-lam_arch:+8.4f}")
print("\nReconcile: does screw lam_arch = mu1-log a (col2)? Does pair track log a-mu1 (col4)? Is lam_full~0?")
