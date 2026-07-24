# Connes Step 2 (LEDGER §97), ANALYTIC handle: Davis-Kahan. If k_lambda is an approx eigenvector of the Weil/Suzuki
# operator A_lambda with residual r = (H-mu)c_k (mu=Rayleigh quotient), and 'gap' = dist(mu, rest of spectrum), then
#   sin(angle(k_lambda, theta_x)) <= ||r|| / gap      (Davis-Kahan sin-theta).
# So ||theta_x - k_lambda|| (unit) is CONTROLLED by the residual/gap -- an object that does NOT require resolving the
# ill-conditioned eigenvector (S97). If ||r||/gap is small + STABLE across resolution, that IS the analytic bound
# Connes needs (his Fig fpro1: spectra eps(lambda)~1-chi_2 already track). We compute it.
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto
sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
from k_lambda_construct import y2_matrix  # not needed but keeps import parity
from stage3_theta_vs_k import k_lambda_on

def weil_matrix_and_basis(a, M=8000, K=120):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gvals=np.array([screw_g(d,Lam,Nmax) for d in diffs])
    idx=np.arange(M+1); G=gvals[(idx[:,None]-idx[None,:])+M]
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg); Vv=np.sin(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    Wp=Vp*w[:,None]
    Qmat=Wp.T@G@Wp; Qmat=(Qmat+Qmat.T)/2
    H=Qmat/a                      # standard eigenproblem H c = lambda c (mass = a I, c orthonormal in eucl)
    return x,H,Vv,w

def project_k(x,Vv,w,lmbda):
    """coefficients of k_lambda in the sine basis: c_j = <sin_j|k>/<sin_j|sin_j>. Mass ||v||^2 = a*sum c^2."""
    k=k_lambda_on(x,lmbda)
    num=(Vv*w[:,None]).T@k                # <sin_j|k>
    den=(Vv*w[:,None]*Vv).sum(0)          # <sin_j|sin_j>
    c=num/den
    return c,k

print("DAVIS-KAHAN residual bound for Connes Step 2:  sin(angle) <= ||r||/gap.  (M=8000,K=120)")
print(f"{'lambda':>7} {'a':>6} {'eig0':>11} {'gap(eig1-eig0)':>14} {'||r||':>10} {'||r||/gap(bound)':>16} {'sin(meas)':>10}")
for lmbda in [1.8, 2.0, 2.2, 2.45, 2.7]:
    a=np.log(lmbda)
    x,H,Vv,w=weil_matrix_and_basis(a)
    ev,ec=np.linalg.eigh(H)
    c,k=project_k(x,Vv,w,lmbda)
    c=c/np.linalg.norm(c)
    mu=float(c@H@c)
    r=H@c-mu*c
    rn=np.linalg.norm(r)
    gap=ev[1]-ev[0]
    bound=rn/gap
    # measured sin(angle) between k (coeff c) and true ground eigenvector ec[:,0]
    ov=abs(c@ec[:,0]); sin_meas=np.sqrt(max(0.0,1-ov**2))
    print(f"{lmbda:>7.2f} {a:>6.3f} {ev[0]:>11.3e} {gap:>14.4e} {rn:>10.3e} {bound:>16.4e} {sin_meas:>10.4f}",flush=True)
print("\nREAD: if ||r||/gap (the ANALYTIC bound) is small and controls sin(meas), Davis-Kahan gives Connes Step 2.")
print("Key: is ||r||/gap STABLE/shrinking with lambda (unlike the ill-conditioned eigenvector distance S97)?")
