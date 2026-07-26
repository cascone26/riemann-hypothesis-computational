#!/usr/bin/env python3
"""D1 — Dual/SOS certificate extraction (the Viazovska bet). The certified Groskin tau (PD) has an SOS
certificate = its LDL^T: tau = L^T D L, so tau(v) = Σ_i d_i (L v)_i^2 — a manifest sum of squares. We mined
the pivots D (§211, the deep well) but NEVER the L factor (the SOS linear forms). QUESTION: does the
certificate structure (pivots D + linear forms L, esp. the CRITICAL near-null one) converge to a RECOGNIZABLE
form as c→∞? If yes → conjecture the limit certificate (the 'magic function' analog) → verify positivity as an
IDENTITY (escapes the product-emergent cancellation). If structureless/divergent → honest new fact.

Structure-mining uses float64 midpoints of the Arb entries (the L-factor + large pivots are O(1), reliable;
only the single tiniest near-null pivot is float-noisy — irrelevant for STRUCTURE).
"""
import importlib.util, os
import numpy as np
GPATH=os.path.expanduser('~/projects/riemann/external/groskin_2607.02828/code/arb_ldlt_certify.py')
spec=importlib.util.spec_from_file_location('gc',GPATH); g=importlib.util.module_from_spec(spec); spec.loader.exec_module(g)
def tof(x):
    try: return float(x)
    except Exception:
        s=str(x);
        if '+/-' in s: s=s.split('+/-')[0].strip('[] ')
        return float(s)

def ldl(M):
    n=M.shape[0]; d=np.zeros(n); L=np.eye(n)
    for i in range(n):
        s=M[i,i]-sum(L[i,k]**2*d[k] for k in range(i)); d[i]=s
        for j in range(i+1,n):
            L[j,i]=(M[j,i]-sum(L[j,k]*L[i,k]*d[k] for k in range(i)))/d[i] if d[i]!=0 else 0
    return d,L

N=20; prec=1200
results={}
for c in [13,29,50,80]:
    A,DIM=g.build_arb_tau(c,N,prec)
    M=np.array([[tof(A[i,j]) for j in range(DIM)] for i in range(DIM)])
    d,L=ldl(M)
    # critical near-null: index of min |pivot| (the marginal SOS term)
    kmin=int(np.argmin(np.abs(d)))
    # the critical linear form = row kmin of the SOS map (sqrt(d) * L-column influence). Use L[:,kmin]?
    # tau(v)=Σ d_i (Σ_j L[i,j] v_j)^2. SOS term i = sqrt(d_i)*(L row i). Critical term i=kmin (smallest d).
    crit = L[kmin,:].copy()  # the linear form coefficients (in Fourier index n=-N..N)
    results[c]=dict(d=d, L=L, kmin=kmin, crit=crit, DIM=DIM)
    logmin=np.log(np.abs(d[kmin]))
    print(f"c={c:3d} DIM={DIM}: min-pivot at index {kmin} (freq n={kmin-N:+d}), log|d_min|={logmin:.3f}")

print("\n=== PIVOT WELL: log|d_i| vs frequency n, across c (does the well SHAPE stabilize?) ===")
print(f"{'n':>4}", *[f"c={c:>4d}" for c in results])
for i in range(0,N*2+1,3):
    n=i-N
    row=[f"{np.log(np.abs(results[c]['d'][i])):+7.2f}" if i<results[c]['DIM'] else "   -   " for c in results]
    print(f"{n:>4}", *row)

print("\n=== CRITICAL near-null SOS LINEAR FORM: coefficients (freq n) of the marginal square, per c ===")
print("(if this converges to a recognizable shape as c grows => the extremal witness / magic-function candidate)")
for c in results:
    r=results[c]; cr=r['crit']; cr=cr/np.max(np.abs(cr))  # normalize
    dom=np.argsort(np.abs(cr))[::-1][:6]
    print(f" c={c:3d} kmin-freq n={r['kmin']-N:+d}: dominant coeffs", ", ".join(f"n={i-r['DIM']//2:+d}:{cr[i]:+.3f}" for i in sorted(dom)))

# HEIGHT-SPACE convergence: index n -> height t = 2*pi*n/log(c). Resample each critical form as f(t) on a
# common t-grid, then compare shapes. This corrects the index-shift (peak tracks a fixed height as c grows).
import numpy as np
print("\n=== HEIGHT-SPACE analysis: critical form peak vs first-zero prediction n_peak = gamma1*log(c)/(2pi) ===")
g1=14.134725
for c in results:
    r=results[c]; cr=r['crit']; N_=r['DIM']//2
    n_peak = int(np.argmax(np.abs(cr))) - N_
    pred = g1*np.log(c)/(2*np.pi)
    print(f" c={c:3d}: observed |n_peak|={abs(n_peak):2d}  vs  gamma1-prediction {pred:.2f}  (t at peak = {abs(n_peak)*2*np.pi/np.log(c):.2f})")

print("\n=== SHAPE CONVERGENCE in HEIGHT space (resample f(t), common grid t in [0,30], cosine sim) ===")
tgrid=np.linspace(0.5,30,300)
def height_form(c,r):
    cr=r['crit']; N_=r['DIM']//2; L=np.log(c)
    # magnitudes at heights t_n = 2*pi*n/L for n=1..N_ (positive side; form is real, use |.| envelope)
    tn=np.array([2*np.pi*n/L for n in range(1,N_+1)])
    val=np.array([cr[N_+n] for n in range(1,N_+1)])  # coeff at +n
    valm=np.array([cr[N_-n] for n in range(1,N_+1)])  # coeff at -n
    env=np.abs(val)+np.abs(valm)
    # interpolate envelope onto tgrid
    return np.interp(tgrid, tn, env, left=0, right=0)
cs=sorted(results)
forms={c:height_form(c,results[c]) for c in cs}
for a,b in zip(cs,cs[1:]):
    fa=forms[a]/ (np.linalg.norm(forms[a])+1e-30); fb=forms[b]/(np.linalg.norm(forms[b])+1e-30)
    print(f" c={a} vs c={b}: |cosine sim| in HEIGHT space = {abs(float(fa@fb)):.4f}  (→1 = converging)")
print("\n height-space critical-form envelope (t, then c=13/29/50/80 values, normalized):")
for i in range(0,300,25):
    t=tgrid[i]; print(f"  t={t:5.1f}", *[f"{forms[c][i]/np.max(forms[c]):.3f}" for c in cs])
print("\nREAD: if HEIGHT-space sim→1 and the envelope peaks at the zeros' heights (14.1, 21.0, 25.0...) => the")
print("critical SOS form is the zero-detecting extremal witness with a STABLE limit => the guessable magic object.")
