#!/usr/bin/env python3
"""§431 (core falsification, extended instrument): directly test 'zero-side PSD ⟺ zeros real'. Take real zeros,
move ONE to ρ=1/2+β+iγ_k (with its FE quadruple {1/2±β±iγ}), recompute the zero-side 'Gram' with F at complex
argument, check when λmin goes NEGATIVE (PSD breaks). Measure the breaking displacement β* vs the zero's height γ.
This validates the mechanism AND tests §422: does breaking get EASIER (smaller β*) or HARDER at high γ?"""
import numpy as np, l2_fejer as LF
Z=np.load('zeros150.npy'); W=LF.W
def fej_c(x):  # (sin(Wx)/(Wx))^2 at complex x
    y=W*x; return (np.sin(y)/y)**2
def Frow_c(c,r):  # F_c at complex r
    return fej_c(r-c)+fej_c(r+c)
def zeroside(cs, kmove, beta):
    """zero-side Gram: on-line zeros contribute 2F(γ)F(γ); the moved zero contributes its FE quadruple."""
    n=len(cs); M=np.zeros((n,n))
    for i,g in enumerate(Z):
        if i==kmove: continue
        Fv=np.array([Frow_c(c,g).real for c in cs])
        M+=2.0*np.outer(Fv,Fv)
    # moved zero ρ=1/2+β+iγ: explicit-formula argument (ρ-1/2)/i = γ - iβ ; FE quadruple gives ±β,±γ
    g=Z[kmove]
    for arg in (g-1j*beta, g+1j*beta):     # the two FE-partners (±β); conj handled by real part
        Fv=np.array([Frow_c(c,arg) for c in cs])
        M+=np.real(np.outer(Fv,Fv))
    return 0.5*(M+M.T)
BAND=(14.0,44.0); cs=LF.probes(BAND,4.0)
print(f"§431 off-line break: move one zero to 1/2+β+iγ, find β* where zero-side λmin<0. band{BAND} n={len(cs)}")
print(f"{'zero γ':>9} {'k':>3}  β* (breaking displacement)   λmin(β=0)  λmin(β=0.3)")
for k in [0,3,7,12,20,30,50,80,120,148]:
    if k>=len(Z): continue
    g=Z[k]
    l0=float(np.linalg.eigvalsh(zeroside(cs,k,0.0))[0])
    # bisect for β* where λmin crosses 0
    lo,hi=0.0,3.0
    if float(np.linalg.eigvalsh(zeroside(cs,k,hi))[0])>=0: bstar=np.inf
    else:
        for _ in range(30):
            mid=(lo+hi)/2
            if float(np.linalg.eigvalsh(zeroside(cs,k,mid))[0])>=0: lo=mid
            else: hi=mid
        bstar=(lo+hi)/2
    l3=float(np.linalg.eigvalsh(zeroside(cs,k,0.3))[0])
    print(f"{g:9.3f} {k:3d}   β*={bstar:8.4f}              {l0:+.3e}  {l3:+.3e}",flush=True)
print("""
reading: λmin(β=0)≥0 (real zeros => PSD, validates mechanism). β* = displacement that breaks PSD. If β* DECREASES
with γ (high zeros break with smaller displacement), off-line detection is HARDER at high γ in the sense that the
signal per unit β is... actually smaller β* = MORE sensitive. Compare to §422 (high γ off-line signal weak): expect
β* to GROW with γ (need larger displacement to break) if the instrument's sensitivity to high-γ off-line drops.
""")
