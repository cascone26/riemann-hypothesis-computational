#!/usr/bin/env python3
"""§401 (angle 5, CORRECTED): the GENUINE F_q Weil form uses the ROSATI involution alpha -> q/alpha
(functional-equation partner), NOT complex conjugation. Trace/Rosati form:
   A[r,s] = sum_i alpha_i^r * (q/alpha_i)^s
On RH (|alpha_i|^2 = q): q/alpha_i = conj(alpha_i) => A is a Hermitian PSD Gram (Weil positivity).
OFF-line (|alpha_i|^2 != q): q/alpha_i != conj(alpha_i) => A[r,s] = q^s * sum_i alpha_i^{r-s} = q^s p_{r-s}
(power sums) — a HANKEL form that becomes INDEFINITE. RH <=> Rosati form PSD <=> |alpha|=sqrt q.
This is the correct indefinite certificate (Castelnuovo-Severi / Hodge-index in eigenvalue form).
Test: signature of (A+A^T)/2 on-line vs off-line, genus 2 (2g=4 eigenvalues, conj+FE pairs)."""
import numpy as np
q = 25; g = 2; sq = np.sqrt(q)
th1, th2 = 0.7, 2.1
def eigs(b1=0.0, b2=0.0):
    # 2g=4 Frobenius eigenvalues; FE pairs (alpha, q/alpha). b = radial off-line (0=RH).
    a1 = sq*np.exp(b1)*np.exp(1j*th1)
    a2 = sq*np.exp(b2)*np.exp(1j*th2)
    return [a1, np.conj(q/np.conj(a1)) if False else q/a1*np.exp(0),  # keep it simple: {a1, conj(a1) scaled}
            a2, q/a2]
def rosati(al, M):
    al = np.array(al)
    A = np.zeros((M, M), dtype=complex)
    for r in range(M):
        for s in range(M):
            A[r, s] = np.sum(al**r * (q/al)**s)
    return np.real((A + A.conj().T)/2)
def inertia(Ms, tol=1e-9):
    ev = np.linalg.eigvalsh(Ms); s = tol*max(1, abs(ev).max())
    return int((ev>s).sum()), int((ev<-s).sum()), int((abs(ev)<=s).sum())
# eigenvalues: on RH must be closed under conj AND under q/.; use two conj-pairs on circle radius sqrt q,
# displaced radially by b (breaking |alpha|=sqrt q while keeping FE partner q/alpha in the set).
def eigset(b):
    a1 = sq*np.exp(b)*np.exp(1j*th1)   # |a1|=sqrt q * e^b
    a2 = sq*np.exp(b)*np.exp(-1j*th1)  # conj direction, same radius
    # FE partners q/a1, q/a2 have radius sqrt q * e^{-b}; full set = {a1,a2,q/a1,q/a2} (RH: b=0 => all on circle)
    return [a1, a2, q/a1, q/a2]
print(f"GENUS-2 F_q={q} ROSATI/Weil form (alpha->q/alpha). RH <=> PSD. off-line b breaks |alpha|=sqrt q.")
print(f"{'b':>6} {'M':>3} {'inertia A':>12} {'lmin':>11} {'lmax':>10}")
for b in (0.0, 0.02, 0.05, 0.1, 0.2):
    for M in (5, 7):
        A = rosati(eigset(b), M)
        ev = np.linalg.eigvalsh(A)
        print(f"{b:6.3f} {M:3d} {str(inertia(A)):>12} {ev[0]:+11.3e} {ev[-1]:+10.2e}", flush=True)
print("\nreading: b=0 (RH) => A PSD (inertia (rank,0,.)). b>0 (off circle) => A develops NEGATIVE eigenvalues")
print("(indefinite). # negative per displaced pair = the Hodge-index signature deficit. THIS is the genuine")
print("F_q certificate (Rosati positivity = RH); compare its structure to the number-field pole+arch-prime.")
