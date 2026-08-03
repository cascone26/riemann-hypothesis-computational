"""Non-abelian probe: does the elementary readout on the eta(z)eta(23z) zeros read the
S3 Artin representation — distinguishing principal from non-principal SPLIT primes?"""
import numpy as np
from scipy.special import loggamma as slg
def load_pari(f):
    return np.array([float(l.replace(" E","e")) for l in open(f) if l.strip()])
zeros = np.sort(load_pari("h23_zeros.txt"))
ph = load_pari("h23_argLp.txt")
N = len(zeros)
# gamma factor weight-1 level 23: Gamma_C(s) with conductor 23: theta = Im lgG(1/2+it... shift 0? weight k=1: Gamma_C(s + (k-1)/2) = Gamma_C(s)
th = np.imag(slg(0.5+1j*zeros)) - zeros*np.log(2*np.pi/np.sqrt(23))
k = np.arange(1, N+1) - 0.5
S = k - th/np.pi
print("N=%d | S-band [%.2f,%.2f] drift %.3f" % (N, S.min(), S.max(), abs(S[-8:].mean()-S[:8].mean())))
d = ph - np.pi*S
z = np.exp(1j*d).mean()
sd = np.std(np.mod(d-np.angle(z)+np.pi, 2*np.pi)-np.pi)
print("position-only R=%.6f sd=%.5f c=%+.5f" % (abs(z), sd, np.angle(z)))
w = np.exp(-1j*ph)
G2 = (w*np.exp(1j*np.log(2)*zeros)).mean()
rot = np.exp(-1j*(np.angle(G2)))  # b(2)=+1 anchor (non-principal)
noise = np.mean([abs((w*np.exp(1j*u*zeros)).mean()) for u in (0.95,1.5,2.2,2.9)])
# channels: (n, b_true, tag)
tests = [(2,+1,"nonpr"), (3,+1,"nonpr"), (13,+1,"nonpr"), (4,+1,"p^2 npr"), (9,+1,"p^2 npr"),
         (6,+1,"npr*npr"), (5,0,"inert"), (7,0,"inert"), (25,-1,"inert p^2"), (49,-1,"inert p^2"),
         (23,-1,"ram"), (59,-2,"PRINCIPAL"), (101,-2,"PRINCIPAL")]
print("readout (noise %.4f), anchor b(2)=+1:" % noise)
ok = 0
for n, b, tag in tests:
    G = ((w*np.exp(1j*np.log(n)*zeros)).mean()*rot).real
    good = (b == 0 and abs(G) < 2.5*noise) or (b != 0 and np.sign(G) == np.sign(b))
    ok += bool(good)
    print(f"  n={n:<4} b={b:+d} ({tag:<10}) G={G:+.4f}  {'anchor' if n==2 else ('OK' if good else 'MISS')}")
print("total %d/%d" % (ok, len(tests)))
# THE decisive ratio: principal (b=-2) vs nonprincipal (b=+1): sign flip within SPLIT class
