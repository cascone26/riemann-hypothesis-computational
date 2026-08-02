import numpy as np, os, sys

# WEIL FORM — CORRECT build via the screw kernel (reuse suzuki machinery), NEW analysis:
# full spectrum + eigenvectors of the Galerkin matrix Q_{kl}=int int g(x-y)v_k'(x)v_l'(y).
# Q PSD <=> RH (Weil/Suzuki). Questions:
#  (1) full eigenvalue spectrum (not just min): PSD? where does the near-null sit?
#  (2) the near-null eigenvector as a function on [-a,a]: does it have NODES at points related
#      to the zeros / is it the "zero-node function" (474)?
#  (3) razor scaling: min-eig vs a and K (saturation).

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "suzuki_screw_lambda_a.py")).read()
cut = src.find("if __name__")
exec(src[:cut] if cut > 0 else src, globals())

def Qmat_and_grid(a, M=3000, K=40):
    Nmax = int(np.floor(np.exp(2 * a))) + 1
    Lam = vonmangoldt_upto(min(Nmax, 30_000_000))
    x = np.linspace(-a, a, M + 1); h = x[1] - x[0]
    diffs = np.arange(-M, M + 1) * h
    gvals = np.array([screw_g(d, Lam, Nmax) for d in diffs])
    idx = np.arange(M + 1)
    G = gvals[(idx[:, None] - idx[None, :]) + M]
    kk = np.arange(1, K + 1)
    arg = np.pi * np.outer((x + a), kk) / (2 * a)
    Vp = (kk * np.pi / (2 * a)) * np.cos(arg)
    w = np.full(M + 1, h); w[0] = w[-1] = 0.5 * h
    Wp = Vp * w[:, None]
    Q = Wp.T @ G @ Wp; Q = (Q + Q.T) / 2
    return Q / a, x, kk

def sine_recon(vec, x, kk, a):
    # reconstruct the test function u(x)=sum vec_k sin(k pi (x+a)/2a) on the grid
    B = np.sin(np.pi * np.outer((x + a), kk) / (2 * a))
    return B @ vec

print("=== (1) full spectrum of the Weil Galerkin form (PSD <=> RH) ===")
for a in [1.5, 2.5, 3.5]:
    Q, x, kk = Qmat_and_grid(a, K=45)
    ev, V = np.linalg.eigh(Q)
    print(f"a={a}: min={ev[0]:+.4e} 2nd={ev[1]:+.4e} max={ev[-1]:.3e} | n_neg={np.sum(ev< -1e-8*ev[-1])}/{len(ev)}", flush=True)

print("\n=== (2) near-null eigenvector structure at a=3.0 (nodes vs primes/scale) ===")
Q, x, kk = Qmat_and_grid(3.0, K=50)
ev, V = np.linalg.eigh(Q)
u0 = sine_recon(V[:, 0], x, kk, 3.0)   # lowest-eigenvalue test function
# node locations (sign changes) of u0(x)
sgn = np.sign(u0); nodes = x[:-1][sgn[:-1] * sgn[1:] < 0]
print(f"lowest-eig={ev[0]:+.3e}; near-null test-fn has {len(nodes)} nodes on [-3,3]")
print(f"node positions: {np.round(nodes,3)}")
print(f"e^|node| (prime-power scale?): {np.round(np.exp(np.abs(nodes)),2)}")

print("\n=== (3) razor scaling min-eig vs (a, K) ===")
for a in [2.0, 3.0]:
    for K in [30, 50, 70]:
        Q, x, kk = Qmat_and_grid(a, K=K)
        ev = np.linalg.eigvalsh(Q)
        print(f"  a={a} K={K}: min-eig={ev[0]:+.4e} (min/max={ev[0]/ev[-1]:+.2e})", flush=True)
print("DONE", flush=True)
