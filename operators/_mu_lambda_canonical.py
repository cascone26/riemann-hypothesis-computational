import numpy as np, math
import operators.connes_completeness_highT as C
# Canonical CCM/Connes Weil form tau = arch + boundary - prime (NO ad-hoc tuning).
# mu_lambda = min eigenvalue. Track SCALE-FREE ratio |mu|/||tau|| vs lambda=sqrt(LAMBDA_SQ),
# with N scaled to resolve. CCM Cor 3.8: mu_lambda -> 0 <=> RH. Honest: this is the crux quantity.

def mu_of(lsq, N):
    L = float(np.log(lsq)); VM = C.vm_table(lsq)
    tau, idx = C.build_tau(N, L, VM)
    ev = np.linalg.eigvalsh(tau)
    return ev[0], np.max(np.abs(ev)), np.median(np.abs(ev))

print("Canonical Weil-form min eigenvalue mu_lambda vs lambda (N scaled to resolve)")
print(f"{'LSQ':>5} {'lam':>6} {'N':>4} {'mu_min':>12} {'|mu|/||tau||':>12} {'|mu|/median':>12}")
rows=[]
for lsq in [14, 30, 60, 120, 250, 500]:
    lam = math.sqrt(lsq)
    N = int(20 + 6*math.log(lsq))   # scale dimension with prime reach
    mu, nrm, med = mu_of(lsq, N)
    r1 = abs(mu)/nrm; r2 = abs(mu)/med
    rows.append((lam, r1, r2))
    print(f"{lsq:>5} {lam:>6.2f} {N:>4} {mu:>12.4e} {r1:>12.4e} {r2:>12.4e}", flush=True)

lams=np.array([r[0] for r in rows]); r1s=np.array([r[1] for r in rows])
good = r1s>1e-14
if good.sum()>=3:
    p=np.polyfit(np.log(lams[good]), np.log(r1s[good]),1)[0]
    print(f"\nscale-free ratio |mu|/||tau|| ~ lambda^{p:.2f}  (<0 => collapsing toward 0 = Cor 3.8 / RH-consistent)")
