#!/usr/bin/env python3
"""Analyze mE_overnight.txt: filter by conditioning (Gpivmin healthy AND nprobes≈capacity so neither
over- nor under-resolved), fit m_E(W) and test the m_E·D law. Prints the well-conditioned window and
the verdict on m_E·D (const / rising / other)."""
import sys, numpy as np
fn = sys.argv[1] if len(sys.argv)>1 else "mE_overnight.txt"
rows=[]
for L in open(fn):
    if L.startswith('#') or not L.strip(): continue
    p=L.split()
    if len(p)<7 or 'DEGENERATE' in L or 'ERROR' in L: continue
    try: rows.append((float(p[0]),float(p[1]),int(p[2]),float(p[3]),float(p[4]),float(p[5]),float(p[6])))
    except: continue
if not rows: print("no valid rows yet"); sys.exit()
print(f"{len(rows)} points. Columns: 4W cap nprobes Gpiv m_E D mE·D  | resolution flag")
print("4W    cap   npr  Gpiv     m_E      D      mE·D    | flag")
well=[]
for W,cap,npr,gp,mE,D,mED in rows:
    # trust Gpivmin (real conditioning signal; §808 proved 8@cap6.31 is fine w/ healthy pivot)
    under = npr < cap - 0.5   # genuinely under-resolved (fewer probes than Shannon #)
    illc = gp < 0.05          # ill-conditioned G (over-capacity tell)
    flag = "UNDER-CAP" if under else ("ill-cond(over-cap)" if illc else "GOOD")
    if flag=="GOOD": well.append((W,mE,D,mED))
    print(f"{W:.2f}  {cap:.2f}  {npr:2d}  {gp:.2e} {mE:.5f} {D:.3f}  {mED:.5f}  | {flag}")
if len(well)>=3:
    Ws=np.array([w[0] for w in well]); mEs=np.array([w[1] for w in well]); mEDs=np.array([w[3] for w in well])
    print(f"\nWELL-CONDITIONED points: {len(well)} (4W {Ws.min():.2f}-{Ws.max():.2f})")
    cv=np.std(mEDs)/np.mean(mEDs)
    print(f"m_E·D over good points: mean={np.mean(mEDs):.4f} std={np.std(mEDs):.4f} CV={cv:.3f}")
    # slope of m_E·D vs W
    if len(well)>=2:
        sl=np.polyfit(Ws,mEDs,1)[0]
        print(f"m_E·D slope vs 4W = {sl:+.4f}  → {'~CONST (CV<10%)' if cv<0.10 else ('RISING' if sl>0.005 else 'varying')}")
        slm=np.polyfit(Ws,np.log(mEs),1)[0]
        print(f"m_E log-slope = {slm:.3f} (m_E ~ e^{{{slm:.2f}·4W}})")
else:
    print(f"\nonly {len(well)} well-conditioned points so far — need ≥3 for the law verdict")
