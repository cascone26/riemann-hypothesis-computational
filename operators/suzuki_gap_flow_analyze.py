#!/usr/bin/env python3
"""suzuki_gap_flow_analyze.py -- §241 analysis: plateaus (stickiness), handoffs, theta-robustness.
Reads results/suzuki_gap_flow.json. For each theta:
 (1) PLATEAU census: points where |dz_k/da| < 0.3 * (z/a) [= 30% of free-scaling speed] sustained over
     >=2 consecutive steps -> plateau position = mean z there. Compare to gamma_k.
 (2) HANDOFF census: local minima of mingap(a) -> avoided crossings.
 (3) theta-robustness: plateau positions across theta values (should agree at the gammas if stickiness
     is a property of the ZEROS, not of the extension phase).
"""
import json, numpy as np
GAMMAS=[14.1347,21.0220,25.0109,30.4249,32.9351,37.5862,40.9187,43.3271,48.0052,49.7738,52.9703,56.4462]
d=json.load(open("results/suzuki_gap_flow.json"))
all_plateaus={}
for th,rows in d.items():
    if not rows: continue
    avals=[r["a"] for r in rows]
    # build branch matrix by continuity
    branches=[]
    for r in rows:
        zs=r["zeros"]
        if not branches:
            branches=[[z] for z in zs]; continue
        last=[next((v for v in reversed(b) if v is not None),1e9) for b in branches]
        used=set()
        for b in branches: b.append(None)
        for z in sorted(zs):
            cand=sorted(range(len(last)),key=lambda j:abs(z-last[j]) if j not in used else 9e9)
            j=cand[0]
            if abs(z-last[j])<1.5 and j not in used: branches[j][-1]=z; used.add(j)
            else: branches.append([None]*(len(branches[0])-1)+[z])
    # pad
    L=max(len(b) for b in branches)
    for b in branches: b+= [None]*(L-len(b))
    plateaus=[]
    for bi,b in enumerate(branches):
        run=[]
        for i in range(1,len(avals)):
            if b[i] is None or b[i-1] is None:
                run=[]; continue
            da=avals[i]-avals[i-1]
            if da<=0: continue
            v=(b[i]-b[i-1])/da; free=b[i]/avals[i]
            if abs(v)<0.3*free: run.append((avals[i],b[i]))
            else:
                if len(run)>=2: plateaus.append((run[0][0],run[-1][0],float(np.mean([z for _,z in run]))))
                run=[]
        if len(run)>=2: plateaus.append((run[0][0],run[-1][0],float(np.mean([z for _,z in run]))))
    all_plateaus[th]=plateaus
    mg=[(r["a"],r["mingap"]) for r in rows if r["mingap"]==r["mingap"]]
    dips=[(mg[i][0],round(mg[i][1],2)) for i in range(1,len(mg)-1) if mg[i][1]<mg[i-1][1] and mg[i][1]<mg[i+1][1] and mg[i][1]<2.6]
    print(f"\n=== theta={th} ===")
    print(f" plateaus (a_start,a_end,z): {[(round(a,2),round(b,2),round(z,2)) for a,b,z in sorted(plateaus,key=lambda p:p[2])]}")
    for a0,a1,z in sorted(plateaus,key=lambda p:p[2]):
        g=min(GAMMAS,key=lambda g:abs(g-z))
        print(f"   plateau z={z:6.2f}  [a {a0:.2f}-{a1:.2f}]  nearest gamma={g:.4f}  offset={z-g:+.3f}")
    print(f" avoided-crossing dips (a,mingap): {dips}")
print("\n=== THETA-ROBUSTNESS: plateau positions by theta ===")
for th,pl in all_plateaus.items():
    print(f" theta={th}: {sorted(round(z,2) for _,_,z in pl)}")
print("\nREAD: if plateau positions agree across theta AND sit at gammas (offsets shrinking with a),")
print("stickiness is a property of the zeros -> the staircase frame is real; log §241.")
