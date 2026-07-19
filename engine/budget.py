#!/usr/bin/env python3
"""USAGE-SAFETY layer — the guardrails so the machine can NEVER blow Jacob's weekly Max quota.
Core principle: RUNNING methods is FREE (local numpy/mpmath, 24/7). Only CREATIVE GENERATION costs.
So we run free forever and spend Claude RARELY, HARD-CAPPED, and predictably. A broken loop or a
runaway can never eat the week."""
import json, os, time, datetime, subprocess
ST=os.path.join(os.path.dirname(__file__),"budget_state.json")

# ADAPTIVE limits — daily Claude allowance scales with AVAILABLE HEADROOM (protect the week, but
# use plenty when it's safe). Free models (Gemini/Groq) handle unlimited breadth at zero cost regardless.
MIN_PAID_INTERVAL_S  = 45*60      # >=45min between paid calls (was 4h — free does the volume)
CIRCUIT_BREAK_FAILS  = 4          # consecutive bad generations -> paid OFF for the day
WEEK_USAGE_CEILING   = 60         # >= this all-usage% -> FREE-ONLY (protect rest of week)
def paid_cap_today():
    """Daily Claude allowance = f(current week usage%). Generous when safe, zero when tight."""
    u=_week_usage_pct()
    if u is None:    return 6      # unknown -> moderate
    if u < 15:       return 20     # early/low week: explore hard
    if u < 30:       return 12
    if u < 45:       return 6
    if u < WEEK_USAGE_CEILING: return 2
    return 0                       # protect the rest of the week

def _load():
    try: return json.load(open(ST))
    except: return {"day":"","paid_today":0,"last_paid":0,"consec_fail":0}
def _save(s): json.dump(s,open(ST,"w"),indent=2)
def _today(): return datetime.date.today().isoformat()  # note: engine runs as daemon, real date OK here

def _week_usage_pct():
    """Best-effort read of current all-usage% from the usage tracker; None if unavailable."""
    try:
        r=subprocess.run(["node",os.path.expanduser("~/tools/usage-tracker.js"),"status"],
                         capture_output=True,text=True,timeout=20)
        import re
        m=re.search(r'(\d+(?:\.\d+)?)\s*%',r.stdout)
        return float(m.group(1)) if m else None
    except Exception: return None

def can_spend_paid():
    """Gate for a PAID (Claude) call. Returns (ok:bool, reason:str)."""
    s=_load()
    if s["day"]!=_today(): s={"day":_today(),"paid_today":0,"last_paid":s.get("last_paid",0),"consec_fail":0}; _save(s)
    cap=paid_cap_today()
    if s["paid_today"]>=cap: return False,f"adaptive daily cap {cap} reached (usage-aware)"
    if time.time()-s.get("last_paid",0) < MIN_PAID_INTERVAL_S:
        return False,f"min {MIN_PAID_INTERVAL_S//3600}h interval not elapsed"
    if s.get("consec_fail",0)>=CIRCUIT_BREAK_FAILS: return False,"circuit-breaker OPEN (too many bad generations today)"
    u=_week_usage_pct()
    if u is not None and u>=WEEK_USAGE_CEILING: return False,f"week usage {u}% >= ceiling {WEEK_USAGE_CEILING}% -> FREE-ONLY"
    return True,"ok"

def record_paid(success):
    s=_load()
    if s["day"]!=_today(): s={"day":_today(),"paid_today":0,"last_paid":0,"consec_fail":0}
    s["paid_today"]=s.get("paid_today",0)+1; s["last_paid"]=time.time()
    s["consec_fail"]=0 if success else s.get("consec_fail",0)+1
    _save(s)

if __name__=="__main__":
    ok,why=can_spend_paid(); print(f"can_spend_paid: {ok} ({why})"); print("state:",json.load(open(ST)) if os.path.exists(ST) else "none")
