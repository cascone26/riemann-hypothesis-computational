#!/usr/bin/env python3
"""RH external-frontier watcher. Checks arXiv for new papers that could supply the ONE open piece
our verified forge chain needs (Weil positivity / even-simple / scaling-site Hodge index), and pings
cobo-cascone-rh. The forge master theorem closes RH the instant this input is proven — so catch it fast.
Run daily. State in seen.json (no re-ping)."""
import json, os, sys, urllib.request, urllib.parse, re, html
HERE=os.path.dirname(os.path.abspath(__file__)); SEEN=os.path.join(HERE,'seen.json')
# authors + keyword queries most likely to carry the missing piece
QUERIES=[
 'au:Connes_A AND (abs:zeta OR abs:Weil OR abs:prolate OR abs:Riemann)',
 'au:Consani_C AND (abs:zeta OR abs:Weil OR abs:Riemann)',
 'abs:"Weil positivity" AND abs:Riemann',
 'abs:"zeta cycles"',
 'abs:prolate AND abs:"Riemann"',
 'abs:"scaling site" AND abs:Riemann',
]
def fetch(q):
    url='http://export.arxiv.org/api/query?'+urllib.parse.urlencode({'search_query':q,'sortBy':'submittedDate','sortOrder':'descending','max_results':8})
    try:
        with urllib.request.urlopen(url,timeout=30) as r: x=r.read().decode()
    except Exception as e: return []
    out=[]
    for m in re.finditer(r'<entry>(.*?)</entry>', x, re.S):
        e=m.group(1)
        idm=re.search(r'<id>(.*?)</id>',e); tim=re.search(r'<title>(.*?)</title>',e,re.S); dm=re.search(r'<published>(.*?)</published>',e)
        if idm and tim:
            out.append({'id':idm.group(1).strip(),'title':html.unescape(re.sub(r'\s+',' ',tim.group(1)).strip()),'date':(dm.group(1)[:10] if dm else '')})
    return out
def main():
    seen=set()
    if os.path.exists(SEEN):
        try: seen=set(json.load(open(SEEN)))
        except: seen=set()
    new=[]
    for q in QUERIES:
        for e in fetch(q):
            if e['id'] not in seen:
                seen.add(e['id']); new.append(e)
    json.dump(sorted(seen), open(SEEN,'w'))
    # first run: baseline only (record, don't ping the backlog)
    baseline = not os.path.exists(SEEN+'.init')
    if baseline:
        open(SEEN+'.init','w').write('1'); print(f'baseline recorded: {len(seen)} papers seen, no ping'); return
    if new:
        msg='RH WATCH: %d new arXiv paper(s) on the open frontier (Weil positivity / even-simple / scaling-site). Check if any supplies the forge chain\'s missing input:\\n'%len(new)+'\\n'.join(f"- {e['date']} {e['title']} ({e['id']})" for e in new[:6])
        os.system('node ~/tools/notify/notify.js '+repr(msg)+' cobo-cascone-rh 2>/dev/null')
        print(f'PINGED: {len(new)} new')
    else: print('no new frontier papers')
if __name__=='__main__': main()
