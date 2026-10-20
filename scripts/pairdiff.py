#!/usr/bin/env python3
"""pairdiff.py labexe "<SFW or empty>" "<SFW2>" -- files... : per-file kpi for two weight sets; prints worst/best ratios + summary"""
import sys,os,subprocess,glob,concurrent.futures as cf,json
exe=os.path.abspath(sys.argv[1]); w1=sys.argv[2]; w2=sys.argv[3]; files=sorted(sum([glob.glob(g) for g in sys.argv[sys.argv.index('--')+1:]],[]))
def run(w):
    env=dict(os.environ)
    if w: env['SFW']=w
    chunks=[files[j::12] for j in range(12)]
    def one(ch):
        out=subprocess.run([exe]+ch,capture_output=True,text=True,env=env).stdout
        return {l.split()[0]:float(l.split()[2][4:]) for l in out.splitlines() if ' w=' in l and 'kpi=' in l}
    r={}
    with cf.ThreadPoolExecutor(12) as ex:
        for x in ex.map(one,chunks): r.update(x)
    return r
a=run(w1); b=run(w2)
rows=sorted(((b[f]-a[f])/max(a[f],1e-9),f,a[f],b[f]) for f in files if f in a and f in b)
print('n',len(rows),'mean A',sum(a.values())/len(a),'mean B',sum(b.values())/len(b),'oracle',sum(max(a[f],b[f]) for f in a)/len(a))
print('B worse than A (worst 8):'); [print(' %.3f %s A=%.4f B=%.4f'%r) for r in rows[:8]]
print('B better than A (best 8):'); [print(' %.3f %s A=%.4f B=%.4f'%r) for r in rows[-8:]]
json.dump(dict(A=a,B=b),open('build/pairdiff_last.json','w'))
