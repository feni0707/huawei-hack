#!/usr/bin/env python3
# simcmp.py A B list...  — парное сравнение двух стендов scripts/simc.cpp по классам (класс = имя папки теста).
#   A, B: "build/simcA" или "build/simcA|1=0.5 2=3" (ручки K); list — .txt со списком тестов или сами .in.
#   Печатает классы, где есть разница: среднее B−A, ст. ошибка, победы/поражения; open — взвешенный Σ по materials/tests.
#   Результаты кешируются по md5(бинарник + ручки) в build/simcmp_cache. J=10 — потоков. SFCUTS=20 — сглаживание.
import subprocess, sys, re, os, collections, math, json, hashlib
from concurrent.futures import ThreadPoolExecutor
CACHE='build/simcmp_cache'
os.makedirs(CACHE,exist_ok=True)
files=[]
for a in sys.argv[3:]:
    if a.endswith('.txt'): files+= [l.strip() for l in open(a) if l.strip()]
    else: files.append(a)
J=int(os.environ.get('J','10'))
def run(spec, fl):
    b,_,k=spec.partition('|'); env=dict(os.environ,K=k)
    out=subprocess.run([b]+fl,env=env,capture_output=True,text=True).stdout
    return {m.group(1):(float(m.group(2)), 'WA' in m.group(0)) for m in re.finditer(r'(\S+)\s+w=\S+\s+kpi=(\S+).*',out)}
def runall(spec):
    b,_,k=spec.partition('|')
    h=hashlib.md5((open(b,'rb').read()+k.encode()+os.environ.get('SFCUTS','').encode())).hexdigest()
    cf=os.path.join(CACHE,h+'.json'); d={}
    if os.path.exists(cf): d=json.load(open(cf))
    todo=[f for f in files if f not in d]
    if todo:
        parts=[todo[i::J*4] for i in range(J*4)]
        with ThreadPoolExecutor(J) as ex: rs=list(ex.map(lambda p: run(spec,p) if p else {}, parts))
        for r in rs: d.update(r)
        json.dump(d,open(cf,'w'))
    return d
a=runall(sys.argv[1]); b=runall(sys.argv[2])
c=collections.defaultdict(list); wa=0
for f in files:
    if b[f][1]: wa+=1; print("WA",f)
    if f.startswith('materials/'): cls='OPEN'+f.split('/')[-1]
    else: cls=f.split('/')[-2]
    c[cls].append(b[f][0]-a[f][0])
W={1:0.1,15:10,12:4,13:4,14:4}
oa=ob=0
for f in files:
    if f.startswith('materials/tests/'):
        t=int(f.split('/')[-1][:-3]); w=W.get(t,1); oa+=w*a[f][0]; ob+=w*b[f][0]
tot=0;nn=0
for k,v in sorted(c.items()):
    n=len(v); m=sum(v)/n; sd=math.sqrt(sum((x-m)**2 for x in v)/max(1,n-1))
    tot+=m; nn+=1
    if k.startswith('OPEN') and abs(m)<1e-9: continue
    if not k.startswith('OPEN') and all(abs(x)<1e-9 for x in v): continue
    print(f"{k:16s} n={n:3d} mean={m*100:+.3f}% se={sd/math.sqrt(n)*100:.3f}% w/l={sum(x>1e-9 for x in v)}/{sum(x<-1e-9 for x in v)}")
print(f"open A={oa:.4f} B={ob:.4f} d={ob-oa:+.4f}  classmean={tot/nn*100:+.4f}%  WA={wa}")
