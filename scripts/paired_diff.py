#!/usr/bin/env python3
# paired.py "SFA|SFWA" "SFB|SFWB" glob... : парная разница B-A по классам: среднее, ст. ошибка, доля побед
import subprocess, sys, glob, re, os, collections, math
from concurrent.futures import ThreadPoolExecutor
binp=os.environ.get('BIN','build/sim_sf3')
files=[f for g in sys.argv[3:] for f in sorted(glob.glob(g))]
def run(spec, fl):
    sf,_,sw=spec.partition('|'); env=dict(os.environ,SF=sf)
    if sw: env['SFW']=sw
    out=subprocess.run([binp]+fl,env=env,capture_output=True,text=True).stdout
    return {m.group(1):float(m.group(2)) for m in re.finditer(r'(\S+)\s+w=\S+\s+kpi=(\S+)',out)}
def runall(spec):
    parts=[files[i::8] for i in range(8)]
    with ThreadPoolExecutor(8) as ex: rs=list(ex.map(lambda p: run(spec,p), parts))
    d={}
    for r in rs: d.update(r)
    return d
a=runall(sys.argv[1]); b=runall(sys.argv[2])
c=collections.defaultdict(list)
for f in files: c[f.split('/')[-2]].append(b[f]-a[f])
for k,v in sorted(c.items()):
    n=len(v); m=sum(v)/n; sd=math.sqrt(sum((x-m)**2 for x in v)/max(1,n-1)); 
    print(f"{k:14s} n={n:3d} mean={m*100:+.3f}% se={sd/math.sqrt(n)*100:.3f}% wins={sum(x>0 for x in v)} losses={sum(x<0 for x in v)}")
