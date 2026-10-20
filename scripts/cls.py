#!/usr/bin/env python3
"""Per-class comparison of simc-style executables: cls.py exe1[,ENV=..] exe2 ... -- dir1 dir2 ...
Each dir contains class subdirs with *.in. Prints mean KPI per class (and gap to UB) per exe + paired diff vs first."""
import sys,os,subprocess,collections,math,glob,concurrent.futures as cf
args=sys.argv[1:]
i=args.index('--'); exes=args[:i]; dirs=args[i+1:]
def ub(path):
    f=open(path).read().split(); n=int(f[1]); c=list(map(int,f[2:2+n])); ln=[0];live=0;pk=0
    for x in c:
        if x>0: ln.append(x);live+=x;pk=max(pk,live)
        else: live-=ln[-x]
    return live/(((pk+19)//20)*20)
files=[]
for d in dirs:
    for cd in sorted(glob.glob(os.path.join(d,'*'))):
        if os.path.isdir(cd):
            for f in sorted(glob.glob(os.path.join(cd,'*.in'))): files.append((os.path.basename(cd),f))
def runexe(spec):
    parts=spec.split(','); exe=os.path.abspath(parts[0]); env=dict(os.environ)
    for p in parts[1:]:
        k,v=p.split('=',1); env[k]=v
    res={}
    chunks=[files[j::8] for j in range(8)]
    def one(ch):
        out=subprocess.run([exe]+[f for _,f in ch],capture_output=True,text=True,env=env).stdout
        r={}
        for line in out.splitlines():
            t=line.split()
            if len(t)>=3 and t[2].startswith('kpi='): r[t[0]]=(float(t[2][4:]),'WA' in line)
        return r
    with cf.ThreadPoolExecutor(8) as ex:
        for r in ex.map(one,chunks): res.update(r)
    return res
results=[runexe(e) for e in exes]
ubs={f:ub(f) for _,f in files}
cls=collections.defaultdict(list)
for c,f in files: cls[c].append(f)
W={'two_phase':10,'two_phase_real':10,'alternating':4,'lifo_mix':4,'phases':4,'sorted_groups':2,'periodic':2,'random_bursts':2}
hdr='class'.ljust(16)+'n'.rjust(4)+'   UB   '+''.join(f' {os.path.basename(e.split(",")[0])[:11]:>11}' for e in exes)+'   paired diff vs first (±se)'
print(hdr)
tot=[0.0]*len(exes)
for c,fs in sorted(cls.items()):
    m=[sum(r[f][0] for f in fs)/len(fs) for r in results]
    line=c.ljust(16)+str(len(fs)).rjust(4)+f' {sum(ubs[f] for f in fs)/len(fs):.4f}'+''.join(f' {x:11.4f}' for x in m)
    for k in range(1,len(exes)):
        d=[results[k][f][0]-results[0][f][0] for f in fs]; mu=sum(d)/len(d); sd=math.sqrt(sum((x-mu)**2 for x in d)/max(len(d)-1,1))/math.sqrt(len(d))
        line+=f'   {mu*100:+.3f}%±{sd*100:.3f}'
    wa=sum(1 for r in results for f in fs if r[f][1])
    if wa: line+=f'  WA={wa}'
    print(line)
    for k in range(len(exes)): tot[k]+=W.get(c,1)*m[k]
print('obj'.ljust(16),'    ','      ',''.join(f' {x:11.4f}' for x in tot))
