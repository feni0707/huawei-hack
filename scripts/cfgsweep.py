#!/usr/bin/env python3
"""cfgsweep.py exe -- dir... ; configs from CFGS env (name=SF-string;name2=...) -> per-class means per config + oracle best-of"""
import sys,os,subprocess,glob,collections,math,json,concurrent.futures as cf
args=sys.argv[1:]; i=args.index('--'); exe=os.path.abspath(args[0]); dirs=args[i+1:]
cfgs=[c.split('=',1) for c in os.environ['CFGS'].split(';')] if 'CFGS' in os.environ else [('default','')]
files=[]
for d in dirs:
    for cd in sorted(glob.glob(os.path.join(d,'*'))):
        if os.path.isdir(cd):
            for f in sorted(glob.glob(os.path.join(cd,'*.in'))): files.append((os.path.basename(cd),f))
def run(cfg):
    name,sf=cfg; env=dict(os.environ)
    # cfg string "SF:...|SFW:..." or plain SF
    for part in sf.split('|'):
        if part.startswith('SFW:'): env['SFW']=part[4:]
        elif part.startswith('ENV:'):
            k,v=part[4:].split('=',1); env[k]=v
        elif part: env['SF']=part
    res={}
    chunks=[files[j::12] for j in range(12)]
    def one(ch):
        out=subprocess.run([exe]+[f for _,f in ch],capture_output=True,text=True,env=env).stdout
        r={}
        for line in out.splitlines():
            t=line.split()
            if len(t)>=3 and t[2].startswith('kpi='): r[t[0]]=float(t[2][4:])
        return r
    with cf.ThreadPoolExecutor(12) as ex:
        for r in ex.map(one,chunks): res.update(r)
    return res
R={name:run((name,sf)) for name,sf in cfgs}
json.dump(R,open('build/cfgsweep_last.json','w'))
cls=collections.defaultdict(list)
for c,f in files: cls[c].append(f)
names=[n for n,_ in cfgs]
print('class'.ljust(16),'n'.rjust(3),''.join(f'{n[:9]:>10}' for n in names),'   oracle  gain')
W={'two_phase':10,'two_phase_real':10,'alternating':4,'lifo_mix':4,'phases':4,'sorted_groups':2,'periodic':2,'random_bursts':2}
tot=collections.Counter()
for c,fs in sorted(cls.items()):
    m=[sum(R[n][f] for f in fs)/len(fs) for n in names]
    orc=sum(max(R[n][f] for n in names) for f in fs)/len(fs)
    print(c.ljust(16),str(len(fs)).rjust(3),''.join(f'{x:10.4f}' for x in m),f' {orc:.4f} {100*(orc-m[0]):+.3f}%')
    for n,x in zip(names,m): tot[n]+=W.get(c,1)*x
    tot['oracle']+=W.get(c,1)*orc
print('obj'.ljust(20),''.join(f'{tot[n]:10.3f}' for n in names),f' {tot["oracle"]:.3f}')
