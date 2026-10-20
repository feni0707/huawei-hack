import subprocess,os,glob,sys,collections,concurrent.futures as cf
# wsweep.py <labexe> "<name>=<i>:<v>,<i>:<v>;..." -- dirs...  : override LW entries; prints per-class mean (and obj)
exe=os.path.abspath(sys.argv[1]); spec=sys.argv[2]; dirs=sys.argv[sys.argv.index('--')+1:]
base="1 -38.8 -22.4 13.5 47 -203 107 -7.25 2.78 2.26 -8.96 23.4 1.05e+03 69 0 46.7 27.9 50 0 0 0 0 0 0 0 20 50".split()
files=[]
for d in dirs:
    for cd in sorted(glob.glob(os.path.join(d,'*'))):
        if os.path.isdir(cd):
            for f in sorted(glob.glob(os.path.join(cd,'*.in'))): files.append((os.path.basename(cd),f))
    for f in sorted(glob.glob(os.path.join(d,'*.in'))): files.append(('open',f))
W={'two_phase':10,'two_phase_real':10,'alternating':4,'lifo_mix':4,'phases':4,'sorted_groups':2,'periodic':2,'random_bursts':2}
def run(name,ovr):
    w=base[:]
    for kv in ovr.split(','):
        if kv: i,v=kv.split(':'); w[int(i)]=v
    env=dict(os.environ); env['SFW']=','.join(w)
    chunks=[files[j::12] for j in range(12)]
    def one(ch):
        out=subprocess.run([exe]+[f for _,f in ch],capture_output=True,text=True,env=env).stdout
        return {l.split()[0]:float(l.split()[2][4:]) for l in out.splitlines() if ' w=' in l and 'kpi=' in l}
    r={}
    with cf.ThreadPoolExecutor(12) as ex:
        for x in ex.map(one,chunks): r.update(x)
    return r
res={}
for part in spec.split(';'):
    name,ovr=part.split('=',1); res[name]=run(name,ovr)
names=list(res); cls=collections.defaultdict(list)
for c,f in files: cls[c].append(f)
print('class'.ljust(18),'n'.rjust(3),''.join(f'{n[:10]:>11}' for n in names))
tot=collections.Counter()
for c,fs in sorted(cls.items()):
    m=[sum(res[n][f] for f in fs)/len(fs) for n in names]
    print(c.ljust(18),str(len(fs)).rjust(3),''.join(f'{x:11.4f}' for x in m))
    for n,x in zip(names,m): tot[n]+=W.get(c,1)*x
print('obj'.ljust(22),''.join(f'{tot[n]:11.3f}' for n in names))
