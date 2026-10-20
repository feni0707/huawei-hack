#!/usr/bin/env python3
"""sep-CMA-ES (чистый Python, без numpy) над весами LW линейной оценки scorefit_lab (SFW).
 es_weights.py TAG GENS TRAIN_GLOBS... [-- VAL_GLOBS...]   env: SIG, LAM, SEED, SF (доп. ручки лабы), IDX (индексы LW через запятую)
 Лог: build/es/<TAG>.log ; лучший вектор печатается в конце и пишется в build/es/<TAG>.best"""
import sys,os,glob,subprocess,math,random,concurrent.futures as cf,time
args=sys.argv[1:]
TAG=args[0]; GENS=int(args[1]); rest=args[2:]
if '--' in rest: i=rest.index('--'); tr=rest[:i]; va=rest[i+1:]
else: tr=rest; va=[]
TRAIN=sorted(sum([glob.glob(g) for g in tr],[])); VAL=sorted(sum([glob.glob(g) for g in va],[]))
exe=os.path.abspath('build/bin/lab.exe')
BASE=[float(x) for x in "1 -38.8 -22.4 13.5 47 -203 107 -7.25 2.78 2.26 -8.96 23.4 1.05e+03 69 0 46.7 27.9 50 0 0 0 0 0 0 0 20 50".split()]
if os.environ.get('W0'): BASE=[float(x) for x in os.environ['W0'].split()]
SC=[0,5,5,2,5,50,50,5,2,5,2,5,200,50,0,5,3,3,50,50,0,0,5,5,200,5,30]
IDX=[int(x) for x in os.environ['IDX'].split(',')] if os.environ.get('IDX') else [i for i in range(27) if SC[i]>0 and i not in (14,18,19,22,23,24)]
n=len(IDX)
def tow(x):
    w=BASE[:]
    for k,i in enumerate(IDX): w[i]=BASE[i]+SC[i]*x[k]
    return ','.join(f'{v:.6g}' for v in w)
def evalfiles(sfw,files):
    env=dict(os.environ); env['SFW']=sfw
    if 'SF' in os.environ: env['SF']=os.environ['SF']
    out=subprocess.run([exe]+files,capture_output=True,text=True,env=env).stdout
    r={}
    for l in out.splitlines():
        t=l.split()
        if len(t)>=3 and t[1].startswith('w=') and t[2].startswith('kpi='): r[t[0]]=float(t[2][4:])
    return r
def score(sfw,files):
    r=evalfiles(sfw,files)
    return sum(r.get(f,0.0) for f in files)/len(files), r
lam=int(os.environ.get('LAM','16')); mu=lam//2
w=[math.log(mu+0.5)-math.log(i+1) for i in range(mu)]; sw=sum(w); w=[x/sw for x in w]; mueff=1/sum(x*x for x in w)
cs=(mueff+2)/(n+mueff+5); damps=1+2*max(0,math.sqrt((mueff-1)/(n+1))-1)+cs
cc=(4+mueff/n)/(n+4+2*mueff/n)
c1=2/((n+1.3)**2+mueff)*(n+1.5)/3; cmu=min(1-c1,2*(mueff-2+1/mueff)/((n+2)**2+mueff)*(n+1.5)/3)
chiN=math.sqrt(n)*(1-1/(4*n)+1/(21*n*n))
m=[0.0]*n; sigma=float(os.environ.get('SIG','0.4')); C=[1.0]*n; pc=[0.0]*n; ps=[0.0]*n
rnd=random.Random(int(os.environ.get('SEED','1')))
os.makedirs('build/es',exist_ok=True)
log=open(f'build/es/{TAG}.log','a')
def P(*a):
    print(*a,file=log,flush=True)
b0,r0=score(tow(m),TRAIN); P('base train',round(b0,5),len(TRAIN),'n',n)
if VAL: v0,_=score(tow(m),VAL); P('base val',round(v0,5),len(VAL))
best=(b0,m[:])
ex=cf.ThreadPoolExecutor(int(os.environ.get('J','12')))
t0=time.time()
for g in range(GENS):
    D=[math.sqrt(c) for c in C]
    Z=[[rnd.gauss(0,1) for _ in range(n)] for _ in range(lam)]
    Y=[[Z[k][j]*D[j] for j in range(n)] for k in range(lam)]
    X=[[m[j]+sigma*Y[k][j] for j in range(n)] for k in range(lam)]
    fs=list(ex.map(lambda x: score(tow(x),TRAIN)[0], X))
    order=sorted(range(lam),key=lambda k:-fs[k])
    if fs[order[0]]>best[0]: best=(fs[order[0]],X[order[0]][:])
    yw=[sum(w[i]*Y[order[i]][j] for i in range(mu)) for j in range(n)]
    m=[m[j]+sigma*yw[j] for j in range(n)]
    ps=[(1-cs)*ps[j]+math.sqrt(cs*(2-cs)*mueff)*yw[j]/D[j] for j in range(n)]
    nps=math.sqrt(sum(x*x for x in ps))
    hs=nps/math.sqrt(1-(1-cs)**(2*(g+1)))/chiN<1.4+2/(n+1)
    pc=[(1-cc)*pc[j]+(math.sqrt(cc*(2-cc)*mueff)*yw[j] if hs else 0) for j in range(n)]
    C=[(1-c1-cmu)*C[j]+c1*(pc[j]**2+(0 if hs else cc*(2-cc)*C[j]))+cmu*sum(w[i]*Y[order[i]][j]**2 for i in range(mu)) for j in range(n)]
    sigma*=math.exp((cs/damps)*(nps/chiN-1))
    P(f'gen {g} best_gen {fs[order[0]]:.5f} mean {sum(fs)/lam:.5f} best {best[0]:.5f} sigma {sigma:.3f} t {time.time()-t0:.0f}s')
    if VAL and (g%10==9 or g==GENS-1):
        vb,_=score(tow(best[1]),VAL); vm,_=score(tow(m),VAL); P(f'  VAL best {vb:.5f} mean-vec {vm:.5f}')
P('BEST',best[0],tow(best[1]))
open(f'build/es/{TAG}.best','w').write(tow(best[1])+'\n')
print(tow(best[1]))
