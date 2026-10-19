#!/usr/bin/env python3
# Пример: W0="..." SIG=0.4 python3 scripts/cmaes_weights.py run1 150   (нужен build/sim_sf из scorefit_lab.cpp -DSIM,
# наборы build/lifo/{a,b} — gen_tests.py --class lifo_mix --count 90 --seed 200/201, build/big/s100..s111)
# CMA-ES над весами линейной оценки (P14=6) в build/sim_sf. Лог: build/cma/<tag>.log
import numpy as np, subprocess, glob, re, os, sys, collections, json
from concurrent.futures import ThreadPoolExecutor
TAG=sys.argv[1] if len(sys.argv)>1 else 'run'
GENS=int(sys.argv[2]) if len(sys.argv)>2 else 120
BASE_SF=os.environ.get('BASE_SF','14=6 15=0')
W0=np.array([1,0,0,0,0,0,0,0,0,0,0,0,1000,150,0,0,0,3,0,0,0,0,0,0,0,0,0],float)
SC=np.array([0,5,5,2,5,50,50,5,2,5,2,5,200,50,0,5,3,3,50,50,0,0,5,5,200,5,30],float)
if os.environ.get('SC'): SC=np.array(list(map(float,os.environ['SC'].split())),float)
if os.environ.get('W0'): W0=np.array(list(map(float,os.environ['W0'].split())),float)
if os.environ.get('SIG'): sigma0=float(os.environ['SIG'])
idx=[i for i in range(len(W0)) if SC[i]>0]
def files(sets):
    out=[]
    for g in sets: out+=sorted(glob.glob(g))
    return out
TRAIN=files(os.environ.get('TRAIN_EXTRA','').split()+['build/lifo/a/*/*.in']+[f'build/big/s10{i}/{c}/*.in' for i in range(6) for c in ['alternating','two_phase','phases','sorted_groups','periodic']])
VAL=files(os.environ.get('VAL_EXTRA','').split()+['build/lifo/b/*/*.in']+[f'build/big/{s}/{c}/*.in' for s in ['s106','s107','s108','s109','s110','s111'] for c in ['alternating','two_phase','phases','sorted_groups','periodic']])
CW={'two_phase_real':10,'two_phase':10,'alternating':4,'lifo_mix':4,'phases':4,'sorted_groups':1,'periodic':1}
def chunks(l,k): return [l[i::k] for i in range(k)]
def evalw(w, fl, nsplit=1):
    sw=' '.join(f'{x:.6g}' for x in w)
    env=dict(os.environ,SF=BASE_SF,SFW=sw,SFQ='')
    del env['SFQ']
    res={}
    for part in chunks(fl,nsplit):
        out=subprocess.run(['build/sim_sf']+part,env=env,capture_output=True,text=True).stdout
        for m in re.finditer(r'(\S+)\s+w=\S+\s+kpi=(\S+)',out): res[m.group(1)]=float(m.group(2))
    c=collections.defaultdict(list)
    for f,v in res.items(): c[f.split('/')[-2]].append(v)
    av={k:sum(v)/len(v) for k,v in c.items()}
    return sum(CW[k]*av[k] for k in av), av
def tow(x):
    w=W0.copy(); w[idx]+=SC[idx]*x; return w
n=len(idx); lam=12; mu=lam//2
wts=np.log(mu+0.5)-np.log(np.arange(1,mu+1)); wts/=wts.sum(); mueff=1/np.sum(wts**2)
cc=(4+mueff/n)/(n+4+2*mueff/n); cs=(mueff+2)/(n+mueff+5); c1=2/((n+1.3)**2+mueff)
cmu=min(1-c1,2*(mueff-2+1/mueff)/((n+2)**2+mueff)); damps=1+2*max(0,np.sqrt((mueff-1)/(n+1))-1)+cs
chiN=np.sqrt(n)*(1-1/(4*n)+1/(21*n*n))
m=np.zeros(n); sigma=float(os.environ.get('SIG','0.6')); pc=np.zeros(n); ps=np.zeros(n); C=np.eye(n)
rng=np.random.default_rng(int(os.environ.get('SEED','1')))
log=open(f'build/cma/{TAG}.log','a')
ex=ThreadPoolExecutor(int(os.environ.get('J','9')))
b0,av0=evalw(tow(m),TRAIN,4); print('base train',b0,av0,file=log,flush=True)
best=(b0,m.copy())
for g in range(GENS):
    D2,B=np.linalg.eigh(C); D=np.sqrt(np.maximum(D2,1e-20))
    z=rng.standard_normal((lam,n)); y=z@np.diag(D)@B.T; X=m+sigma*y
    fs=list(ex.map(lambda x: evalw(tow(x),TRAIN)[0], X))
    order=np.argsort(fs)[::-1]   # максимизация
    if fs[order[0]]>best[0]: best=(fs[order[0]],X[order[0]].copy())
    yw=np.sum(wts[:,None]*y[order[:mu]],axis=0)
    m=m+sigma*yw
    Cih=B@np.diag(1/D)@B.T
    ps=(1-cs)*ps+np.sqrt(cs*(2-cs)*mueff)*(Cih@yw)
    hs=np.linalg.norm(ps)/np.sqrt(1-(1-cs)**(2*(g+1)))/chiN<1.4+2/(n+1)
    pc=(1-cc)*pc+hs*np.sqrt(cc*(2-cc)*mueff)*yw
    C=(1-c1-cmu)*C+c1*(np.outer(pc,pc)+(1-hs)*cc*(2-cc)*C)+cmu*sum(wts[i]*np.outer(y[order[i]],y[order[i]]) for i in range(mu))
    sigma*=np.exp((cs/damps)*(np.linalg.norm(ps)/chiN-1))
    print(f'gen {g} best_gen {fs[order[0]]:.4f} mean_gen {np.mean(fs):.4f} best {best[0]:.4f} sigma {sigma:.3f} w {" ".join(f"{v:.3g}" for v in tow(best[1]))}',file=log,flush=True)
    if g%10==9 or g==GENS-1:
        vb,avb=evalw(tow(best[1]),VAL,4); vm,avm=evalw(tow(m),VAL,4)
        print(f'  VAL best {vb:.4f} {json.dumps({k:round(v,4) for k,v in avb.items()})}  mean {vm:.4f}',file=log,flush=True)
v0,avv=evalw(tow(np.zeros(n)),VAL,4); print('base VAL',v0,avv,file=log,flush=True)
