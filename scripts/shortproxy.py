import sys,collections,random
t=int(sys.argv[1]) if len(sys.argv)>1 else 15
f=open(f'materials/tests/{t}.in').read().split();n=int(f[1]);c=list(map(int,f[2:2+n]))
half=n//2
ln=[0];birth=[0];death={}
for i,x in enumerate(c):
    if x>0: ln.append(x);birth.append(i)
    else: death[-x]=i
# shorts: born in phase 1 (<half) and length<=7; alive at phase switch; victims = die at >= half
sh=[p for p in range(1,len(ln)) if birth[p]<half and ln[p]<=7 and (p not in death or death[p]>=half)]
vic={p for p in sh if p in death}
print('shorts alive at switch',len(sh),'victims',len(vic),'frac',len(vic)/len(sh), 'victim cells',sum(ln[p] for p in vic),'total short cells',sum(ln[p] for p in sh))
cnt=collections.Counter(ln[p] for p in sh); print('len dist',sorted(cnt.items()))
vl=collections.Counter(ln[p] for p in vic); print('victim frac by len',{k:round(vl[k]/cnt[k],2) for k in sorted(cnt)})
def dead(seq):
    d=0;run=0;tot=0
    for p in seq+[None]:
        if p is not None and p in vic: run+=ln[p]
        else:
            if 0<run<8: d+=run
            if run>=8: tot+=run
            run=0
    return d,tot
# alive-at-switch pids in pid order (what best fit ~ gives)
print('PID order dead cells',dead(sh))
# zones by length
z=[]
for L in range(1,8): z+= [p for p in sh if ln[p]==L]
print('zones by length',dead(z))
# zones by length with blocks of 16 slots alternating? (skip)
# random order
r=sh[:];random.seed(1);random.shuffle(r);print('random order',dead(r))
# sorted by length desc interleaved pairing big-small
big=[p for p in sh if ln[p]>=5];small=[p for p in sh if ln[p]<5]
def interleave(a,b):
    out=[];i=j=0
    while i<len(a) or j<len(b):
        if i<len(a): out.append(a[i]);i+=1
        if j<len(b): out.append(b[j]);j+=1
    return out
print('interleave big/small (pid order each)',dead(interleave(big,small)))
# complement pairing: for each big of length L, partner small of length >=8-L
import bisect
pools={L:[p for p in sh if ln[p]==L] for L in range(1,8)}
# oracle: knowing victims, best grouping dead lower bound
vc=sum(ln[p] for p in vic); print('victim cells',vc)
print('--- correlation')
shs=set(sh)
order=sh
def cond(seq,name):
    a=b=c2=d=0
    for x,y in zip(seq,seq[1:]):
        vx=x in vic; vy=y in vic
        if vx: a+=1; b+= vy
        else: c2+=1; d+= vy
    print(name,'P(v|prev v)=%.3f P(v|prev nv)=%.3f'%(b/max(a,1),d/max(c2,1)))
cond(sh,'pid-neighbors')
for L in (1,4,7):
    cond([p for p in sh if ln[p]==L],f'same-len L={L}')
# victims density over pid index (deciles)
import math
dec=[0]*20
for k,p in enumerate(sh):
    dec[k*20//len(sh)] += p in vic
print('victims per 5% of shorts stack:',dec)
print('--- patterns')
for L in (4,):
    seq=[p for p in sh if ln[p]==L]
    s=''.join('V' if p in vic else '.' for p in seq)
    # print the top-25% part and the area around 30-40%
    print(L,len(seq)); print(s[-1700:-1200]); print(s[int(len(s)*0.3):int(len(s)*0.3)+300])
print('--- real-victim arrangements')
st={L:[p for p in sh if ln[p]==L] for L in range(1,8)}
def pairz(stA,stB):
    out=[]
    for k in range(max(len(stA),len(stB))):
        if k<len(stA): out.append(stA[k])
        if k<len(stB): out.append(stB[k])
    return out
pz=[]
for a,b in ((1,7),(2,6),(3,5)): pz+=pairz(st[a],st[b])
pz+=st[4]
print('pair zones (stack-aligned)',dead(pz))
z=[]
for L in range(1,8): z+=st[L]
print('zones',dead(z))
# pair adjacent-in-time (pid-near) partner: pair each X with next Y arriving
import itertools
def timepair():
    out=[];wait={L:[] for L in range(1,8)}
    for p in sh:
        L=ln[p];Y=8-L
        if Y!=L and wait[Y]:
            q=wait[Y].pop(0); out.append(q);out.append(p)
        else:
            wait[L].append(p)
            if Y==L and len(wait[L])==2:
                a,b=wait[L];wait[L]=[];out.append(a);out.append(b)
    for L in wait: out+=wait[L]
    return out
print('time-paired bins',dead(timepair()))
