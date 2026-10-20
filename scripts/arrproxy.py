import sys,collections,random
def load(p):
    f=open(p).read().split();n=int(f[1]);return list(map(int,f[2:2+n]))
def prep(path):
    c=load(path);n=len(c);half=n//2
    ln=[0];birth=[0];death={}
    for i,x in enumerate(c):
        if x>0: ln.append(x);birth.append(i)
        else: death[-x]=i
    sh=[p for p in range(1,len(ln)) if birth[p]<half and ln[p]<=7 and (p not in death or death[p]>=half)]
    vic={p for p in sh if p in death}
    return ln,sh,vic
def dead(seq,ln,vic):
    d=0;run=0
    for p in seq+[None]:
        if p is not None and p in vic: run+=ln[p]
        else:
            if 0<run<8: d+=run
            run=0
    return d
def blocks_sorted(sh,ln,B,rev):
    out=[]
    for i in range(0,len(sh),B):
        b=sh[i:i+B]; b=sorted(b,key=lambda p:(ln[p],p),reverse=rev); out+=b
    return out
def zones(sh,ln):
    out=[]
    for L in range(1,8): out+=[p for p in sh if ln[p]==L]
    return out
def zigzag(sh,ln,B):
    out=[]
    for i in range(0,len(sh),B):
        b=sorted(sh[i:i+B],key=lambda p:(ln[p],p)); lo=0;hi=len(b)-1;z=[]
        while lo<=hi:
            z.append(b[hi]);hi-=1
            if lo<=hi: z.append(b[lo]);lo+=1
        out+=z
    return out
for path in sys.argv[1:]:
    ln,sh,vic=prep(path)
    tot=sum(ln[p] for p in vic)
    res={'pid':dead(sh,ln,vic),'zones':dead(zones(sh,ln),ln,vic)}
    for B in (7,14,28,56,112): res[f'bsort{B}']=dead(blocks_sorted(sh,ln,B,False),ln,vic)
    for B in (14,56): res[f'bsortD{B}']=dead(blocks_sorted(sh,ln,B,True),ln,vic); res[f'zigzag{B}']=dead(zigzag(sh,ln,B),ln,vic)
    print(path.split('/')[-1],'victim cells',tot,' '.join(f'{k}:{v}' for k,v in res.items()))
