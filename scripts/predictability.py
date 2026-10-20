import sys,collections,bisect
def load(p):
    f=open(p).read().split();n=int(f[1]);return list(map(int,f[2:2+n]))
def tokens(c):
    ln=[0];stacks=collections.defaultdict(list);pid=0;tok=[]
    for x in c:
        if x>0:
            pid+=1;ln.append(x);stacks[x].append(pid);tok.append(('P',x))
        else:
            p=-x;L=ln[p];st=stacks[L];j=bisect.bisect_left(st,p);r=len(st)-1-j;st.pop(j);tok.append(('D',L,min(r,15)))
    return tok
def run(path,maxo=6):
    c=load(path);tok=tokens(c);n=len(tok)
    ctx=[collections.defaultdict(collections.Counter) for _ in range(maxo+1)]
    hit=0;cnt=0;hit_o=[0]*(maxo+1)
    hist=[]
    # exact next-token prediction by longest matching context (PPM-C style, no blending): predict argmax at longest order with data
    for i,t in enumerate(tok):
        pred=None;used=-1
        for o in range(min(maxo,len(hist)),-1,-1):
            key=tuple(hist[len(hist)-o:]) if o else ()
            d=ctx[o].get(key)
            if d:
                pred=d.most_common(1)[0][0];used=o;break
        if i>=50:
            cnt+=1; hit+= (pred==t)
        for o in range(0,min(maxo,len(hist))+1):
            key=tuple(hist[len(hist)-o:]) if o else ()
            ctx[o][key][t]+=1
        hist.append(t)
    # also type-only (P vs D) and placement-length accuracy
    return hit/max(cnt,1),n
for p in sys.argv[1:]:
    a,n=run(p);print(p.split('/')[-1],n,'next-token acc %.3f'%a)
