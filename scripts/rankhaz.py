import sys,collections,bisect
def load(p):
    f=open(p).read().split();n=int(f[1]);return list(map(int,f[2:2+n]))
def run(path,Ws=(10,50,300),sample=7):
    c=load(path);n=len(c)
    ln=[0];birth=[0];death={}
    for i,x in enumerate(c):
        if x>0: ln.append(x);birth.append(i)
        else: death[-x]=i
    stacks=collections.defaultdict(list)  # length -> alive pids ascending
    pid=0
    tot=collections.Counter();die={W:collections.Counter() for W in Ws}
    for i,x in enumerate(c):
        if x>0:
            pid+=1;stacks[x].append(pid)
        else:
            p=-x;stacks[ln[p]].remove(p)
        if i%sample==0:
            for L,st in stacks.items():
                m=len(st)
                for r in range(min(m,12)):
                    p=st[m-1-r]
                    k=r
                    tot[k]+=1
                    d=death.get(p,10**9)
                    for W in Ws:
                        if d-i<=W: die[W][k]+=1
    print(path)
    for W in Ws:
        print(' W=%d'%W,' '.join('r%d:%.3f'%(k,die[W][k]/max(tot[k],1)) for k in range(12)))
for p in sys.argv[1:]: run(p)
