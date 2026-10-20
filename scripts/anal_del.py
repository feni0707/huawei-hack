import sys,collections
def load(t):
    f=open(f'materials/tests/{t}.in').read().split();n=int(f[1]);return list(map(int,f[2:2+n]))
def run(t,win=None):
    c=load(t);n=len(c)
    ln=[0];birth=[0];alive_by_len=collections.defaultdict(list)  # list of pids in birth order
    import bisect
    pid=0;live=0
    stats=[]
    W=win or max(200,n//20)
    cur=collections.Counter()
    prev_place=None
    for i,x in enumerate(c):
        if x>0:
            pid+=1;ln.append(x);birth.append(i);alive_by_len[x].append(pid);cur['P']+=1;cur['Lsum']+=x;prev_place=pid
            cur[f'len{x}']+=0
        else:
            p=-x;L=ln[p];lst=alive_by_len[L]
            j=bisect.bisect_left(lst,p);rank=len(lst)-1-j
            lst.pop(j)
            cur['D']+=1
            cur['rank0']+= rank==0
            cur['rank<=2']+= rank<=2
            cur['immediate']+= (p==prev_place)
            cur['ageSum']+= i-birth[p]
            cur['rankSum']+=rank
            prev_place=None
        if (i+1)%W==0:
            stats.append((i+1,cur));cur=collections.Counter()
    print(f"test {t} n={n} pids={pid}")
    for i,cu in stats:
        D=cu['D'];P=cu['P']
        print(f"  upto {i:6d} P={P:4d} D={D:4d} avgPlaceLen={cu['Lsum']/max(P,1):5.1f} del_rank0={cu['rank0']/max(D,1):.2f} rank<=2={cu['rank<=2']/max(D,1):.2f} immediate={cu['immediate']/max(D,1):.2f} meanRank={cu['rankSum']/max(D,1):.1f} meanAge={cu['ageSum']/max(D,1):.0f}")
for t in sys.argv[1:]:
    run(int(t))
