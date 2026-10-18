import sys, collections
def run(path, k=10):
    t=open(path).read().split(); n=int(t[1]); cmds=list(map(int,t[2:2+n]))
    alive={}; pid=0
    print("==", path, "n", n)
    step=max(1,n//k)
    for lo in range(0,n,step):
        pl=[]; dl=[]; dage=[]
        for i in range(lo,min(n,lo+step)):
            c=cmds[i]
            pass
        # recompute with state: need full pass; do single pass storing per-window lists
    # single pass
    alive={}; pid=0; wins=[([],[],[]) for _ in range(k+1)]
    for i,c in enumerate(cmds):
        w=min(k,i//step)
        if c>0:
            pid+=1; alive[pid]=c; wins[w][0].append(c)
        else:
            L=alive.pop(-c); wins[w][1].append(L); wins[w][2].append(pid-(-c))
    for w,(pl,dl,da) in enumerate(wins):
        if not pl and not dl: continue
        def st(a): 
            return "min%d max%d avg%.1f"%(min(a),max(a),sum(a)/len(a)) if a else "-"
        med = sorted(da)[len(da)//2] if da else 0
        print(f" win{w}: places {len(pl):5d} [{st(pl)}] dels {len(dl):5d} [{st(dl)}] median age(pid) {med}")
for p in sys.argv[1:]: run(p)
