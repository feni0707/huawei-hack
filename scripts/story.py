import sys
t=int(sys.argv[1]); a=int(sys.argv[2]) if len(sys.argv)>2 else 0; b=int(sys.argv[3]) if len(sys.argv)>3 else 10**9
f=open(f'materials/tests/{t}.in').read().split();n=int(f[1]);c=list(map(int,f[2:2+n]))
ln=[0];pid=0;out=[];alive=set()
i=0
tok=[]
for i,x in enumerate(c):
    if x>0:
        pid+=1;ln.append(x);alive.add(pid);tok.append(f"+{x}")
    else:
        p=-x;alive.discard(p)
        # relative: how many pids ago?
        tok.append(f"-{ln[p]}[{pid-p}]")
# compress runs of identical tokens
res=[];last=None;cnt=0
for k,x in enumerate(tok):
    if k<a or k>=b: continue
    if x==last: cnt+=1
    else:
        if last: res.append(last+(f"x{cnt}" if cnt>1 else ""))
        last=x;cnt=1
res.append(last+(f"x{cnt}" if cnt>1 else ""))
print(' '.join(res))
