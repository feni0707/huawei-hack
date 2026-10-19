import sys
def load(t):
    f=open(f'materials/tests/{t}.in').read().split()
    n=int(f[1]);return list(map(int,f[2:2+n]))
c=load(15)
pre=sum(1 for x in c[:4999] if x>0)
ln=[0]+[x for x in c if x>0]
dth={};bth={};pid=0
for i,x in enumerate(c):
    if x>0: pid+=1;bth[pid]=i
    else: dth[-x]=i
D=[0]+list(map(int,open(sys.argv[1]).read().split('\n')[3].split()))
cap=int(sys.argv[2]) if len(sys.argv)>2 else 20760
far=[(p,D[p],ln[p],bth[p],dth.get(p)) for p in range(1,pre+1) if D[p]+ln[p]>cap]
print(len(far),far[:25])
print('max prefix end',max(D[p]+ln[p] for p in range(1,pre+1)), 'last death among prefix', max(dth[p] for p in range(1,pre+1) if p in dth), 'n',len(c))
