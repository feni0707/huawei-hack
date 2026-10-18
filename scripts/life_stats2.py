import sys, collections
t = open(sys.argv[1]).read().split()
w, n = float(t[0]), int(t[1]); cmds = list(map(int, t[2:2+n]))
alive = {}; pid = 0
ageS=[]; ageL=[]; rankS=[]; rankL=[]
import bisect
alive_sorted=[]  # list of alive pids sorted (they are appended in increasing order)
for i,c in enumerate(cmds):
    if c>0:
        pid+=1; alive[pid]=c; alive_sorted.append(pid)
    else:
        p=-c; L=alive.pop(p)
        j=bisect.bisect_left(alive_sorted,p)
        frac=(len(alive_sorted)-1-j)/max(1,len(alive_sorted))  # 0=newest, ->1 oldest
        alive_sorted.pop(j)
        (ageS if L<=7 else ageL).append(pid-p)
        (rankS if L<=7 else rankL).append(frac)
def q(a, ps=(0.05,.25,.5,.75,.9,.99)):
    a=sorted(a); return [a[min(len(a)-1,int(p*len(a)))] for p in ps]
print("short age pids q:", q(ageS)); print("long age pids q:", q(ageL))
print("short rank(0=newest) q:", [round(x,3) for x in q(rankS)]); print("long rank q:", [round(x,3) for x in q(rankL)])
# histogram of rank fraction in 10 bins
for name,r in (("short",rankS),("long",rankL)):
    h=[0]*10
    for x in r: h[min(9,int(x*10))]+=1
    print(name,"rank hist (newest..oldest):",h)
