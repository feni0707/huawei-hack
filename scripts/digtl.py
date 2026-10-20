import sys,re,subprocess,os
# usage: digtl.py test.in [exe]  -> timeline of trench excess over running-peak bound at each dig
t=sys.argv[1]; exe=os.path.abspath(sys.argv[2] if len(sys.argv)>2 else 'build/bin/simc_sf.exe')
f=open(t).read().split(); n=int(f[1]); c=list(map(int,f[2:2+n]))
env=dict(os.environ); env['DIGS']='1'
out=subprocess.run([exe,t],capture_output=True,text=True,env=env).stdout
digs=[(int(m.group(1)),int(m.group(2)),int(m.group(3))) for m in re.finditer(r'dig ci=(\d+) L=\d+ N=\d+ T=(\d+) occ=(\d+)',out)]
live=0;ln=[0];peak=0;pk=[0]*n
for i,x in enumerate(c):
    if x>0: ln.append(x);live+=x
    else: live-=ln[-x]
    peak=max(peak,live);pk[i]=peak
prev=-1
for ci,T,occ in digs:
    b=(pk[ci]+ (ln[-1] if False else 0)+19)//20*20
    # bound includes the pipe being placed: live after placing
    ex=T-((pk[ci]+19)//20*20)
    if ex!=prev: print(f"ci={ci} T={T} bound={(pk[ci]+19)//20*20} excess={ex}"); prev=ex
print('final',digs[-1][1] if digs else 0,'bound',(peak+19)//20*20)
