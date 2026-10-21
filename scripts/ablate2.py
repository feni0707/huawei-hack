import subprocess,os,glob,sys,concurrent.futures as cf
exe=os.path.abspath(sys.argv[1]); groups=sys.argv[2:]
base="1 -38.8 -22.4 13.5 47 -203 107 -7.25 2.78 2.26 -8.96 23.4 1.05e+03 69 0 46.7 27.9 50 0 0 0 0 0 0 0 20 50".split()
def run(w,files):
    env=dict(os.environ); env['SFW']=','.join(w)
    out=subprocess.run([exe]+files,capture_output=True,text=True,env=env).stdout
    k=[float(l.split()[2][4:]) for l in out.splitlines() if ' w=' in l]
    return sum(k)/len(k)
G={g:sorted(glob.glob(g)) for g in groups}
cfgs=[('base',base)]+[(f'z{i}',base[:i]+['0']+base[i+1:]) for i in range(1,27) if float(base[i])!=0]+[('bf',['1']+['0']*26)]
with cf.ThreadPoolExecutor(8) as ex:
    res={(n,g):ex.submit(run,w,fs) for n,w in cfgs for g,fs in G.items()}
print('cfg'.ljust(6),''.join(f'{os.path.basename(os.path.dirname(g.rstrip("/*.in")))[:10]:>11}' for g in G))
for n,w in cfgs:
    print(n.ljust(6),''.join(f'{res[(n,g)].result():11.4f}' for g in G))
