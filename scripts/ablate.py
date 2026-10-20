import subprocess,os,glob,sys,concurrent.futures as cf
# ablate.py <labexe> <dir-glob> : zero each LW weight in turn and report mean kpi over files
exe=os.path.abspath(sys.argv[1]); files=sorted(sum([glob.glob(g) for g in sys.argv[2:]],[]))
base="1 -38.8 -22.4 13.5 47 -203 107 -7.25 2.78 2.26 -8.96 23.4 1.05e+03 69 0 46.7 27.9 50 0 0 0 0 0 0 0 20 50".split()
def run(w):
    env=dict(os.environ); env['SFW']=','.join(w)
    chunks=[files[j::8] for j in range(8)]
    def one(ch):
        out=subprocess.run([exe]+ch,capture_output=True,text=True,env=env).stdout
        return [float(l.split()[2][4:]) for l in out.splitlines() if 'kpi=' in l and ' w=' in l]
    with cf.ThreadPoolExecutor(8) as ex: ks=sum(ex.map(one,chunks),[])
    return sum(ks)/max(len(ks),1)
b=run(base); print('base',round(b,4),len(files),'files')
for i in range(27):
    if float(base[i])==0: continue
    w=base[:]; w[i]='0'
    print(i,base[i],'zeroed ->',round(run(w),4))
