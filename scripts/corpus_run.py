#!/usr/bin/env python3
"""corpus_run.py run <name> <exe> [ENV=V ...]   — прогнать simc-exe по всему корпусу build/corpus (+ open), сохранить build/corpus_res/<name>.json
   corpus_run.py cmp <nameA> <nameB> [--top K] [--split even|odd|all] — сравнение по семействам и худшие тесты для B относительно A.
Семейство теста = <набор без seed>/<класс>, например gen/steady, g/alternating, res/residual, sc/scen, rb/pairs."""
import sys, os, glob, json, re, subprocess, collections, math, concurrent.futures as cf
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
def corpus():
    fs = sorted(glob.glob('build/corpus/*/*/*.in')) + sorted(glob.glob('materials/tests/*.in'))
    return [f.replace('\\', '/') for f in fs]
def family(f):
    if f.startswith('materials'): return 'open/open'
    p = f.split('/'); return re.sub(r'\d+$', '', p[2]) + '/' + p[3]
def ub(f):
    t = open(f).read().split(); n = int(t[1]); c = list(map(int, t[2:2 + n])); ln = [0]; live = 0; pk = 0
    for x in c:
        if x > 0: ln.append(x); live += x; pk = max(pk, live)
        else: live -= ln[-x]
    return live / (((pk + 19) // 20) * 20) if pk else 1.0
def run(name, exe, envs):
    files = corpus(); env = dict(os.environ)
    for e in envs: k, v = e.split('=', 1); env[k] = v
    exe = os.path.abspath(exe)
    chunks = [files[i::24] for i in range(24)]
    def one(ch):
        out = subprocess.run([exe] + ch, capture_output=True, text=True, env=env).stdout
        r = {}
        for l in out.splitlines():
            t = l.split()
            if len(t) >= 4 and t[2].startswith('kpi='): r[t[0].replace('\\', '/')] = (float(t[2][4:]), int(t[3][2:]), 'WA' in l)
        return r
    res = {}
    with cf.ThreadPoolExecutor(int(os.environ.get('J', '14'))) as ex:
        for r in ex.map(one, chunks): res.update(r)
    os.makedirs('build/corpus_res', exist_ok=True)
    json.dump(res, open(f'build/corpus_res/{name}.json', 'w'))
    miss = [f for f in files if f not in res]
    wa = [f for f, v in res.items() if v[2]]
    print(f'{name}: {len(res)} tests, missing {len(miss)}, WA {len(wa)}', miss[:3], wa[:3])
def cmp(a, b, top, split):
    A = json.load(open(f'build/corpus_res/{a}.json')); B = json.load(open(f'build/corpus_res/{b}.json'))
    files = [f for f in corpus() if f in A and f in B]
    if split != 'all':
        files = [f for i, f in enumerate(files) if (i % 2 == 0) == (split == 'even')]
    fam = collections.defaultdict(list)
    for f in files: fam[family(f)].append(f)
    print(f'{"family":28s} {"n":>4} {a[:10]:>10} {b[:10]:>10}   diff% ±se   B>A/B<A')
    tot = []
    for k in sorted(fam):
        fs = fam[k]; d = [B[f][0] - A[f][0] for f in fs]; mu = sum(d) / len(d)
        se = math.sqrt(sum((x - mu) ** 2 for x in d) / max(1, len(d) - 1) / len(d))
        tot += d
        print(f'{k:28s} {len(fs):4d} {sum(A[f][0] for f in fs)/len(fs):10.4f} {sum(B[f][0] for f in fs)/len(fs):10.4f}  {100*mu:+.3f}±{100*se:.3f}  {sum(x>1e-9 for x in d)}/{sum(x<-1e-9 for x in d)}')
    mu = sum(tot) / len(tot)
    print(f'{"ALL":28s} {len(tot):4d}  mean diff {100*mu:+.4f}%   B>A {sum(x>1e-9 for x in tot)}  B<A {sum(x<-1e-9 for x in tot)}')
    w = sorted(files, key=lambda f: B[f][0] - A[f][0])
    print(f'-- worst for {b} vs {a}:')
    for f in w[:top]: print(f'  {B[f][0]-A[f][0]:+.4f}  {f}  {a}={A[f][0]:.4f} {b}={B[f][0]:.4f}')
    print(f'-- best for {b} vs {a}:')
    for f in w[-top:]: print(f'  {B[f][0]-A[f][0]:+.4f}  {f}  {a}={A[f][0]:.4f} {b}={B[f][0]:.4f}')
if __name__ == '__main__':
    if sys.argv[1] == 'run': run(sys.argv[2], sys.argv[3], sys.argv[4:])
    else:
        top = 15; split = 'all'
        if '--top' in sys.argv: top = int(sys.argv[sys.argv.index('--top') + 1])
        if '--split' in sys.argv: split = sys.argv[sys.argv.index('--split') + 1]
        cmp(sys.argv[2], sys.argv[3], top, split)
