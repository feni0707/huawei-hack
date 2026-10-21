#!/usr/bin/env python3
"""kn_tune.py test.in out.kn --from F [-j J] [--iters N] [--resume] : подбор ручек движка (mk_sfkn.py) под один тест.
(1+J)-стратегия: J кандидатов — текущий вектор с 1–3 изменёнными ручками (гаусс, масштаб ~30 % значения), каждый — полный
прогон (dsearch, KNEVAL; ручки действуют с команды F — там посылка узнаёт тест); лучший принимается при строгом улучшении
(T, emax, площадь). out.kn — строки "индекс q" (ручка = q/1e4), только отличные от значений по умолчанию.
Классы — через env (CLSK/CLSB/CLSTAB), как у dsearch. Нужен build/bin/dsearch4.exe (env DSEXE)."""
import os, sys, subprocess, tempfile, random, time, concurrent.futures as cf
args = sys.argv[1:]
test, outp = args[0], args[1]
def opt(n, d): return type(d)(args[args.index(n) + 1]) if n in args else d
FROM = opt('--from', 0); J = opt('-j', 8); ITERS = opt('--iters', 100000); SEED = opt('--seed', 1)
exe = os.path.abspath(os.environ.get('DSEXE', 'build/bin/dsearch4.exe'))
tmp = tempfile.mkdtemp()
dflt = {}
for l in subprocess.run([exe, test], capture_output=True, text=True, env=dict(os.environ, KNDEF='1')).stdout.splitlines():
    k, v = l.split(); dflt[int(k)] = round(float(v) * 1e4)
# настраиваемые: веса 0..19, 25, 26 (14, 20..24 не используются) и параметры участков/зон 27..42
TUNE = [k for k in range(0, 20) if k != 14] + [25, 26] + list(range(27, 43))
INT = set(range(27, 34)) | {35, 37, 42}
LIM = {33: (2, 64), 34: (0, 1), 36: (0, 1), 38: (0, 1), 40: (0, 1), 41: (0, 1)}
cur = dict(dflt)
if '--resume' in args and os.path.exists(outp):
    for l in open(outp):
        k, q = map(int, l.split()); cur[k] = q
def write(vec, path):
    with open(path, 'w') as f:
        for k in sorted(vec):
            if vec[k] != dflt[k]: f.write(f'{k} {vec[k]}\n')
def evalv(vec, tag):
    p = os.path.join(tmp, f'k{tag}.kn'); write(vec, p)
    out = subprocess.run([exe, test], capture_output=True, text=True, env=dict(os.environ, KNEVAL='1', KNFILE=p, KNFROM=str(FROM))).stdout
    for l in out.splitlines():
        if l.startswith('K '):
            q = l.split(); return (int(q[1]), int(q[2]), float(q[3]))
    return (1 << 60, 0, 0)
rng = random.Random(SEED)
def mutate(vec):
    v = dict(vec)
    for _ in range(rng.choice([1, 1, 2, 3])):
        k = rng.choice(TUNE)
        x = v[k] / 1e4
        sc = max(abs(dflt[k] / 1e4), abs(x), 1.0) * 0.3
        x += rng.gauss(0, sc)
        if k in LIM: x = min(max(x, LIM[k][0]), LIM[k][1])
        if k in INT: x = max(1, round(x))
        v[k] = round(x * 1e4)
    return v
base = evalv(cur, 'base')
print('start', base, 'changed', sum(1 for k in cur if cur[k] != dflt[k]), flush=True)
t0 = time.time()
with cf.ThreadPoolExecutor(J) as ex:
    for it in range(ITERS):
        cands = [mutate(cur) for _ in range(J)]
        res = list(ex.map(lambda kv: evalv(kv[1], kv[0]), enumerate(cands)))
        bi = min(range(J), key=lambda i: res[i])
        if res[bi] < base:
            base = res[bi]; cur = cands[bi]; write(cur, outp)
            ch = {k: cur[k] / 1e4 for k in cur if cur[k] != dflt[k]}
            print(f'it {it} {time.time() - t0:.0f}s: {base} changed {ch}', flush=True)
        elif it % 20 == 0:
            print(f'it {it} {time.time() - t0:.0f}s: {base}', flush=True)
print('final', base)
