#!/usr/bin/env python3
"""ov_hole.py test.in cur.d --minpid P [-j J] [--rounds R] [--perhole K] [--maxc C] : прицельный поиск переопределений позиций.
Раунд: дамп текущей раскладки (dsearch, DUMP) -> момент te наибольшего начала хвоста -> дыры ниже хвоста в момент te и время,
с которого каждая свободна непрерывно (t_h) -> кандидаты: трубы, рождённые в (t_h, te], живые в te, не длиннее дыры
(сначала легшие в хвост, длинные раньше), позиции — края дыры -> оценка параллельно (dsearch, CANDD) -> лучший строго
улучшающий (T, emax, площадь над THR) принимается и дописывается в cur.d. Классы и THR — через env, как у dsearch.
Нужен build/bin/dsearch4.exe (или env DSEXE); ручки — KNFILE/KNFROM."""
import os, sys, subprocess, tempfile, heapq, time, concurrent.futures as cf
args = sys.argv[1:]
test, curp = args[0], args[1]
def opt(n, d): return type(d)(args[args.index(n) + 1]) if n in args else d
MINPID = opt('--minpid', 1); J = opt('-j', 8); R = opt('--rounds', 1000); PERHOLE = opt('--perhole', 3); MAXC = opt('--maxc', 400)
exe = os.path.abspath(os.environ.get('DSEXE', 'build/bin/dsearch4.exe'))
tmp = tempfile.mkdtemp()
tok = open(test).read().split(); n = int(tok[1]); cmds = list(map(int, tok[2:2 + n]))
L = [0]; birth = [0]; death = [1 << 30]
for q, x in enumerate(cmds):
    if x > 0: L.append(x); birth.append(q); death.append(1 << 30)
    elif x < 0: death[-x] = q
P = len(L) - 1

def run(env_extra):
    return subprocess.run([exe, test], capture_output=True, text=True, env=dict(os.environ, MINPID=str(MINPID), **env_extra)).stdout

def dump():
    dp = os.path.join(tmp, 'dump.txt'); op = os.path.join(tmp, 'pruned.d')
    out = run({'SECS': '0', 'DFILE': curp, 'OUT': op, 'DUMP': dp})
    fin = [l for l in out.splitlines() if l.startswith('final')][0].split()
    base = (int(fin[2]), int(fin[4]), float(fin[6]))
    D = [0] * (P + 1)
    for l in open(dp):
        p, d, _ = map(int, l.split()); D[p] = d
    os.replace(op, curp)   # прунинг: несработавшие переопределения выброшены
    return base, D

def analyze(D):
    # начало хвоста по времени (максимум концов живых труб), te — первый максимум
    hp = []; ts = []; cur = 0; pid = 0; best = -1; te = 0
    tsBefore = [0] * (P + 1)
    for q, x in enumerate(cmds):
        if x > 0:
            pid += 1; tsBefore[pid] = hp and -hp[0][0] or 0
            while hp and death[hp[0][1]] <= q - 1 and death[hp[0][1]] != 1 << 30: heapq.heappop(hp)
            tsBefore[pid] = -hp[0][0] if hp else 0
            heapq.heappush(hp, (-(D[pid] + x), pid))
        while hp and death[hp[0][1]] <= q: heapq.heappop(hp)
        e = -hp[0][0] if hp else 0
        if e > best: best = e; te = q
    # занятость и «свободна с» в момент te
    T = best + 64
    lastFree = [0] * (T + 64); occ = [0] * (T + 64)
    pid = 0
    for q, x in enumerate(cmds[:te + 1]):
        if x > 0:
            pid += 1
            for c in range(D[pid], D[pid] + x): occ[c] = pid
        else:
            p = -x
            for c in range(D[p], D[p] + L[p]): occ[c] = 0; lastFree[c] = q
    holes = []; c = 0
    while c < best:
        if occ[c]: c += 1; continue
        e = c
        while e < best and not occ[e]: e += 1
        holes.append((c, e - c, max(lastFree[c:e])))
        c = e
    alive_te = [p for p in range(MINPID, P + 1) if birth[p] <= te < death[p]]
    return te, best, holes, alive_te, tsBefore

def candidates(D, te, holes, alive_te, tsBefore):
    byLen = {}
    for p in alive_te: byLen.setdefault(L[p], []).append(p)
    cands = []
    for a, s, th in holes:
        pool = []
        for l in range(min(s, 40), 0, -1):
            for p in byLen.get(l, []):
                if birth[p] > th and not (a <= D[p] < a + s):
                    tail = D[p] >= tsBefore[p]
                    pool.append((0 if tail else 1, -l, p))
        pool.sort()
        for _, _, p in pool[:PERHOLE]:
            cands.append((0 if D[p] >= tsBefore[p] else 1, -s, p, a))
            if s > L[p]: cands.append((0 if D[p] >= tsBefore[p] else 1, -s, p, a + s - L[p]))
    cands.sort()
    seen = set(); res = []
    for _, _, p, d in cands:
        if (p, d) in seen: continue
        seen.add((p, d)); res.append((p, d))
    return res

def evalc(cands):
    chunks = [cands[i::J] for i in range(J)]
    def one(k):
        if not chunks[k]: return []
        cp = os.path.join(tmp, f'c{k}.txt')
        with open(cp, 'w') as f:
            for p, d in chunks[k]: f.write(f'{p} {d}\n')
        out = run({'CANDD': cp, 'DFILE': curp})
        r = []
        for l in out.splitlines():
            if l.startswith('E '):
                q = l.split(); r.append((int(q[1]), int(q[2]), int(q[3]), int(q[4]), int(q[5]), int(q[6]), float(q[7])))
        return r
    res = []
    with cf.ThreadPoolExecutor(J) as ex:
        for r in ex.map(one, range(J)): res += r
    return res

t0 = time.time()
tried = set()
for rd in range(R):
    base, D = dump()
    te, emax, holes, alive_te, tsb = analyze(D)
    allc = [c for c in candidates(D, te, holes, alive_te, tsb) if (c, base) not in tried]
    cands = allc[:MAXC]
    print(f'round {rd}: base {base}, te {te}, holes {len(holes)} ({sum(h[1] for h in holes)} m), candidates {len(cands)} of {len(allc)}, {time.time() - t0:.0f} s', flush=True)
    if not cands: break
    for c in cands: tried.add((c, base))
    res = evalc(cands)
    ok = [r for r in res if r[4] >= 0 and (r[4], r[5], r[6]) < base]
    if not ok:
        print('  нет улучшений', flush=True); continue   # следующая порция кандидатов при той же раскладке
    ok.sort(key=lambda r: (r[4], r[5], r[6]))
    p, d, gs, gl, T, em, ar = ok[0]
    lines = [l for l in open(curp) if l.strip() and int(l.split()[0]) != p]
    lines.append(f'{p} {d} {gs} {gl}\n')
    lines.sort(key=lambda l: int(l.split()[0]))
    open(curp, 'w').writelines(lines)
    print(f'  принято pid {p} -> {d}: ({T}, {em}, {ar:.0f}); улучшающих {len(ok)}', flush=True)
base, D = dump()
print('final', base)
