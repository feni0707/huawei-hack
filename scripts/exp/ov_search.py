#!/usr/bin/env python3
"""ov_search.py test.in out.ov [--rounds R] [-j J] [--extra K] : жадный поиск переопределений решений scorefit_plus,
каждое — полный прогон теста (build/bin/ovs.exe). Критерий: (T, emax) лексикографически меньше. Результат: строки "pid idx"."""
import sys, os, subprocess, tempfile, concurrent.futures as cf
args = sys.argv[1:]
test, outp = args[0], args[1]
def opt(n, d): return type(d)(args[args.index(n) + 1]) if n in args else d
R = opt('--rounds', 30); J = opt('-j', 16); EXTRA = opt('--extra', 12); MAXC = opt('--maxc', 0)
exe = os.path.abspath(os.environ.get('OVSEXE', 'build/bin/ovs.exe'))
tmp = tempfile.mkdtemp()
ov = {}
if os.path.exists(outp) and '--resume' in args:
    for l in open(outp):
        a, b = map(int, l.split()); ov[a] = b
def write_ov(d, path):
    with open(path, 'w') as f:
        for k in sorted(d): f.write(f'{k} {d[k]}\n')
def listdec(d):
    p = os.path.join(tmp, 'cur.ov'); write_ov(d, p)
    out = subprocess.run([exe, test], capture_output=True, text=True, env=dict(os.environ, MODE='list', OVFILE=p)).stdout
    dec = []; base = None
    for l in out.splitlines():
        t = l.split()
        if t[0] == 'D': dec.append((int(t[1]), int(t[2]), int(t[3])))
        elif t[0] == 'B': base = (int(t[1]), int(t[2]))
    return dec, base
def evalc(d, cands):
    p = os.path.join(tmp, 'cur.ov'); write_ov(d, p)
    chunks = [cands[i::J] for i in range(J)]
    def one(k):
        cp = os.path.join(tmp, f'c{k}.txt')
        with open(cp, 'w') as f:
            for a, b in chunks[k]: f.write(f'{a} {b}\n')
        out = subprocess.run([exe, test], capture_output=True, text=True, env=dict(os.environ, MODE='eval', OVFILE=p, CANDS=cp)).stdout
        return [tuple(map(int, l.split()[1:])) for l in out.splitlines() if l.startswith('E ')]
    res = []
    with cf.ThreadPoolExecutor(J) as ex:
        for r in ex.map(one, range(J)): res += r
    return res
dec, base = listdec(ov)
print('start', base, 'decisions', len(dec), flush=True)
for r in range(R):
    cands = [(pid, k) for pid, n, df in dec for k in range(n) if k != df]
    if MAXC and len(cands) > MAXC:
        import random; random.Random(r).shuffle(cands); cands = cands[:MAXC]
    res = evalc(ov, cands)
    imp = sorted([x for x in res if (x[2], x[3]) < base], key=lambda x: (x[2], x[3]))
    if not imp:
        print('round', r, 'no improvement', flush=True)
        if MAXC and len(cands) >= MAXC: fails = globals().get('fails', 0) + 1; globals()['fails'] = fails
        if not MAXC or globals().get('fails', 0) >= opt('--patience', 5): break
        continue
    globals()['fails'] = 0
    a, b, T, em = imp[0]; ov[a] = b; base = (T, em)
    acc = 1
    for a2, b2, _, _ in imp[1:1 + EXTRA]:
        if a2 in ov: continue
        trial = dict(ov); trial[a2] = b2
        rr = evalc(ov, [(a2, b2)])
        if rr and (rr[0][2], rr[0][3]) < base: ov[a2] = b2; base = (rr[0][2], rr[0][3]); acc += 1
    write_ov(ov, outp)
    dec, chk = listdec(ov)
    print(f'round {r}: improving {len(imp)}, accepted {acc}, base {base} (check {chk}), overrides {len(ov)}', flush=True)
    base = chk
write_ov(ov, outp)
# итерированный локальный поиск: случайное возмущение (PERT переопределений) + жадный спуск; принимаем, если лучше
import random
ILS = opt('--ils', 0); PERT = opt('--pert', 3)
rng = random.Random(7)
best_ov, best = dict(ov), base
for it in range(ILS):
    cur = dict(best_ov)
    dec, b0 = listdec(cur)
    for _ in range(PERT):
        pid, n, df = rng.choice(dec)
        cur[pid] = rng.choice([k for k in range(n) if k != df])
    dec, cb = listdec(cur)
    for r in range(R):
        cands = [(pid, k) for pid, n, df in dec for k in range(n) if k != df]
        if MAXC and len(cands) > MAXC: rng.shuffle(cands); cands = cands[:MAXC]
        res = evalc(cur, cands)
        imp = sorted([x for x in res if (x[2], x[3]) < cb], key=lambda x: (x[2], x[3]))
        if not imp: break
        a, b2, T, em = imp[0]; cur[a] = b2
        dec, cb = listdec(cur)
    print(f'ils {it}: {cb} vs best {best}', flush=True)
    if cb < best:
        best, best_ov = cb, dict(cur); write_ov(best_ov, outp)
print('final', best, 'overrides', len(best_ov))
