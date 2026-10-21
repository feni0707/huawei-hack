#!/usr/bin/env python3
"""ov_rd.py test.in out.ov [-j J] [--steps S] [--resume] [--late F] : случайный спуск по переопределениям решений.
Шаг: J случайных кандидатов (pid, вариант) параллельно, каждый — полный прогон (OVSEXE, MODE=eval); лучший принимается,
если (T, emax) строго меньше. Для тестов, где почти любая правка что-то меняет (15), это дешевле, чем выбирать лучшего
из большой выборки (ov_search.py). --late F: доля кандидатов из последней трети решений (там прогон короче не станет,
но влияние на итог прямее). Результат: строки "pid idx"."""
import sys, os, subprocess, tempfile, random, time, concurrent.futures as cf
args = sys.argv[1:]
test, outp = args[0], args[1]
def opt(n, d): return type(d)(args[args.index(n) + 1]) if n in args else d
J = opt('-j', 12); STEPS = opt('--steps', 100000); SEED = opt('--seed', 1)
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
def eval1(d, cand, k):
    p = os.path.join(tmp, f'b{k}.ov'); write_ov(d, p)
    cp = os.path.join(tmp, f'c{k}.txt')
    with open(cp, 'w') as f: f.write(f'{cand[0]} {cand[1]}\n')
    out = subprocess.run([exe, test], capture_output=True, text=True, env=dict(os.environ, MODE='eval', OVFILE=p, CANDS=cp)).stdout
    for l in out.splitlines():
        if l.startswith('E '): q = list(map(int, l.split()[1:])); return (q[2], q[3])
    return None
rng = random.Random(SEED)
dec, base = listdec(ov)
print('start', base, 'decisions', len(dec), 'overrides', len(ov), flush=True)
t0 = time.time(); acc = 0
with cf.ThreadPoolExecutor(J) as ex:
    for step in range(STEPS):
        cands = []
        while len(cands) < J:
            pid, n, df = rng.choice(dec)
            k = rng.choice([x for x in range(n) if x != df])
            if (pid, k) not in cands: cands.append((pid, k))
        snap = dict(ov)
        res = list(ex.map(lambda ik: eval1(snap, ik[1], ik[0]), enumerate(cands)))
        best = None
        for c, r in zip(cands, res):
            if r is not None and r < base and (best is None or r < best[1]): best = (c, r)
        if best:
            prev = dict(ov)
            ov[best[0][0]] = best[0][1]
            dec, chk = listdec(ov)
            if chk != best[1]: print('  расхождение проверки', best, chk, flush=True)
            if chk < base:
                base = chk; acc += 1; write_ov(ov, outp)
            else:
                ov = prev; dec, base = listdec(ov)
        if step % 10 == 0 or best:
            print(f'step {step}: base {base}, accepted {acc}, overrides {len(ov)}, {time.time() - t0:.0f} s', flush=True)
write_ov(ov, outp)
print('final', base, 'overrides', len(ov))
