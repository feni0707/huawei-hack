#!/usr/bin/env python3
"""Тесты-«остатки»: префикс открытого теста + чужое продолжение (ловушка для решений, подогнанных под открытые тесты, и
проверка, насколько гибкой остаётся раскладка после типичного начала).
Продолжения: (a) массовое удаление доли q живых и укладка длин из набора до прежнего пика (UB = 1);
(b) удаление всех труб одной чётности PID и спрос на длину 2·L; (c) смена распределения длин на случайное с удалениями.
  gen_residual.py --seed S --count K --out DIR   -> DIR/residual/<k>.in"""
import argparse, os, random

def load(t):
    f = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'materials', 'tests', f'{t}.in')).read().split()
    n = int(f[1]); return list(map(int, f[2:2 + n]))

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--seed', type=int, default=1); ap.add_argument('--count', type=int, default=50)
    ap.add_argument('--out', default='build/residual_gen'); ap.add_argument('--tests', default='')   # --tests 5,9: только эти источники
    a = ap.parse_args()
    d = os.path.join(a.out, 'residual'); os.makedirs(d, exist_ok=True)
    rng = random.Random(a.seed)
    opens = {t: load(t) for t in (map(int, a.tests.split(',')) if a.tests else range(5, 16))}
    for k in range(a.count):
        t = rng.choice(list(opens)); src = opens[t]
        p = rng.randint(1, len(src) - 1)
        cmds = src[:p]
        ln = [0]; alive = []; pos = {}
        def place(L):
            ln.append(L); pid = len(ln) - 1; pos[pid] = len(alive); alive.append(pid); cmds.append(L)
        def remove(q):
            i = pos.pop(q); last = alive.pop()
            if last != q: alive[i] = last; pos[last] = i
            cmds.append(-q)
        for x in src[:p]:
            if x > 0: ln.append(x); pos[len(ln) - 1] = len(alive); alive.append(len(ln) - 1)
            else:
                q = -x; i = pos.pop(q); last = alive.pop()
                if last != q: alive[i] = last; pos[last] = i
        live = sum(ln[q] for q in alive)
        peak = live
        kind = rng.choice('aabc')
        budget = 200000 - len(cmds) - 2
        if kind == 'a':
            qf = rng.uniform(0.5, 1.0)
            for q in rng.sample(alive, int(len(alive) * qf)):
                if budget <= 0: break
                remove(q); live -= ln[q]; budget -= 1
            sizes = rng.choice([[40], [20, 40], [rng.randint(8, 40)], [rng.randint(1, 40) for _ in range(3)]])
            while live < peak and budget > 0:
                L = rng.choice(sizes); L = min(L, max(1, peak - live)) if rng.random() < 0.5 else L
                place(L); live += L; budget -= 1
        elif kind == 'b':
            par = rng.randint(0, 1)
            vict = [q for q in alive if q % 2 == par]
            for q in vict:
                if budget <= 0: break
                remove(q); live -= ln[q]; budget -= 1
            for _ in range(min(budget, rng.randint(50, 3000))):
                L = min(40, 2 * rng.choice([ln[q] for q in vict] or [20])); place(L); live += L
        else:
            for _ in range(min(budget // 2, rng.randint(200, 20000))):
                if alive and rng.random() < 0.45:
                    q = rng.choice(alive) if rng.random() < 0.5 else alive[-1]
                    remove(q); live -= ln[q]
                else:
                    L = rng.randint(1, 40); place(L); live += L
        if not alive: place(1)
        with open(os.path.join(d, f'{k}.in'), 'w') as f:
            f.write(f'1 {len(cmds)}\n' + '\n'.join(map(str, cmds)) + '\n0\n')
    print('done', d)
main()
