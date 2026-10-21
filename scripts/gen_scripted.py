#!/usr/bin/env python3
"""Генератор семейства «scripted»: сценарные тесты в духе открытых 8/12/13/14 (механика, не константы).
Тест = последовательность фаз; фазы:
  ladder   — лестница длин L1<L2<..: для каждой ~10-25 циклов «+L -L(только что положенная)», изредка +L +L (запас растёт)
  pair     — чередование двух-трёх длин: +a -b +b -a ... ; удаляется живая своей длины с малым рангом от конца (0..5), изредка старая
  ramp     — рост: случайные длины (равномерно/узкая полоса/кластеры), редкие удаления случайной (с уклоном к недавним) живой
  ones     — длинная серия +1 (с примесью 2,3), потом всплески удаления/добавления «подобных» длин
  swap     — стационар: удалить живую длины из полосы (ранг от конца/возраст случайны) и добавить близкую длину
  burst    — рост большими длинами 15..26 и в конце «замены» старых
Использование: gen_scripted.py --out DIR --count 40 --seed 1 [--min 1500 --max 25000]
Формат: «вес N» затем N команд (0 в конце не пишется — как в открытых: N = число команд без 0? см. add_end)."""
import argparse, os, random

class G:
    def __init__(s, rnd, limit):
        s.r = rnd; s.limit = limit; s.cmds = []; s.plen = [0]; s.live = {}  # L -> list pid (по возрастанию)
        s.alive = 0
    def full(s): return len(s.cmds) >= s.limit
    def add(s, L):
        if s.full(): return None
        s.plen.append(L); pid = len(s.plen) - 1
        s.live.setdefault(L, []).append(pid); s.cmds.append(L); s.alive += 1; return pid
    def rem(s, pid):
        if s.full(): return
        L = s.plen[pid]; lst = s.live[L]; lst.remove(pid) if lst[-1] != pid else lst.pop()
        s.cmds.append(-pid); s.alive -= 1
    def rem_rank(s, L, rank):
        lst = s.live.get(L)
        if not lst: return False
        rank = min(rank, len(lst) - 1); s.rem(lst[len(lst) - 1 - rank]); return True
    def rem_random(s, recent_bias=0.0, lo=1, hi=40):
        Ls = [L for L in range(lo, hi + 1) if s.live.get(L)]
        if not Ls: return False
        L = s.r.choice(Ls)
        n = len(s.live[L])
        if s.r.random() < recent_bias: rank = min(n - 1, int(s.r.expovariate(0.35)))
        else: rank = s.r.randrange(n)
        return s.rem_rank(L, rank)
    def rem_any_pid(s, recent_bias):
        # случайная живая по PID: с уклоном к недавним
        alls = [p for lst in s.live.values() for p in lst]
        if not alls: return False
        alls.sort()
        if s.r.random() < recent_bias: i = len(alls) - 1 - min(len(alls) - 1, int(s.r.expovariate(0.05)))
        else: i = s.r.randrange(len(alls))
        s.rem(alls[i]); return True

def geo_rank(r, p_stay=0.6, cap=5):
    k = 0
    while k < cap and r.random() > p_stay: k += 1
    return k

def ph_ladder(g, n):
    r = g.r; L = r.randint(5, 9); step = r.choice([1, 1, 1, 2]); reps = r.randint(10, 25); pextra = r.choice([0.05, 0.1, 0.15])
    end = g.cmds.__len__() + n
    while len(g.cmds) < end and not g.full():
        for _ in range(reps):
            if r.random() < pextra: g.add(L); g.add(L)
            else: g.add(L)
            g.rem_rank(L, 0) if r.random() < 0.9 else None
        L += step
        if L > 39: L = r.randint(5, 9)
        reps = r.randint(10, 25)

def ph_pair(g, n):
    r = g.r; k = r.choice([2, 2, 3]); base = r.choice([5, 5, 6, 7, 8, 9, 10])
    Ls = sorted(r.sample([base, base * 2, base * 3, base + 1, base + 5, base * 2 + 1, 16], k))
    Ls = sorted(set(min(39, x) for x in Ls)); maxrank = r.choice([2, 4, 6]); pold = r.choice([0.0, 0.02, 0.05])
    # кольцо: +L_i -L_{i+1 mod}: «+5 -10 +10 -5» ; стартуем с небольшим запасом
    for L in Ls:
        for _ in range(r.randint(2, 6)): g.add(L)
    end = len(g.cmds) + n; i = 0
    while len(g.cmds) < end and not g.full():
        La = Ls[i % len(Ls)]; Lb = Ls[(i + 1) % len(Ls)]
        g.add(La)
        if r.random() < 0.97:
            if r.random() < pold: g.rem_rank(Lb, r.randint(5, 40))
            else: g.rem_rank(Lb, geo_rank(r, 0.55, maxrank))
        if r.random() < 0.02: g.add(r.choice(Ls))
        i += 1

def draw_len(r, mode, c):
    if mode == 0: return r.randint(5, 39)
    if mode == 1: return max(1, min(40, int(r.gauss(c, 3))))
    if mode == 2: return r.choice(c)
    return r.randint(1, 40)

def ph_ramp(g, n):
    r = g.r; mode = r.choice([0, 1, 1, 2, 2]); c = r.randint(8, 30) if mode == 1 else r.sample([5, 10, 15, 20, 9, 14, 6, 13, 16], 3)
    pdel = r.choice([0.0, 0.05, 0.15, 0.3, 0.45]); rb = r.choice([0.3, 0.7, 0.95]); sty = r.choice(['len', 'pid'])
    tail = r.random() < 0.5
    end = len(g.cmds) + n
    while len(g.cmds) < end and not g.full():
        if r.random() < pdel and g.alive > 20:
            if sty == 'len': g.rem_random(rb)
            else: g.rem_any_pid(rb)
        else:
            if mode == 2 and r.random() < 0.85: g.add(r.choice(c))
            else: g.add(draw_len(r, mode, c))

def ph_ones(g, n):
    r = g.r; end = len(g.cmds) + n; burst = r.random() < 0.6
    nb = int(n * (0.25 if burst else 0.02))
    while len(g.cmds) < end - nb and not g.full():
        x = r.random(); g.add(1 if x < 0.7 else (2 if x < 0.85 else (3 if x < 0.97 else r.randint(4, 12))))
    # всплески удаления / добавления похожих
    while len(g.cmds) < end and not g.full():
        if r.random() < 0.15 and g.alive > 10:
            g.rem_random(0.6, 1, 12)
            if r.random() < 0.8: g.add(r.randint(1, 12))
        else: g.add(1 if r.random() < 0.8 else r.randint(2, 12))

def ph_swap(g, n):
    r = g.r; lo = r.randint(10, 24); hi = lo + r.randint(4, 14); hi = min(hi, 40); end = len(g.cmds) + n
    pgrow = r.choice([0.3, 0.5, 0.7])
    while len(g.cmds) < end and not g.full():
        x = r.random()
        if x < pgrow or g.alive < 30: g.add(r.randint(lo, hi))
        else:
            L = r.randint(lo, hi)
            if g.live.get(L):
                if r.random() < 0.5: g.rem_rank(L, r.randrange(len(g.live[L])))
                else: g.rem_rank(L, geo_rank(r, 0.5, 6))
                g.add(min(40, max(1, L + r.choice([-1, 0, 0, 1]))))

PH = {'ladder': ph_ladder, 'pair': ph_pair, 'ramp': ph_ramp, 'ones': ph_ones, 'swap': ph_swap}

def gen(seed, nmin, nmax):
    r = random.Random(seed); total = int(nmin * (nmax / nmin) ** r.random())
    g = G(r, total); k = r.choice([1, 2, 2, 3, 4]); names = list(PH)
    first = r.choice(['ramp', 'ladder', 'pair', 'ones'])
    seq = [first] + [r.choice(names) for _ in range(k - 1)]
    # доли
    ws = [r.random() + 0.3 for _ in seq]; sw = sum(ws)
    for nm, w in zip(seq, ws):
        if g.full(): break
        PH[nm](g, max(50, int(total * w / sw)))
    while not g.full(): ph_ramp(g, total - len(g.cmds))
    return g.cmds[:total], seq

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', required=True); ap.add_argument('--count', type=int, default=40)
    ap.add_argument('--seed', type=int, default=1); ap.add_argument('--min', type=int, default=1500); ap.add_argument('--max', type=int, default=25000)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    for i in range(a.count):
        cmds, seq = gen(a.seed * 100003 + i, a.min, a.max)
        # проверка корректности
        alive = set(); n = 0
        for c in cmds:
            if c > 0: n += 1; alive.add(n)
            else: assert -c in alive; alive.discard(-c)
        with open(os.path.join(a.out, '%02d.in' % (i + 1)), 'w') as f:
            f.write('4 %d\n' % len(cmds) + ''.join('%d\n' % c for c in cmds) + '0\n')
        print('%02d n=%d phases=%s' % (i + 1, len(cmds), '+'.join(seq)))
main()
