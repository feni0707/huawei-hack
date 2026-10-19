"""Ещё классы «на устойчивость» (нет в открытых тестах; ищем, где scorefit хуже best fit):
python3 scripts/gen_robust3.py [--seed S] [--count K] [--out DIR] -> DIR/<class>/<i>.in (по умолчанию build/robust3)

  fifo        — удаляется самая старая живая (очередь), длины из случайной смеси
  uniform_del — удаляется случайная живая (равномерно)
  size_life   — время жизни зависит от длины (короткие живут мало / долго — случайно для теста)
  sawtooth    — набрать K труб, удалить случайные 60–95 %, повторить с другим распределением длин
  desc_groups — группы равных длин по убыванию 40 → 1, удаляется каждая вторая / третья
  pairs       — пары (длинная, короткая), позже удаляется один из элементов пары (всегда один и тот же)
  regimes     — марковская смена режимов длин (диапазоны), удаляются свежие
  queue_mix   — удаления то из начала (FIFO), то из конца (LIFO)
  rev_two_phase, stack_swings, parity_flip, triplets, two_lengths, del_half_group,
  walk_len, bimodal_life, pid_mod, len_class_death, old_window, long_then_ones — см. комментарии функций
"""
import argparse, os, random, sys
import importlib.util
spec = importlib.util.spec_from_file_location("g", os.path.join(os.path.dirname(__file__), "gen_tests.py"))
g = importlib.util.module_from_spec(spec); spec.loader.exec_module(g)
B = g.Builder


def len_sampler(r):
    kind = r.randrange(4)
    if kind == 0:
        lo, hi = sorted(r.sample(range(1, 41), 2)); return lambda: r.randint(lo, hi)
    if kind == 1:
        vals = r.sample(range(1, 41), r.randint(2, 6)); return lambda: r.choice(vals)
    if kind == 2:
        return lambda: min(40, 1 + int(r.expovariate(1 / r.uniform(2, 12))))
    return lambda: r.randint(1, 40)


def fifo(r, n):
    b = B(n); ln = len_sampler(r); p = r.uniform(0.3, 0.6)
    while b.room():
        b.place(ln())
        while b.room() and len(b.alive) > 1 and r.random() < p:
            b.remove(min(b.alive))
    return b.cmds


def uniform_del(r, n):
    b = B(n); ln = len_sampler(r); p = r.uniform(0.3, 0.6)
    while b.room():
        b.place(ln())
        while b.room() and len(b.alive) > 1 and r.random() < p:
            b.remove(r.choice(b.alive))
    return b.cmds


def size_life(r, n):
    # время жизни (в укладках) ~ Exp(mean(L)); mean растёт или падает с длиной
    b = B(n); ln = len_sampler(r); inc = r.random() < 0.5
    import heapq
    q = []; t = 0
    while b.room():
        L = ln(); p = b.place(L); t += 1
        m = (5 + 40 * L) if inc else (5 + 40 * (41 - L))
        heapq.heappush(q, (t + r.expovariate(1 / m) * r.uniform(0.6, 1.0), p))
        while q and q[0][0] <= t and b.room():
            _, x = heapq.heappop(q)
            if b.is_alive(x): b.remove(x)
    return b.cmds


def sawtooth(r, n):
    b = B(n)
    while b.room(10):
        ln = len_sampler(r); K = r.randint(50, 1500)
        for _ in range(K):
            if not b.room(2): break
            b.place(ln())
        frac = r.uniform(0.6, 0.95)
        for x in r.sample(b.alive, int(len(b.alive) * frac)):
            if not b.room(): break
            b.remove(x)
    return b.cmds


def desc_groups(r, n):
    b = B(n); step = r.choice([2, 3]); L = 40
    while b.room(10) and L >= 1:
        k = r.randint(60, 300); pids = []
        for _ in range(k):
            if not b.room(2): break
            pids.append(b.place(L))
        off = r.randrange(step)
        for i, x in enumerate(pids):
            if i % step == off and b.room(): b.remove(x)
        L -= r.randint(1, 3)
    return b.cmds


def pairs(r, n):
    b = B(n); first = r.random() < 0.5; ln = len_sampler(r); pend = []
    while b.room(3):
        a = b.place(ln()); c = b.place(r.randint(1, 5))
        pend.append(a if first else c)
        while pend and b.room() and r.random() < 0.45:
            x = pend.pop(max(0, len(pend) - 1 - int(r.expovariate(1 / 20))))
            if b.is_alive(x): b.remove(x)
    return b.cmds


def regimes(r, n):
    b = B(n); ranges = [sorted(r.sample(range(1, 41), 2)) for _ in range(3)]; cur = 0
    while b.room():
        if r.random() < 0.002: cur = r.randrange(3)
        lo, hi = ranges[cur]; b.place(r.randint(lo, hi))
        while b.room() and len(b.alive) > 1 and r.random() < 0.42:
            b.remove_recent(r, 30)
    return b.cmds


def queue_mix(r, n):
    b = B(n); ln = len_sampler(r); pf = r.uniform(0.2, 0.8)
    while b.room():
        b.place(ln())
        while b.room() and len(b.alive) > 1 and r.random() < 0.45:
            b.remove(min(b.alive) if r.random() < pf else max(b.alive))
    return b.cmds


CLASSES = dict(fifo=fifo, uniform_del=uniform_del, size_life=size_life, sawtooth=sawtooth,
               desc_groups=desc_groups, pairs=pairs, regimes=regimes, queue_mix=queue_mix)


def rev_two_phase(r, n):
    # зеркало two_phase: сначала длинные (удаляются свежие), потом только короткие, иногда умирают старые длинные
    b = B(n); a = r.randint(5, 8); split = int(n * r.uniform(0.4, 0.6))
    while b.room():
        k = len(b.cmds)
        b.place(r.randint(a + 1, 40) if k < split else r.randint(1, a))
        while b.room() and len(b.alive) > 1 and r.random() < 0.45:
            if k >= split and r.random() < 0.2: b.remove(r.choice(b.alive))
            else: b.remove_recent(r, 30)
    return b.cmds


def stack_swings(r, n):
    # стек: набрать много, снять много с конца (большие качели живой массы)
    b = B(n); ln = len_sampler(r)
    while b.room(10):
        for _ in range(r.randint(20, 800)):
            if not b.room(2): break
            b.place(ln())
        for _ in range(int(len(b.alive) * r.uniform(0.3, 0.9))):
            if not b.room() or not b.alive: break
            b.remove(max(b.alive))
    return b.cmds


def parity_flip(r, n):
    # группы равных длин; умирает каждая вторая, но какая чётность — меняется от группы к группе случайно
    b = B(n); L = r.randint(1, 10)
    while b.room(300):
        pids = [b.place(L) for _ in range(r.randint(50, 250))]
        off = r.randrange(2)
        for i, x in enumerate(pids):
            if i % 2 == off and b.room(): b.remove(x)
        L = max(1, min(40, L + r.randint(-3, 5)))
    return b.cmds


def triplets(r, n):
    # тройки (a, b, c); позже умирает средняя
    b = B(n); ln = len_sampler(r); pend = []
    while b.room(4):
        b.place(ln()); m = b.place(ln()); b.place(ln())
        pend.append(m)
        while pend and b.room() and r.random() < 0.4:
            x = pend.pop(max(0, len(pend) - 1 - int(r.expovariate(1 / 30))))
            if b.is_alive(x): b.remove(x)
    return b.cmds


def two_lengths(r, n):
    # две близкие длины (L, L+1): точные попадания редки
    b = B(n); L = r.randint(3, 30)
    while b.room():
        b.place(L + (r.random() < 0.5))
        while b.room() and len(b.alive) > 1 and r.random() < 0.45:
            b.remove_recent(r, r.choice([5, 30, 200]))
    return b.cmds


def del_half_group(r, n):
    # группы равных длин, удаляется первая или вторая половина группы целиком
    b = B(n); L = r.randint(1, 5)
    while b.room(300):
        pids = [b.place(L) for _ in range(r.randint(40, 250))]
        h = len(pids) // 2
        for x in (pids[:h] if r.random() < 0.5 else pids[h:]):
            if b.room(): b.remove(x)
        L = min(40, L + r.randint(0, 3))
    return b.cmds


CLASSES.update(rev_two_phase=rev_two_phase, stack_swings=stack_swings, parity_flip=parity_flip,
               triplets=triplets, two_lengths=two_lengths, del_half_group=del_half_group)


def walk_len(r, n):
    # длина — случайное блуждание ±1..3, удаляются свежие и изредка старые
    b = B(n); L = r.randint(1, 40)
    while b.room():
        L = max(1, min(40, L + r.randint(-3, 3))); b.place(L)
        while b.room() and len(b.alive) > 1 and r.random() < 0.43:
            b.remove_recent(r, 20) if r.random() < 0.85 else b.remove(r.choice(b.alive))
    return b.cmds


def bimodal_life(r, n):
    # каждая труба с вероятностью p умирает почти сразу (1–5 команд), иначе живёт вечно; признаков нет
    b = B(n); ln = len_sampler(r); p = r.uniform(0.3, 0.7); pend = []
    while b.room():
        x = b.place(ln())
        if r.random() < p: pend.append((len(b.cmds) + r.randint(1, 5), x))
        while pend and pend[0][0] <= len(b.cmds) and b.room():
            _, y = pend.pop(0)
            if b.is_alive(y): b.remove(y)
    return b.cmds


def pid_mod(r, n):
    # группы, удаляется каждая k-я (k = 4 или 5) — признака mod k у потоков нет
    b = B(n); k = r.choice([4, 5]); L = r.randint(1, 6)
    while b.room(300):
        pids = [b.place(L) for _ in range(r.randint(60, 300))]
        off = r.randrange(k)
        for i, x in enumerate(pids):
            if i % k == off and b.room(): b.remove(x)
        L = min(40, L + r.randint(0, 3))
    return b.cmds


def len_class_death(r, n):
    # набор длин меняется по эпохам; при смене эпохи трубы одной старой длины постепенно вымирают
    b = B(n); lens = r.sample(range(1, 41), 4); by = {}
    while b.room(10):
        for _ in range(r.randint(300, 1500)):
            if not b.room(2): break
            L = r.choice(lens); by.setdefault(L, []).append(b.place(L))
            if r.random() < 0.3:
                v = r.choice(lens)
                if by.get(v):
                    x = by[v].pop(r.randrange(len(by[v])))
                    if b.is_alive(x) and b.room(): b.remove(x)
        dead = lens.pop(r.randrange(len(lens)))
        for x in by.pop(dead, []):
            if b.is_alive(x) and b.room() and r.random() < 0.9: b.remove(x)
        lens.append(r.choice([v for v in range(1, 41) if v not in lens]))
    return b.cmds


def old_window(r, n):
    # все длины; удаления: свежие или случайные из окна в глубине (как 2-я фаза 15, но без смены длин)
    b = B(n); ln = len_sampler(r); pw = r.uniform(0.1, 0.4)
    while b.room():
        b.place(ln())
        while b.room() and len(b.alive) > 2 and r.random() < 0.45:
            if r.random() < pw:
                k = len(b.alive); srt = sorted(b.alive)
                b.remove(srt[int(k * r.uniform(0.5, 0.9))])
            else: b.remove_recent(r, 30)
    return b.cmds


def long_then_ones(r, n):
    # длинные, случайная половина умирает, дальше только единички (и немного удалений)
    b = B(n); split = int(n * r.uniform(0.3, 0.5))
    while b.room() and len(b.cmds) < split:
        b.place(r.randint(10, 40))
    for x in r.sample(b.alive, len(b.alive) // 2):
        if b.room(): b.remove(x)
    while b.room():
        b.place(1)
        if r.random() < 0.2 and b.room() and len(b.alive) > 1: b.remove_recent(r, 10)
    return b.cmds


CLASSES.update(walk_len=walk_len, bimodal_life=bimodal_life, pid_mod=pid_mod, len_class_death=len_class_death,
               old_window=old_window, long_then_ones=long_then_ones)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--count", type=int, default=12)
    ap.add_argument("--out", default="build/robust3")
    ap.add_argument("--class", dest="cls", default=None)
    a = ap.parse_args()
    sizes = [2000, 10000, 30000]
    for name, fn in CLASSES.items():
        if a.cls and name != a.cls: continue
        d = os.path.join(a.out, name); os.makedirs(d, exist_ok=True)
        for i in range(a.count):
            r = random.Random(a.seed * 1000003 + i * 7919 + sum(map(ord, name)) * 131)
            cmds = fn(r, sizes[i % len(sizes)])
            g.validate(cmds)
            g.write_test(os.path.join(d, f"{i + 1}.in"), cmds, 1)
    print("ok")


if __name__ == "__main__":
    main()
