#!/usr/bin/env python3
"""Генератор тестов в формате materials/tests/*.in для проверки обобщения (закрытые тесты 16–30 не видны).

Классы повторяют характер открытых тестов (см. docs/TESTS_ANALYSIS.md):
  sorted_groups — группы одинаковых длин по возрастанию, после группы удаляется каждая вторая (как 5, 9)
  alternating   — почти строгое чередование «уложить / удалить недавнюю», доминируют 2–3 длины (как 8, 13, 12)
  random_bursts — случайные длины (короткие 1..7 часты), пачки удалений свежих и старых труб (как 11, 15)
  phases        — длинная фаза коротких труб без удалений, смешанные фазы, в конце длинные трубы (как 14)
  periodic      — периодический паттерн коротких длин, выборочные удаления по смещениям в периоде,
                  затем дозаполнение «уложить, уложить, удалить» (как 10)

Формат: `<weight> <n>`, n команд, завершающий 0. Удаляются только живые PID, длины 1..40, n ≤ 200000.
Детерминированность: один и тот же --seed даёт те же файлы (у каждого теста свой sub-seed).

  python3 scripts/gen_tests.py                                  # все классы, по 5 тестов, build/gen/<class>/
  python3 scripts/gen_tests.py --class random_bursts --size 200000 --count 2
  python3 scripts/gen_tests.py --class alternating --size 10000 --seed 7 --out build/gen_alt
"""
import argparse
import os
import random
import sys

MAX_CMDS = 200_000
MAX_LEN = 40


class Builder:
    """Накапливает команды и следит за живыми PID (удаление — O(1) через swap с последним)."""

    def __init__(self, limit):
        self.cmds = []
        self.limit = limit
        self.pid = 0
        self.alive = []        # живые PID
        self.idx = {}          # PID -> индекс в alive
        self.length = {}       # PID -> длина

    def room(self, k=1):
        return len(self.cmds) + k <= self.limit

    def place(self, L):
        L = max(1, min(MAX_LEN, int(L)))
        self.pid += 1
        self.cmds.append(L)
        self.idx[self.pid] = len(self.alive)
        self.alive.append(self.pid)
        self.length[self.pid] = L
        return self.pid

    def remove(self, p):
        i = self.idx.pop(p)
        last = self.alive.pop()
        if last != p:
            self.alive[i] = last
            self.idx[last] = i
        self.cmds.append(-p)

    def is_alive(self, p):
        return p in self.idx

    def remove_recent(self, rng, mean_age):
        """Удалить живую трубу «недавнего» возраста: PID ≈ последний − Geom(mean_age); иначе случайную живую."""
        if not self.alive:
            return False
        for _ in range(8):
            age = int(rng.expovariate(1.0 / max(mean_age, 0.5)))
            p = self.pid - age
            if p >= 1 and self.is_alive(p):
                self.remove(p)
                return True
        self.remove(rng.choice(self.alive))
        return True

    def remove_old(self, rng):
        if not self.alive:
            return False
        self.remove(rng.choice(self.alive))
        return True


# ---------------------------------------------------------------- классы

def gen_sorted_groups(rng, n):
    variant = rng.choice(["full", "pow2", "subset"])
    if variant == "full":
        lens = list(range(1, rng.choice([20, 30, 40]) + 1))
    elif variant == "pow2":
        lens = [1, 2, 4, 8, 16, 32]
    else:
        lens = sorted(rng.sample(range(1, MAX_LEN + 1), rng.randint(5, 15)))
    group = max(2, int(n / (1.5 * len(lens))))
    b = Builder(n)
    for L in lens:
        if not b.room(group + group // 2):
            break
        pids = [b.place(L) for _ in range(group)]
        parity = rng.randint(0, 1)
        for p in pids[parity::2]:
            b.remove(p)
    return b.cmds


def gen_alternating(rng, n):
    # Доминирующие длины: как в 8/13 (5, 10, 16) или случайные — чтобы не подгоняться под конкретные числа.
    if rng.random() < 0.5:
        dom = [5, 10, 16]
    else:
        dom = sorted(rng.sample(range(4, 21), 3))
    minor = sorted(set(rng.sample(range(3, 25), 4)) - set(dom))
    weights = {L: w for L, w in zip(dom, [0.36, 0.34, 0.15])}
    for L in minor:
        weights[L] = 0.12 / max(1, len(minor))
    choices, wts = list(weights), list(weights.values())
    # Вероятность удаления после укладки: ≈0.9 как 8/13 (мало живых в конце), ≈0.3 как 12.
    p_del = rng.choice([0.9, 0.94, 0.6, 0.32])
    mean_age = rng.choice([2, 4, 8])
    b = Builder(n)
    cur = rng.choices(choices, wts)[0]
    while b.room():
        if b.alive and b.cmds[-1] > 0 and rng.random() < p_del:   # удаления не идут подряд
            if rng.random() < 0.92:
                b.remove_recent(rng, mean_age)
            else:
                b.remove_old(rng)
            continue
        if rng.random() < 0.03:
            L = rng.randint(1, 39)                     # редкая «чужая» длина
        else:
            if rng.random() > 0.9:                     # длины идут сериями, иногда переключаются
                cur = rng.choices(choices, wts)[0]
            L = cur
        b.place(L)
    return b.cmds


def _random_len(rng, spikes):
    r = rng.random()
    if r < 0.68:
        return rng.randint(1, 7)
    if r < 0.75 and spikes:
        return rng.choice(spikes)
    return rng.randint(8, MAX_LEN)


def gen_random_bursts(rng, n):
    spikes = rng.sample(range(8, MAX_LEN + 1), 2)
    p_burst = rng.uniform(0.2, 0.3)                    # ≈ 35 % команд — удаления (как 11/15)
    p_recent = rng.uniform(0.5, 0.7)
    b = Builder(n)
    while b.room():
        b.place(_random_len(rng, spikes))
        if b.alive and rng.random() < p_burst:
            k = 1
            while k < 14 and rng.random() < 0.4:       # размер пачки ~ геометрический, до 14
                k += 1
            for _ in range(k):
                if not b.room() or not b.alive:
                    break
                if rng.random() < p_recent:
                    b.remove_recent(rng, 40)
                else:
                    b.remove_old(rng)
    return b.cmds


def gen_phases(rng, n):
    b = Builder(n)
    n1 = int(n * rng.uniform(0.5, 0.65))               # фаза 1: почти только длина 1, без удалений
    n2 = int(n * rng.uniform(0.2, 0.3))                # фаза 2: смешанные короткие/длинные, редкие удаления
    while b.room() and len(b.cmds) < n1:
        b.place(1 if rng.random() < 0.93 else rng.randint(2, 3))
    while b.room() and len(b.cmds) < n1 + n2:
        if b.alive and rng.random() < 0.06:
            b.remove_recent(rng, 300) if rng.random() < 0.7 else b.remove_old(rng)
        else:
            r = rng.random()
            b.place(rng.randint(1, 3) if r < 0.8 else rng.randint(4, MAX_LEN))
    lo = rng.randint(5, 14)                            # фаза 3: длинные трубы
    while b.room():
        if b.alive and rng.random() < 0.12:
            b.remove_recent(rng, 50) if rng.random() < 0.6 else b.remove_old(rng)
        else:
            b.place(rng.randint(lo, MAX_LEN) if rng.random() < 0.85 else rng.randint(1, 4))
    return b.cmds


def gen_periodic(rng, n):
    # Паттерн коротких длин; в половине случаев сумма периода кратна 20 (как в 10: 3 5 2 5 1 3 1 = 20).
    k = rng.randint(4, 9)
    want20 = rng.random() < 0.5
    while True:
        pat = [rng.randint(1, 5) for _ in range(k)]
        if not want20 or sum(pat) % 20 == 0:
            break
    dels = sorted(rng.sample(range(k), rng.randint(1, max(1, k // 2))))  # смещения, удаляемые в фазе 2
    rest = [i for i in range(k) if i not in dels]
    tail_del = rng.choice(rest) if rest else None                         # смещение, удаляемое в фазе 3
    refill = [rng.randint(1, 5) for _ in range(rng.randint(1, 3))]
    b = Builder(n)
    periods = max(1, int(n * rng.uniform(0.45, 0.55)) // k)
    starts = []
    for _ in range(periods):
        if not b.room(k):
            break
        starts.append(b.pid + 1)
        for L in pat:
            b.place(L)
    for s in starts:                                   # фаза 2: выборочные массовые удаления
        for d in dels:
            if not b.room():
                break
            b.remove(s + d)
    i = 0                                              # фаза 3: «уложить refill…, удалить старую трубу»
    while b.room(len(refill) + 1):
        for L in refill:
            b.place(L)
        if tail_del is not None and i < len(starts):
            b.remove(starts[i] + tail_del)
            i += 1
        else:
            b.remove_old(rng)
    return b.cmds


CLASSES = {
    "sorted_groups": (gen_sorted_groups, [2000, 7500, 9000]),
    "alternating": (gen_alternating, [1500, 6500, 10000]),
    "random_bursts": (gen_random_bursts, [10000, 50000, 200000]),
    "phases": (gen_phases, [10000, 25000]),
    "periodic": (gen_periodic, [3000, 10000]),
}


def validate(cmds):
    assert 1 <= len(cmds) <= MAX_CMDS, len(cmds)
    alive, pid = set(), 0
    for c in cmds:
        if c > 0:
            assert c <= MAX_LEN
            pid += 1
            alive.add(pid)
        else:
            assert c < 0 and -c in alive, c
            alive.remove(-c)


def write_test(path, cmds, weight):
    validate(cmds)
    with open(path, "w") as f:
        f.write(f"{weight:g} {len(cmds)}\n")
        f.write("\n".join(map(str, cmds)))
        f.write("\n0\n")


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--class", dest="cls", default="all", help="класс теста или 'all': " + ", ".join(CLASSES))
    ap.add_argument("--count", type=int, default=5, help="тестов на класс (по умолчанию 5)")
    ap.add_argument("--size", type=int, default=0,
                    help="число команд (≤200000); по умолчанию — по кругу из типичных размеров класса")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--weight", type=float, default=1)
    ap.add_argument("--out", default=os.path.join(root, "build", "gen"),
                    help="корневая папка; тесты пишутся в <out>/<class>/<i>.in")
    ap.add_argument("--clean", action="store_true", help="удалить старые *.in в папке класса перед генерацией")
    args = ap.parse_args()

    names = list(CLASSES) if args.cls == "all" else [args.cls]
    for name in names:
        if name not in CLASSES:
            sys.exit(f"unknown class {name!r}; choose from: {', '.join(CLASSES)}")
    if args.size and not 1 <= args.size <= MAX_CMDS:
        sys.exit(f"--size must be in 1..{MAX_CMDS}")

    for ci, name in enumerate(names):
        fn, sizes = CLASSES[name]
        d = os.path.join(args.out, name)
        os.makedirs(d, exist_ok=True)
        if args.clean:
            for f in os.listdir(d):
                if f.endswith(".in"):
                    os.remove(os.path.join(d, f))
        stats = []
        for i in range(1, args.count + 1):
            size = args.size or sizes[(i - 1) % len(sizes)]
            rng = random.Random(f"{args.seed}:{name}:{i}:{size}")
            cmds = fn(rng, size)
            write_test(os.path.join(d, f"{i}.in"), cmds, args.weight)
            stats.append(len(cmds))
        print(f"{name:14s} {args.count} tests -> {os.path.relpath(d)}  (n: {', '.join(map(str, stats))})")


if __name__ == "__main__":
    main()
