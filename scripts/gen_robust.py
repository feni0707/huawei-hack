"""Классы «на устойчивость», которых нет в открытых тестах: python3 scripts/gen_robust.py -> build/robust/<class>/."""
import random, sys, os
sys.path.insert(0, 'scripts')
import importlib.util
spec = importlib.util.spec_from_file_location("g", "scripts/gen_tests.py"); g = importlib.util.module_from_spec(spec); spec.loader.exec_module(g)
B = g.Builder
def three_phase(r, n):
    b = B(n); cuts = [n//3, 2*n//3]
    while b.room():
        k = len(b.cmds)
        L = r.randint(1,7) if (k < cuts[0] or k >= cuts[1]) else r.randint(8,40)
        b.place(L)
        if b.alive and r.random() < 0.45:
            for _ in range(1 + (r.random() < 0.4)):
                if b.room() and b.alive: b.remove_recent(r, 30) if r.random() < 0.8 else b.remove_old(r)
    return b.cmds
def long_die(r, n):
    b = B(n); longs = []
    while b.room():
        L = r.randint(1,40); p = b.place(L)
        if L > 20: longs.append(p)
        if r.random() < 0.4 and longs and b.room():
            q = longs.pop(r.randrange(len(longs)))
            if b.is_alive(q): b.remove(q)
    return b.cmds
def parity_noise(r, n):
    # группы равных длин, удаляется случайная половина (не чётность) пачкой — ловушка для потоков
    b = B(n); L = 1
    while b.room(200):
        pids = [b.place(L) for _ in range(r.randint(40, 120))]
        dels = r.sample(pids, len(pids)//2)
        for p in sorted(dels): b.remove(p)
        L = min(40, L + r.randint(0, 3))
    return b.cmds
def mixed_groups(r, n):
    # как sorted_groups, но иногда удаляется каждая третья, иногда — случайные
    b = B(n); lens = sorted(r.sample(range(1,41), 12))
    for L in lens:
        if not b.room(700): break
        pids = [b.place(L) for _ in range(r.randint(100, 400))]
        mode = r.choice(['par', 'third', 'rand'])
        if mode == 'par': dels = pids[r.randint(0,1)::2]
        elif mode == 'third': dels = pids[r.randint(0,2)::3]
        else: dels = r.sample(pids, len(pids)//2)
        for p in dels: b.remove(p)
    return b.cmds
for name, fn in [('three_phase', three_phase), ('long_die', long_die), ('parity_noise', parity_noise), ('mixed_groups', mixed_groups)]:
    d = f'build/robust/{name}'; os.makedirs(d, exist_ok=True)
    for i, n in enumerate([10000, 50000, 200000, 20000]):
        r = random.Random(f'{name}{i}')
        c = fn(r, n); g.write_test(f'{d}/{i+1}.in', c, 1)
print('ok')
