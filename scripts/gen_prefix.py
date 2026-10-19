"""Тесты «общий префикс с открытым тестом, потом расхождение» (как 11/15: один генератор, тот же seed, другой n).
usage: python3 scripts/gen_prefix.py OPEN_TEST.in P N MODE SEED OUT
MODE: phase2 — продолжение как 2-я фаза 11/15 (длины 8..c, удаления по стеку своей длины + «ходок» по коротким);
      same   — продолжение тем же распределением длин, что в префиксе, удаления свежих той же длины;
      cut    — префикс и сразу конец.
"""
import random, sys
src, P, N, mode, seed, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], int(sys.argv[5]), sys.argv[6]
toks = open(src).read().split()
w, n0 = toks[0], int(toks[1])
cmds = [int(x) for x in toks[2:2 + n0]]
rng = random.Random(seed)
pre = cmds[:P]
stacks = {l: [] for l in range(1, 41)}
alive_len = {}
pid = 0
for c in pre:
    if c > 0:
        pid += 1; stacks[c].append(pid); alive_len[pid] = c
    elif c < 0:
        p = -c; L = alive_len.pop(p); stacks[L].remove(p)
res = list(pre)
lens_seen = [c for c in pre if c > 0]
c_max = rng.choice([32, 36, 40])
ptr = {l: 0 for l in range(1, 41)}

def place(L):
    global pid
    pid += 1; stacks[L].append(pid); alive_len[pid] = L; res.append(L)

def rem(L, r):
    st = stacks[L]
    if not st: return False
    r = max(0, min(len(st) - 1, int(r)))
    p = st.pop(len(st) - 1 - r); del alive_len[p]; res.append(-p); return True

if mode != "cut":
    while len(res) < N:
        if mode == "phase2":
            place(rng.randint(8, c_max))
            while len(res) < N and rng.random() < 0.5:
                if rng.random() < 0.18:
                    L = rng.randint(1, 7)
                    if not stacks[L]: continue
                    step = 0
                    while rng.random() < 0.45: step += 1
                    r = min(len(stacks[L]) - 1, ptr[L] + step); ptr[L] = r
                    rem(L, r)
                else:
                    cand = [l for l in range(8, 41) if stacks[l]]
                    if not cand: continue
                    L = rng.choice(cand)
                    rem(L, rng.expovariate(1 / 1.5) if rng.random() < 0.84 else len(stacks[L]) * rng.random())
        else:
            L = rng.choice(lens_seen)
            place(L)
            while len(res) < N and rng.random() < 0.45:
                cand = [l for l in range(1, 41) if stacks[l]]
                if not cand: break
                L2 = L if stacks[L] and rng.random() < 0.7 else rng.choice(cand)
                rem(L2, rng.expovariate(1 / 2.0))
res = res[:N] if mode != "cut" else res
with open(out, "w") as f:
    f.write(f"{w} {len(res)}\n")
    f.write("\n".join(map(str, res)) + "\n0\n")
