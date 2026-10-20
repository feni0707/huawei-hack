#!/usr/bin/env python3
"""Случайные сценарии из блоков (длины × модель смерти × фазы) — поиск «чужих» классов, где настроенный scorefit проваливается.
  gen_scen.py --seed S --count K --out DIR [--sizes 2000,10000,40000]  -> DIR/scen/<k>.in  (+ DIR/scen/<k>.json с параметрами)"""
import argparse, json, os, random

def mk_len_sampler(rng):
    kind = rng.choice(["uni", "uni", "set", "bimodal", "geo", "pow"])
    if kind == "uni":
        a = rng.randint(1, 30); b = min(40, a + rng.randint(0, 40)); return ("uni", a, b), lambda r: r.randint(a, b)
    if kind == "set":
        k = rng.randint(2, 6); vals = [rng.randint(1, 40) for _ in range(k)]; wts = [rng.random() + 0.1 for _ in range(k)]
        return ("set", vals, wts), lambda r: r.choices(vals, wts)[0]
    if kind == "bimodal":
        a = rng.randint(1, 6); b = rng.randint(25, 40); pa = rng.uniform(0.3, 0.8)
        return ("bimodal", a, b, pa), lambda r: r.randint(1, a) if r.random() < pa else r.randint(b - 8 if b > 12 else 1, b)
    if kind == "geo":
        p = rng.uniform(0.05, 0.4); return ("geo", p), lambda r: min(40, 1 + int(r.expovariate(p)))
    e = rng.uniform(0.5, 3.0); return ("pow", e), lambda r: max(1, min(40, int(40 * r.random() ** e) + 1))

def mk_death(rng):
    kind = rng.choice(["lifo", "lifo", "fifo", "rand", "sizeexp", "bylen_lifo", "window"])
    if kind == "lifo":  # со стека: ранг от верха геометрический
        q = rng.uniform(0.2, 0.8); return ("lifo", q)
    if kind == "fifo":
        W = rng.randint(20, 2000); return ("fifo", W)
    if kind == "rand":
        return ("rand",)
    if kind == "sizeexp":  # срок жизни зависит от длины: mean = m0 * (L/20)^alpha
        m0 = rng.uniform(30, 3000); alpha = rng.uniform(-1.5, 1.5); return ("sizeexp", m0, alpha)
    if kind == "bylen_lifo":  # выбираем длину по текущим долям, внутри длины — вершина стека
        q = rng.uniform(0.3, 0.9); return ("bylen_lifo", q)
    w0 = rng.uniform(0.0, 0.7); wl = rng.uniform(0.05, 0.3); return ("window", w0, wl)

def run(rng, n, phases):
    cmds = []; alive = []; ln = {}; pid = 0; pos_in = {}
    bylen = {}
    due = {}
    def place(L):
        nonlocal pid
        pid += 1; cmds.append(L); ln[pid] = L; pos_in[pid] = len(alive); alive.append(pid); bylen.setdefault(L, []).append(pid); return pid
    def remove(p):
        i = pos_in.pop(p); last = alive.pop()
        if last != p: alive[i] = last; pos_in[last] = i
        bylen[ln[p]].remove(p); cmds.append(-p)
    ph_end = [int(n * sum(ph["w"] for ph in phases[:i + 1]) / sum(ph["w"] for ph in phases)) for i in range(len(phases))]
    pi = 0
    while len(cmds) < n - 1:
        while pi < len(phases) - 1 and len(cmds) >= ph_end[pi]: pi += 1
        ph = phases[pi]; dm = ph["death"]; samp = ph["samp"]
        # due-deletions (sizeexp)
        if dm[0] == "sizeexp":
            idxs = [i for i in list(due) if i <= len(cmds)]
            for i in idxs:
                for q in due.pop(i):
                    if q in pos_in and len(cmds) < n - 1: remove(q)
        if rng.random() < ph["pplace"] or not alive:
            L = samp(rng); p = place(L)
            if dm[0] == "sizeexp":
                mean = max(1.0, dm[1] * (L / 20.0) ** dm[2]); due.setdefault(len(cmds) + int(rng.expovariate(1.0 / mean)) + 1, []).append(p)
        else:
            if dm[0] == "sizeexp": continue
            if dm[0] == "lifo":
                k = 0
                while k < len(alive) - 1 and rng.random() > dm[1]: k += 1
                st = sorted(alive)[-1 - k] if len(alive) < 400 else alive[-1 - k]  # приближённо: по порядку PID
                # ключевое: хотим «свежие» — возьмём k-ю с конца по PID
                srt = max(alive) if k == 0 else sorted(alive, reverse=True)[:k + 1][-1]
                remove(srt)
            elif dm[0] == "fifo":
                if len(alive) > dm[1]: remove(min(alive))
                else: continue
            elif dm[0] == "rand":
                remove(rng.choice(alive))
            elif dm[0] == "bylen_lifo":
                Ls = [l for l, v in bylen.items() if v]
                L = rng.choice(Ls); st = bylen[L]; k = 0
                while k < len(st) - 1 and rng.random() > dm[1]: k += 1
                remove(st[-1 - k])
            elif dm[0] == "window":
                srt = sorted(alive); a = int(len(srt) * dm[1]); b = max(a + 1, int(len(srt) * (dm[1] + dm[2])))
                remove(srt[rng.randrange(a, min(b, len(srt)))] if b > a else srt[0])
    cmds.append(0) if False else None
    return cmds

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--seed", type=int, default=1); ap.add_argument("--count", type=int, default=20)
    ap.add_argument("--out", default="build/scen"); ap.add_argument("--sizes", default="2000,10000,40000")
    a = ap.parse_args(); sizes = [int(x) for x in a.sizes.split(",")]
    d = os.path.join(a.out, "scen"); os.makedirs(d, exist_ok=True)
    for k in range(a.count):
        rng = random.Random(f"{a.seed}:{k}")
        n = rng.choice(sizes)
        nph = rng.choice([1, 1, 2, 2, 3])
        phases = []
        for _ in range(nph):
            desc, samp = mk_len_sampler(rng)
            phases.append(dict(w=rng.uniform(0.5, 2), pplace=rng.uniform(0.45, 0.9), death=mk_death(rng), samp=samp, ldesc=desc))
        cmds = run(rng, n, phases)
        with open(os.path.join(d, f"{k}.in"), "w") as f:
            f.write(f"1 {len(cmds)}\n" + "\n".join(map(str, cmds)) + "\n0\n")
        json.dump([dict(w=p["w"], pplace=p["pplace"], death=p["death"], len=p["ldesc"]) for p in phases], open(os.path.join(d, f"{k}.json"), "w"))
main()
