#!/usr/bin/env python3
"""cf_train.py data.txt [--epochs E] [--l2 L] [--tau T] [--target emax|T] : попарное обучение поправки c' = c + w·φ.
Печатает статистику решений, точность выбора на отложенной части (по файлам), вектор w (строка CFW)."""
import sys, math, random, collections
args = sys.argv[1:]
def opt(name, d):
    return type(d)(args[args.index(name) + 1]) if name in args else d
EP = opt('--epochs', 6); L2 = opt('--l2', 1e-3); TAU = opt('--tau', 20.0); TGT = opt('--target', 'emax'); LR = opt('--lr', 0.5)
FILES_HOLD = opt('--hold', 0.25)
dec = collections.OrderedDict()
NF = None
for line in open(args[0]):
    t = line.split()
    if not t or t[0] != 'R': continue
    key = (t[1], int(t[2]))
    em, T = int(t[6]), int(t[7]); f = list(map(float, t[8:]))
    NF = len(f)
    dec.setdefault(key, []).append((int(t[3]), int(t[5]), em, T, f))
# хвост без места (нужна докопка) — не настоящий выбор, если есть промежутки (так решает и scorefit)
for k in list(dec):
    os_ = dec[k]
    keep = [o for o in os_ if not (o[4][28] == 1 and o[4][29] == 0)]
    if len(keep) >= 2: dec[k] = keep
    elif len(keep) < len(os_): del dec[k]
files = sorted(set(k[0] for k in dec))
rnd = random.Random(1); rnd.shuffle(files)
hold = set(files[:int(len(files) * FILES_HOLD)])
def y(o): return o[2] if TGT == 'emax' else o[3] * 1000 + o[2]
tot = 0; withsig = 0; defbad = 0; gain = 0
for k, os_ in dec.items():
    tot += 1
    ys = [y(o) for o in os_]
    if max(ys) == min(ys): continue
    withsig += 1
    d = [o for o in os_ if o[1]]
    if d and y(d[0]) > min(ys): defbad += 1; gain += y(d[0]) - min(ys)
print(f'decisions {tot}, with signal {withsig}, default not best {defbad} ({100*defbad/max(withsig,1):.1f}% of signal), mean gain if oracle {gain/max(defbad,1):.1f} cells')
# пары
def pairs(keys):
    P = []
    for k in keys:
        os_ = dec[k]
        for i in range(len(os_)):
            for j in range(len(os_)):
                a, b = os_[i], os_[j]
                if y(a) < y(b): P.append((a[4], b[4], min(1.0, (y(b) - y(a)) / 40.0)))
    return P
tr_keys = [k for k in dec if k[0] not in hold]; ho_keys = [k for k in dec if k[0] in hold]
P = pairs(tr_keys)
print('train pairs', len(P), 'hold decisions', len(ho_keys))
w = [0.0] * NF
def cost(f, w):   # в исходных единицах: c = 100*f0
    return 100.0 * f[0] + sum(wi * fi for wi, fi in zip(w, f))
def acc(keys, w):
    good = 0; n = 0; reg = 0; dreg = 0
    for k in keys:
        os_ = dec[k]; ys = [y(o) for o in os_]
        if max(ys) == min(ys): continue
        n += 1
        ch = min(os_, key=lambda o: cost(o[4], w))
        reg += y(ch) - min(ys)
        d = [o for o in os_ if o[1]]
        if d: dreg += y(d[0]) - min(ys)
        good += y(ch) == min(ys)
    return good / max(n, 1), reg / max(n, 1), dreg / max(n, 1)
print('hold base: acc %.4f regret %.2f (default-chosen regret %.2f)' % acc(ho_keys, w))
for ep in range(EP):
    rnd.shuffle(P)
    lr = LR / (1 + ep)
    for fa, fb, wt in P:
        d = (cost(fa, w) - cost(fb, w)) / TAU   # хотим d < 0
        if d > 30: g = 1.0
        elif d < -30: g = 0.0
        else: g = 1.0 / (1.0 + math.exp(-d))
        g *= wt / TAU
        for i in range(NF):
            gi = g * (fa[i] - fb[i]) + L2 * w[i]
            if gi: w[i] -= lr * gi
    print('epoch', ep, 'hold acc %.4f regret %.2f' % acc(ho_keys, w)[:2], flush=True)
print('CFW=' + ','.join(f'{x:.4g}' for x in w))
