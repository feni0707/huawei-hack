import sys, os, math
d = sys.argv[1]; tot=0; totw=0
for i in range(1, 16):
    t = open(os.path.join(d, f"{i}.in")).read().split()
    w, n = float(t[0]), int(t[1]); cmds = list(map(int, t[2:2+n]))
    alive = {}; pid = 0; cur = 0; peak = 0
    for c in cmds:
        if c > 0: pid += 1; alive[pid] = c; cur += c; peak = max(peak, cur)
        else: cur -= alive.pop(-c)
    b = cur / (20*math.ceil(peak/20)); tot += b*w; totw += w
    print(f"test {i:2d} w={w:<4} UB_KPI={b:.4f}  weighted={b*w:.3f}")
print(f"sum UB={tot:.3f} of max {totw}")
