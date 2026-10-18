import sys, os, collections
d = sys.argv[1]
for i in range(1, 16):
    t = open(os.path.join(d, f"{i}.in")).read().split()
    w, n = float(t[0]), int(t[1])
    cmds = list(map(int, t[2:2+n]))
    tail = t[2+n:]
    lays = [c for c in cmds if c > 0]; rems = [c for c in cmds if c < 0]
    lens = collections.Counter(lays)
    # simulate alive set, peak alive length, removal age
    alive = {}; pid = 0; cur = 0; peak = 0; ages = []; maxpid_removed_ok = True
    removed = set()
    for idx, c in enumerate(cmds):
        if c > 0:
            pid += 1; alive[pid] = (c, idx); cur += c; peak = max(peak, cur)
        else:
            k = -c
            if k not in alive: maxpid_removed_ok = False; continue
            L, at = alive.pop(k); cur -= L; ages.append(pid - k)
    # removal bursts: count of consecutive removal runs
    runs = []; r = 0
    for c in cmds:
        if c < 0: r += 1
        else:
            if r: runs.append(r)
            r = 0
    if r: runs.append(r)
    ages.sort()
    med = ages[len(ages)//2] if ages else None
    print(f"test {i:2d}: w={w:<4} n={n:6d} lays={len(lays):6d} rems={len(rems):6d} "
          f"len[min={min(lays)},max={max(lays)},mean={sum(lays)/len(lays):.2f}] distinct={len(lens)} "
          f"sumLaid={sum(lays)} finalAlive={cur} peakAlive={peak} "
          f"remAge[med={med},max={ages[-1] if ages else None}] remRuns[cnt={len(runs)},max={max(runs) if runs else 0}] "
          f"validRemovals={maxpid_removed_ok} tail={tail}")
    top = lens.most_common(8)
    print("          top lengths:", top)
