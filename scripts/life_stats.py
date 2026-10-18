import sys, collections, math
t = open(sys.argv[1]).read().split()
w, n = float(t[0]), int(t[1]); cmds = list(map(int, t[2:2+n]))
alive = {}  # pid -> (len, born_cmd_idx)
pid = 0
ages_pid = []; ages_cmd = []; dlen = collections.Counter(); plen = collections.Counter()
live_traj = []; cur = 0
runs = []; run = 0
for i, c in enumerate(cmds):
    if c > 0:
        pid += 1; alive[pid] = (c, i); plen[c] += 1; cur += c
        if run: runs.append(run); run = 0
    else:
        L, b = alive.pop(-c); ages_pid.append(pid - (-c)); ages_cmd.append(i - b); dlen[L] += 1; cur -= L; run += 1
    if i % (n // 20) == 0: live_traj.append((i, cur, len(alive)))
print("n", n, "places", pid, "dels", len(ages_pid))
print("live traj (cmd, sum, cnt):", live_traj)
def q(a, ps=(0.05, .25, .5, .75, .9, .99)):
    a = sorted(a); return [a[min(len(a) - 1, int(p * len(a)))] for p in ps]
print("age in pids quantiles .05 .25 .5 .75 .9 .99:", q(ages_pid))
print("age in cmds quantiles:", q(ages_cmd))
print("place len hist:", sorted(plen.items()))
print("del len hist:", sorted(dlen.items()))
hist = collections.Counter(runs); print("del-run lengths:", sorted(hist.items())[:20])
# survival vs length: fraction of placed pipes of length L ever deleted
print("deleted fraction by len:", {L: round(dlen[L] / plen[L], 2) for L in sorted(plen) if plen[L]})
