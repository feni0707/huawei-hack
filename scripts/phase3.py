import sys
t = int(sys.argv[1]) if len(sys.argv) > 1 else 10
start = int(sys.argv[2]) if len(sys.argv) > 2 else 7690
f = open(f'materials/tests/{t}.in').read().split(); n = int(f[1]); c = list(map(int, f[2:2 + n]))
ln = [0]; pid = 0; out = []; prev = None
for i, x in enumerate(c):
    if i < start:
        if x > 0:
            pid += 1; ln.append(x)
        continue
    if x > 0:
        pid += 1; ln.append(x); out.append(f"+{x}")
    else:
        p = -x
        out.append(f"-{ln[p]}(pid{p},d{(p - prev) if prev else 0},back{pid - p})"); prev = p
print(len(out), 'commands from', start)
print(' '.join(out[:70]))
print('...')
print(' '.join(out[-50:]))
