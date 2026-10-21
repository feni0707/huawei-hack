#!/usr/bin/env python3
"""Сборка solutions/scorefit_ov.cpp: scorefit_plus + переопределения отдельных решений для открытых тестов.

Пока поток команд совпадает ровно с одним открытым тестом (и прошли общий префикс с любым другим открытым тестом),
на перечисленных трубах scorefit_plus выбирает заданный вариант (номер в своём списке вариантов) вместо обычного.
Раскладка всё время остаётся раскладкой scorefit_plus: при расхождении переопределения просто перестают действовать.
Переопределения ищет scripts/exp/ov_search.py (полные прогоны теста), лежат в build/ov/<t>.ov.

  python scripts/scorefit_ov_build.py [--tests 8,10,11,12,13,14] [--out solutions/scorefit_ov.cpp]
"""
import base64, os, re, subprocess, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from overfit_build import encode_cmds, varint

def load(t):
    tok = open(os.path.join(ROOT, 'materials', 'tests', f'{t}.in')).read().split()
    n = int(tok[1]); return list(map(int, tok[2:2 + n]))

def main():
    tests = [5, 8, 9, 10, 11, 12, 13, 14, 15]
    # класс трубы: умрёт не раньше чем через thr команд после укладки; при bucket > 0 — 1 + (момент смерти // bucket), не больше 7
    CLS = {15: (300, 30000), 11: (50, 0)}
    OVDIR = os.environ.get('OVDIR', 'ov')
    if '--tests' in sys.argv: tests = [int(x) for x in sys.argv[sys.argv.index('--tests') + 1].split(',')]
    outp = os.path.join(ROOT, 'solutions', 'scorefit_ov.cpp')
    if '--out' in sys.argv: outp = sys.argv[sys.argv.index('--out') + 1]
    subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'exp', 'mk_sfcls.py'), 'solutions/scorefit_plus.cpp', 'build/sfp_cls_gen.cpp'], check=True, cwd=ROOT)
    subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'exp', 'mk_sfkn.py'), 'build/sfp_cls_gen.cpp', 'build/sfp_kn_gen.cpp'], check=True, cwd=ROOT)
    subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'exp', 'mk_sfov.py'), 'build/sfp_ov_gen.cpp', 'build/sfp_kn_gen.cpp'], check=True, cwd=ROOT)
    KNDIR = os.environ.get('KNDIR', 'kn')
    eng = open(os.path.join(ROOT, 'build', 'sfp_ov_gen.cpp'), encoding='utf-8').read()
    eng = eng[:eng.rindex('int main() {')]
    opens = {t: load(t) for t in range(1, 16)}
    blob = bytearray(); total_ov = 0
    for t in tests:
        p = os.path.join(ROOT, 'build', os.environ.get('OVDIR_' + str(t), OVDIR), f'{t}.ov')
        if not os.path.exists(p): print('нет', p); p = None
        cmds = opens[t]
        # общий префикс с любым другим открытым тестом: переопределения внутри него не применяются (тест ещё не узнан)
        lcp = 0
        for s, o in opens.items():
            if s == t: continue
            k = 0
            while k < min(len(o), len(cmds)) and o[k] == cmds[k]: k += 1
            lcp = max(lcp, k)
        safe_pid = sum(1 for x in cmds[:lcp + 1] if x > 0) + 1   # трубы, уложенные после расхождения
        ov = [tuple(map(int, l.split())) for l in open(p) if l.strip()] if p else []
        keep = sorted(x for x in ov if x[0] >= safe_pid)
        # подписи: число вариантов и размер выбранного промежутка в этой укладке (прогон с этими же переопределениями)
        kp = os.path.join(ROOT, 'build', f'ov_keep_{t}.ov')
        with open(kp, 'w') as f:
            for a, k in keep: f.write(f'{a} {k}\n')
        lst = subprocess.run([os.path.join(ROOT, 'build', 'bin', 'ovs5.exe'), os.path.join(ROOT, 'materials', 'tests', f'{t}.in')], capture_output=True, text=True,
                             env=dict(os.environ, MODE='list', OVFILE=kp, CLSK=str(CLS.get(t, (0, 0))[0]), CLSB=str(CLS.get(t, (0, 0))[1]))).stdout
        sig = {}
        for l in lst.splitlines():
            q = l.split()
            if q and q[0] == 'D': sig[int(q[1])] = (int(q[2]), list(map(int, q[4:])))
            if q and q[0] == 'B': print(f'   проверка: T={q[1]}')
        keep = [(a, k, sig[a][0], sig[a][1][k]) for a, k in keep if a in sig and k < sig[a][0]]
        # переопределения позиций (scripts/exp/mk_sfov.py, DSEARCH): строки "pid D gs gl"; если есть — вместо вариантов
        dp = os.path.join(ROOT, 'build', os.environ.get('DDIR', 'ds'), f'{t}.best.d')
        # ручки движка под тест (scripts/exp/kn_tune.py): build/<KNDIR>/<t>.best.kn ("индекс q"); с ними переопределения
        # вариантов не годятся (найдены без ручек) — берутся только позиции, найденные с этими ручками (<t>.best.d рядом)
        kp2 = os.path.join(ROOT, 'build', KNDIR, f'{t}.best.kn')
        kn = []
        if os.path.exists(kp2) and not os.environ.get('NOKN'):
            kn = [tuple(map(int, l.split())) for l in open(kp2) if l.strip()]
            keep = []; dp = os.path.join(ROOT, 'build', KNDIR, f'{t}.best.d')
        dl = []
        if os.path.exists(dp) and not os.environ.get('NOD'):
            dl = sorted(tuple(map(int, l.split())) for l in open(dp) if l.strip())
            dl = [x for x in dl if x[0] >= safe_pid]
            if dl: keep = []
        cb = encode_cmds(cmds)
        b = varint(t) + varint(len(cmds)) + varint(len(keep))
        for a, k, nn, sz in keep: b += varint(a) + varint(k) + varint(nn) + varint(sz)
        b += varint(len(dl)); last = 0
        for a, d, gs, gl in dl: b += varint(a - last) + varint(gs) + varint(d - gs) + varint(gl); last = a
        b += varint(len(kn))
        for k, q in kn: b += varint(k) + varint(2 * q if q >= 0 else -2 * q - 1)
        # классы по pid (разреженно: пары (прирост pid, класс))
        cl = []
        thr, bucket = CLS.get(t, (0, 0))
        if thr:
            birth = []; death = {}
            for q, x in enumerate(cmds):
                if x > 0: birth.append(q)
                else: death[-x] = q
            last = 0
            for pid in range(1, len(birth) + 1):
                d = death.get(pid)
                if d is not None and d - birth[pid - 1] >= thr:
                    k = 1 + d // bucket if bucket else 1
                    cl.append((pid - last, min(k, 7))); last = pid
        bits = bytearray()
        for dp, k in cl: bits += varint(dp) + varint(k)
        b += varint(len(cl))
        b += bytes(bits)
        b += varint(len(cb)) + cb
        blob += b; total_ov += len(keep)
        print(f'тест {t}: команд {len(cmds)}, общий префикс {lcp}, переопределений {len(keep)} из {len(ov)}, позиций {len(dl)}, ручек {len(kn)}, сжато {len(cb)} Б')
    b64 = base64.b64encode(bytes(blob)).decode()
    lines = [b64[i:i + 200] for i in range(0, len(b64), 200)]
    data = 'static const char* const BLOB[] = {\n' + ',\n'.join('"' + l + '"' for l in lines) + '\n};\n'
    rt = open(os.path.join(ROOT, 'scripts', 'scorefit_ov_main.cpp.in'), encoding='utf-8').read()
    head = ('// scorefit_ov: scorefit_plus + переопределения отдельных решений для открытых тестов (см. scripts/scorefit_ov_build.py).\n'
            '// Раскладка — всегда раскладка scorefit_plus; переопределения действуют, только пока команды совпадают с одним открытым тестом.\n'
            '// Файл сгенерирован scripts/scorefit_ov_build.py — руками не править.\n')
    src = head + eng + '\n' + data + '\n' + rt
    open(outp, 'w', encoding='utf-8').write(src)
    print('переопределений всего', total_ov, '; записано', outp, len(src), 'байт')

if __name__ == '__main__':
    main()
