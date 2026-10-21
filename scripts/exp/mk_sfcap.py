#!/usr/bin/env python3
"""mk_sfcap.py [out] : scorefit.cpp -> копия с ограничением «предварительных» докопок (участки потоков и зоны по длинам).
Ручки (env в стенде, иначе значения по умолчанию): CAPF — участок/зона не больше max(L, CAPF·T) (0 — выкл.);
CAPI — если свободно > CAPI·T (траншея раздута), участок/зона = только нужное (0 — выкл.)."""
import sys
out = sys.argv[1] if len(sys.argv) > 1 else 'scripts/exp/sf_cap.cpp'
s = open(sys.argv[2] if len(sys.argv) > 2 else 'solutions/scorefit.cpp', encoding='utf-8').read()
knobs = r'''
#ifdef SFSRC
#include <cstdlib>
static double capF() { static double v = getenv("CAPF") ? atof(getenv("CAPF")) : CAPF_DEF; return v; }
static double capI() { static double v = getenv("CAPI") ? atof(getenv("CAPI")) : CAPI_DEF; return v; }
static bool capV() { static bool v = getenv("CAPV") != 0; return v; }
static int capR() { static int v = getenv("CAPR") ? atoi(getenv("CAPR")) : CAPR_DEF; return v; }
#else
static double capF() { return CAPF_DEF; }
static double capI() { return CAPI_DEF; }
static bool capV() { return CAPV_DEF; }
static int capR() { return CAPR_DEF; }
#endif
'''
defs = '#ifndef CAPR_DEF\n#define CAPR_DEF 0\n#endif\n#ifndef CAPV_DEF\n#define CAPV_DEF 0\n#endif\n#ifndef CAPF_DEF\n#define CAPF_DEF 0.0\n#endif\n#ifndef CAPI_DEF\n#define CAPI_DEF 0.0\n#endif\n'
i = s.index('typedef long long ll;')
s = s[:i] + defs + s[i:]
i = s.index('static const double LW[27]')
s = s[:i] + knobs + s[i:]
# метод: допустимый размер предварительной докопки (клетки) для трубы L
s = s.replace('    // ---- тень без потоков', r'''    ll capRes(ll want, int L) const {   // ограничение предварительной резервации (участок потока / зона)
        ll lim = want;
        if (capV() && verB) return want;
        if (capR() > 0 && runIdx >= capR()) return want;   // идёт длинная серия равных длин — участок нужен целиком   // потоки подтверждены (B вымирает) — участки нужны целиком
        if (capF() > 0) { ll c = (ll)(capF() * (double)T); if (c < L) c = L; if (c < lim) lim = c; }
        if (capI() > 0 && (double)(T - occ) > capI() * (double)T && L < lim) lim = L;
        return lim;
    }
    // ---- тень без потоков''', 1)
old = '                if (far && chunk - fl > need) need = chunk - fl;'
assert old in s
s = s.replace(old, '                if (far && chunk - fl > need) need = chunk - fl;\n                if (need > L - fl) { ll cr = capRes(need + fl, L) - fl; if (cr < L - fl) cr = L - fl; need = (int)cr; }', 1)
old = '                if (st >= 0 && chunk - tl > need) need = chunk - tl;'
assert old in s
s = s.replace(old, '                if (st >= 0 && chunk - tl > need) need = chunk - tl;\n                if (st >= 0 && need > L - tl) { ll cr = capRes(need + tl, L) - tl; if (cr < L - tl) cr = L - tl; need = (int)cr; }', 1)
old = '            ll C = (ll)ZSLOTS * L;'
assert old in s
s = s.replace(old, '            ll C = (ll)ZSLOTS * L;\n            { ll cr = capRes(C, L); cr = cr / L * L; if (cr < L) cr = L; C = cr; }', 1)

s = s.replace('static const int CHUNK_MIN = 200, CHUNK_PER_L = 10;', 'static int CHUNK_MIN = 200, CHUNK_PER_L = 10;', 1)
s = s.replace('static const int CHUNK_LOW = 40, CHUNK_LOW_AGE = 1000;', 'static int CHUNK_LOW = 40, CHUNK_LOW_AGE = 1000;', 1)
s = s.replace('static double capF() { static double v', 'static void chunkKnobs() { static bool d = false; if (d) return; d = true; if (getenv("CMIN")) CHUNK_MIN = atoi(getenv("CMIN")); if (getenv("CPL")) CHUNK_PER_L = atoi(getenv("CPL")); if (getenv("CLAGE")) CHUNK_LOW_AGE = atoi(getenv("CLAGE")); }\nstatic double capF() { chunkKnobs(); static double v', 1)

s = s.replace('    Alloc() {\n        if (logN[1] == 0) initLogs();', '    Alloc() {\n#ifdef SFSRC\n        chunkKnobs();\n#endif\n        if (logN[1] == 0) initLogs();', 1)

old_ds = "        double ds = (s <= MAXL && cntTot > 0) ? cntR[s] / cntTot : 0;"
assert old_ds in s
s = s.replace(old_ds, """        double ds = (s <= MAXL && cntTot > 0) ? cntR[s] / cntTot : 0;
        if (dsW() > 0 && s <= MAXL + dsW() && cntTot > 0) {   // спрос на дыру s — доля длин (s - w, s] (класс крупных труб, а не одна длина)
            double a = 0; for (int l = s - dsW() + 1; l <= s; l++) if (l >= 1 && l <= MAXL) a += cntR[l];
            ds = a / cntTot * dsK();
        }""", 1)
s = s.replace('static int capR() { static int v', 'static int dsW() { static int v = getenv("DSW") ? atoi(getenv("DSW")) : DSW_DEF; return v; }\nstatic double dsK() { static double v = getenv("DSK") ? atof(getenv("DSK")) : 1.0; return v; }\nstatic int capR() { static int v', 1)
s = s.replace('static int capR() { return CAPR_DEF; }', 'static int capR() { return CAPR_DEF; }\nstatic int dsW() { return DSW_DEF; }\nstatic double dsK() { return DSK_DEF; }', 1)
s = s.replace('#ifndef CAPR_DEF', '#ifndef DSW_DEF\n#define DSW_DEF 0\n#endif\n#ifndef DSK_DEF\n#define DSK_DEF 1.0\n#endif\n#ifndef CAPR_DEF', 1)
open(out, 'w', encoding='utf-8').write(s)
print('written', out)
