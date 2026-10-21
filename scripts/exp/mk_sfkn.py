#!/usr/bin/env python3
"""mk_sfkn.py src out : параметры движка scorefit_plus становятся изменяемыми («ручки»): веса LW[0..26] и параметры
участков и зон (27..42). knGet/knSet — снять/поставить вектор, KN_DEF — значения по умолчанию (захвачены при старте).
Пока ручки не меняют, ответы побитово как у исходника. Значения в файлах/посылке — целые q, ручка = q / 1e4."""
import sys
src, out = sys.argv[1], sys.argv[2]
s = open(src, encoding='utf-8').read()
def rep(old, new):
    global s
    assert s.count(old) == 1, (s.count(old), old[:70])
    s = s.replace(old, new)
rep('static const int CHUNK = 400, CHUNK_SAME = 800;', 'static int CHUNK = 400, CHUNK_SAME = 800;')
rep('static const ll MULTI_T = 50000;', 'static ll MULTI_T = 50000;')
for decl in ['static const int ZSLOTS = 16;', 'static const double ZSHARE = 0.08;', 'static const int ZMAXACT = 8;',
             'static const double ZMAXSHARE = 0.25;', 'static const int ZMINALIVE = 150;', 'static const double ZFAST = 0.7;',
             'static const double EXWARM = 300;', 'static const double ZGUARD_FREE = 0.15;', 'static const double ZLIFO = 0.7;']:
    rep(decl, decl.replace('static const ', 'static '))
i = s.index('static const double LW[27] = {')
j = s.index('\n', i)
line = s[i:j].replace('static const double LW[27]', 'static double LW[27]')
kn = r'''
// ручки (scripts/exp/mk_sfkn.py): 0..26 — LW, 27..42 — участки и зоны; KN_DEF — значения по умолчанию
static const int NKN = 43;
static void knGet(double* v) {
    for (int k = 0; k < 27; k++) v[k] = LW[k];
    v[27] = CHUNK; v[28] = CHUNK_SAME; v[29] = CHUNK_MIN; v[30] = CHUNK_PER_L; v[31] = CHUNK_LOW; v[32] = CHUNK_LOW_AGE;
    v[33] = ZSLOTS; v[34] = ZSHARE; v[35] = ZMAXACT; v[36] = ZMAXSHARE; v[37] = ZMINALIVE; v[38] = ZFAST; v[39] = EXWARM;
    v[40] = ZGUARD_FREE; v[41] = ZLIFO; v[42] = (double)MULTI_T;
}
static void knSet(const double* v) {
    for (int k = 0; k < 27; k++) LW[k] = v[k];
    CHUNK = (int)std::lround(v[27]); CHUNK_SAME = (int)std::lround(v[28]); CHUNK_MIN = (int)std::lround(v[29]);
    CHUNK_PER_L = (int)std::lround(v[30]); CHUNK_LOW = (int)std::lround(v[31]); CHUNK_LOW_AGE = (int)std::lround(v[32]);
    ZSLOTS = (int)std::lround(v[33]); ZSHARE = v[34]; ZMAXACT = (int)std::lround(v[35]); ZMAXSHARE = v[36];
    ZMINALIVE = (int)std::lround(v[37]); ZFAST = v[38]; EXWARM = v[39]; ZGUARD_FREE = v[40]; ZLIFO = v[41];
    MULTI_T = (ll)std::llround(v[42]);
}
static double KN_DEF[NKN];
static const bool knDefOk = (knGet(KN_DEF), true);
'''
s = s[:i] + line + kn + s[j + 1:]
open(out, 'w', encoding='utf-8').write(s)
print('written', out)
