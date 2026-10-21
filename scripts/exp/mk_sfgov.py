#!/usr/bin/env python3
"""mk_sfgov.py : solutions/scorefit.cpp -> scripts/exp/sfgov.cpp (scorefit с «губернатором»: набор теней с разными весами LW,
реальная раскладка ведётся весами лидера). Преобразование текстовое, чтобы легко повторять при обновлении scorefit.cpp."""
import re, sys
src = open('solutions/scorefit.cpp', encoding='utf-8').read()

# 1. веса: глобальные умолчания LW0, у экземпляра — свои веса LW
src = src.replace('static const double LW[27] =', 'static const double LW0[27] =', 1)
src = src.replace('struct Alloc {', 'struct SFAlloc {', 1)
src = src.replace('    ll T = 0;                 // длина траншеи', '    double LW[27];            // веса линейной оценки этого экземпляра (губернатор переключает набор)\n    ll T = 0;                 // длина траншеи', 1)
src = src.replace('    Alloc() {\n        if (logN[1] == 0) initLogs();', '    SFAlloc() {\n        for (int i = 0; i < 27; i++) LW[i] = LW0[i];\n        if (logN[1] == 0) initLogs();', 1)
src = src.replace('    ~Alloc() { delete shadow; }', '    ~SFAlloc() { delete shadow; }', 1)
src = src.replace('Alloc* shadow = 0;', 'SFAlloc* shadow = 0;', 1)
src = src.replace('shadow = new Alloc(*this);', 'shadow = new SFAlloc(*this);', 1)
src = src.replace('double Alloc::logN[Alloc::LOGN], Alloc::logSL[NB][MAXL + 1];', 'double SFAlloc::logN[SFAlloc::LOGN], SFAlloc::logSL[NB][MAXL + 1];', 1)
src = src.replace('    static double nbW(int r) {', '    double nbW(int r) const {', 1)

# 2. вырезаем ввод-вывод/main и добавляем губернатор
i = src.index('// ---------------------------------------------------------------- ввод-вывод')
head, io = src[:i], src[i:]
gov = r'''
// ---------------------------------------------------------------- губернатор
// Реальная раскладка ведётся экземпляром real с весами текущего лидера; тени (по одной на набор весов) проигрывают всю историю
// своей раскладкой. Лидер — набор с наименьшей оценкой «пик фронта + избыток дыр под фронтом» (затухающее среднее), с гистерезисом.
#ifndef GOV_SETS
#define GOV_SETS 3
#endif
#ifdef SFSRC
extern double K[64];   // ручки стенда: K[1] запас (доля T), K[2] вес дыр, K[3] dwell, K[4] затухание, K[5] период оценки
#endif
static const int GOV_MAXT = 1500000;   // при большей траншее тени не держим (память)
struct Alloc {
    SFAlloc real;
    SFAlloc* sh[GOV_SETS];
    double ws[GOV_SETS][27];
    double emaHoles[GOV_SETS], emaReal;
    ll emax[GOV_SETS], emaxReal;
    int cur = 0, since = 0, nplaced = 0;
    bool active = true;
    double kHoles = 1.0, decay = 0.995, margin = 0.004;
    int dwell = 500, every = 32;
    ll T = 0;                                // длина траншеи реальной раскладки (интерфейс стенда)
    const Paged& zl;                         // для стенда (диагностика)
    int tailStart() const { return real.tailStart(); }
    Alloc() : zl(real.zl) {
        for (int k = 0; k < GOV_SETS; k++) { for (int i = 0; i < 27; i++) ws[k][i] = LW0[i]; sh[k] = new SFAlloc(); emaHoles[k] = 0; emax[k] = 0; }
        // набор 1: без «остаток равен ходовой длине» (LW5, LW6) — устойчивее на классах, где спрос не совпадает с остатками
        ws[1][5] = 0; ws[1][6] = 0;
        // набор 2: чистый best fit по остатку (с потоками и зонами)
        for (int i = 1; i < 27; i++) ws[2][i] = 0;
        for (int k = 0; k < GOV_SETS; k++) for (int i = 0; i < 27; i++) sh[k]->LW[i] = ws[k][i];
        emaReal = 0; emaxReal = 0;
#ifdef SFSRC
        if (K[1] > 0) margin = K[1]; if (K[2] != 0) kHoles = K[2] < 0 ? 0 : K[2]; if (K[3] > 0) dwell = (int)K[3]; if (K[4] > 0) decay = K[4]; if (K[5] > 0) every = (int)K[5];
#endif
    }
    ~Alloc() { for (int k = 0; k < GOV_SETS; k++) delete sh[k]; }
    static ll holesUnder(const SFAlloc& a) { return (ll)a.tailStart() - a.occ; }
    void place(int L, ll& N, ll& D) {
        real.place(L, N, D);
        T = real.T;
        if (active) {
            ll n2, d2;
            for (int k = 0; k < GOV_SETS; k++) sh[k]->place(L, n2, d2);
            nplaced++;
            if (real.T > GOV_MAXT) { active = false; for (int k = 0; k < GOV_SETS; k++) { delete sh[k]; sh[k] = 0; } }
            else if (nplaced % every == 0) evaluate();
        }
    }
    ll remove(int pid) {
        if (active) for (int k = 0; k < GOV_SETS; k++) sh[k]->remove(pid);
        return real.remove(pid);
    }
    void evaluate() {
        since += every;
        for (int k = 0; k < GOV_SETS; k++) {
            double h = (double)holesUnder(*sh[k]);
            emaHoles[k] = emaHoles[k] * decay + h * (1 - decay);
            ll fe = (ll)sh[k]->tailStart(); if (fe > emax[k]) emax[k] = fe;
        }
        { emaReal = emaReal * decay + (double)holesUnder(real) * (1 - decay); ll fe = (ll)real.tailStart(); if (fe > emaxReal) emaxReal = fe; }
        if (since < dwell) return;
        // оценка набора k: пик фронта + kHoles·избыток дыр; текущий набор оцениваем по его тени (реальная раскладка — своя история)
        int best = cur; double bs = emax[cur] + kHoles * emaHoles[cur];
        for (int k = 0; k < GOV_SETS; k++) {
            double s = emax[k] + kHoles * emaHoles[k];
            if (s < bs) { bs = s; best = k; }
        }
        double sc = emax[cur] + kHoles * emaHoles[cur];
        if (best != cur && bs < sc - margin * (double)real.T - 20) {
            cur = best; since = 0;
            for (int i = 0; i < 27; i++) real.LW[i] = ws[cur][i];
        }
    }
};

'''
newmain = r'''int main() {
    Alloc* A = new Alloc();
    int c; char ob[64];
    while (readInt(c)) {
        if (c == 0) return 0;   // конец: ничего не печатаем
        int n;
        if (c > 0) {
            ll N, D;
            A->place(c, N, D);
            n = putNum(ob, N); ob[n++] = ' '; n += putNum(ob + n, D);
        } else {
            n = putNum(ob, A->remove(-c));
        }
        ob[n++] = '\n';
        writeOut(ob, n);
    }
    return 0;
}
'''
# в io оставляем всё до 'int main()' (чтение/запись), main заменяем
j = io.index('int main() {')
io_helpers = io[:j]
out = head + io_helpers + gov + newmain
open('scripts/exp/sfgov.cpp', 'w', encoding='utf-8').write(out)
print('written scripts/exp/sfgov.cpp', len(out))
