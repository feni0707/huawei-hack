#!/usr/bin/env python3
"""mk_sfcf.py [src] [out] : копия scorefit с (1) записью кандидатов и их признаков, (2) принудительным выбором варианта,
(3) обученной поправкой к стоимости (env CFW="w0,w1,..." или -DCFW_BAKED), (4) драйвером контрфактических прогонов (-DCFDUMP).
Варианты решения вне потоков/зон: кандидаты-промежутки (как в scorefit) + хвост. Поправка: c' = c + Σ w_i·φ_i."""
import sys
src = sys.argv[1] if len(sys.argv) > 1 else 'solutions/scorefit.cpp'
out = sys.argv[2] if len(sys.argv) > 2 else 'scripts/exp/sf_cf.cpp'
s = open(src, encoding='utf-8').read()

NF = 32
glob = r'''
// ---------------------------------------------------------------- контрфактические варианты (Claude)
static const int NFT = %d;
struct CfOpt { int g, s; double c; double f[NFT]; };
static int gForcePid = -1, gForceIdx = -1, gChosen = -1;
static bool gRecord = false;
static double gJitter = 0; static unsigned long long gJs = 1;
static double jrnd() { gJs ^= gJs << 13; gJs ^= gJs >> 7; gJs ^= gJs << 17; return (double)(gJs %% 1000003) / 1000003.0; }
static std::vector<CfOpt> gOpts;             // варианты последней записанной укладки (промежутки, затем хвост)
static double gCW[NFT]; static bool gCWset = false;
#ifdef CFW_BAKED
static const double gCWb[NFT] = { CFW_BAKED };
#endif
static void cfInit() {
    static bool done = false; if (done) return; done = true;
    for (int i = 0; i < NFT; i++) gCW[i] = 0;
#ifdef CFW_BAKED
    for (int i = 0; i < NFT; i++) gCW[i] = gCWb[i]; gCWset = true;
#endif
#ifdef SFSRC
    if (const char* e = getenv("CFW")) { char* q = (char*)e; for (int i = 0; i < NFT && *q; i++) { char* r; double v = strtod(q, &r); if (r == q) break; gCW[i] = v; q = r; while (*q == ',' || *q == ' ') q++; } gCWset = true; }
#endif
}
''' % NF
i = s.index('typedef long long ll;')
s = s[:i] + '#include <cstdlib>\n#include <cmath>\n#include <vector>\n' + s[i:]
i = s.index('static const double LW[27]')
s = s[:i] + glob + s[i:]

# признаки кандидата (вне потоков) — функция рядом с gapCost
feat = r'''
    // признаки кандидата для обученной поправки (вне потоков; см. CfOpt)
    void gapFeat(int gs, int s, int L, int m, double selfHz, double base, double* f) const {
        for (int i = 0; i < NFT; i++) f[i] = 0;
        int r = s - L;
        int lp = pidLeftOf(gs), rp = pidRightOf(gs + s);
        int cur = (int)pos.size();
        double hl = hazard(lp), hr = hazard(rp);
        bool dead = r > 0 && r < m;
        double dr = (r >= 1 && r <= MAXL && cntTot > 0) ? cntR[r] / cntTot : 0;
        double ds = (s <= MAXL && cntTot > 0) ? cntR[s] / cntTot : 0;
        double ns = s < NB ? gapCnt[s] : 0, nr = r < NB ? gapCnt[r] : 0;
        f[0] = base / 100.0;
        f[1] = r / 40.0; f[2] = r == 0; f[3] = dead; f[4] = dead ? r / 40.0 : 0; f[5] = r >= m;
        f[6] = dr; f[7] = ds; f[8] = hl; f[9] = hr;
        f[10] = log1n((int)ns) / 5; f[11] = (s < NB ? logSL[s][L] : std::log((double)s / L));
        f[12] = (r > 0 && r < 2 * m && !dead);
        f[13] = hzRate(lp); f[14] = hzRate(rp);
        f[15] = lp > 0 ? std::log(1.0 + cur - lp) / 12 : 0; f[16] = rp > 0 ? std::log(1.0 + cur - rp) / 12 : 0;
        f[17] = lp > 0 && len[lp] == L; f[18] = rp > 0 && len[rp] == L;
        f[19] = (double)gs / (double)(T > 0 ? T : 1);
        f[20] = log1n((int)nr) / 5;
        f[21] = cntTot > 0 ? cntR[L] / cntTot : 0;
        f[22] = selfHz; f[23] = L / 40.0; f[24] = s / 40.0;
        f[25] = lp > 0 ? dabs(std::log(1.0 + cur - lp) - std::log(1.0 + cur - (rp > 0 ? rp : lp))) / 12 : 0;
        f[26] = (lp > 0 && rp > 0) ? 1 : 0;
        f[27] = (double)(T - occ) / (double)(T > 0 ? T : 1);
        // 28..31 — хвост (у промежутков 0)
    }
    double candCost(int gs, int s, int L, int st, int m, double selfHz) {
        double c = gapCost(gs, s, L, st, m, selfHz);
        if (gJitter > 0 && c < 1e17) c += gJitter * jrnd();
        if (st >= 0 || T < 60 || c >= 1e17) return c;
        if (!gCWset && !gRecord && gForcePid < 0) return c;
        CfOpt o; o.g = gs; o.s = s;
        gapFeat(gs, s, L, m, selfHz, c, o.f);
        if (gCWset) { double d = 0; for (int i = 0; i < NFT; i++) d += gCW[i] * o.f[i]; c += d; }
        o.c = c;
        if (gRecord || gForcePid >= 0) gOpts.push_back(o);
        return c;
    }
'''
i = s.index('    // ---- зоны по длинам.')
s = s[:i] + feat + s[i:]

# в place: кандидаты через candCost, принудительный выбор, хвост как вариант
old = '        double bestC = 1e18; int bgs = -1, bs = 0;'
assert old in s
s = s.replace(old, '        double bestC = 1e18; int bgs = -1, bs = 0;\n        cfInit(); gOpts.clear();', 1)
s = s.replace('zoneBlocked(g, s, st) ? 1e17 : gapCost(g, s, L, st, m, selfHz);', 'zoneBlocked(g, s, st) ? 1e17 : candCost(g, s, L, st, m, selfHz);')
s = s.replace('zoneBlocked(gk[i], s, st) ? 1e17 : gapCost(gk[i], s, L, st, m, selfHz);', 'zoneBlocked(gk[i], s, st) ? 1e17 : candCost(gk[i], s, L, st, m, selfHz);')
s = s.replace('zoneBlocked(g, bl, st) ? 1e17 : gapCost(g, bl, L, st, m, selfHz);', 'zoneBlocked(g, bl, st) ? 1e17 : candCost(g, bl, L, st, m, selfHz);')
old = '        bool useGap = bgs >= 0 && bestC < 1e17;'
assert old in s
new = r'''        bool useGap = bgs >= 0 && bestC < 1e17;
        bool tailRoom = T - tailStart() >= L;
        if (st < 0 && T >= 60 && (gRecord || gForcePid >= 0 || gCWset)) {
            // хвост как вариант: базовая стоимость LW[13] (если есть место без докопки), иначе 300
            CfOpt o; o.g = -1; o.s = 0; for (int i = 0; i < NFT; i++) o.f[i] = 0;
            double base = tailRoom ? LW[13] : 300; double baseUse = tailRoom ? LW[13] : 1e9;
            o.f[0] = base / 100.0; o.f[28] = 1; o.f[29] = tailRoom; o.f[30] = tailRoom ? 0 : (double)((L - (T - tailStart()) + 19) / 20) / 5;
            o.f[31] = (double)(T - occ) / (double)(T > 0 ? T : 1);
            double c = baseUse; if (gCWset && tailRoom) for (int i = 0; i < NFT; i++) c += gCW[i] * o.f[i];
            o.c = c;
            if (gRecord || gForcePid >= 0) gOpts.push_back(o);
            if (gCWset) {   // обученный выбор: хвост против лучшего промежутка
                if (useGap && c < bestC) useGap = false;
                else if (!useGap && bgs >= 0 && bestC < 1e17 && bestC <= c) useGap = true;
            }
            if (gForcePid == pid && gForceIdx >= 0 && gForceIdx < (int)gOpts.size()) {
                const CfOpt& f = gOpts[gForceIdx];
                if (f.g < 0) useGap = false; else { bgs = f.g; bs = f.s; bestC = -1e18; useGap = true; }
            }
        }
        if (gCWset && st < 0 && T >= 60) goto cfdecided;'''
s = s.replace(old, new, 1)
old = '        if (useGap && st < 0 && T - tailStart() >= L && LW[13] < bestC) useGap = false;'
assert old in s
s = s.replace(old, old + '\n        cfdecided:\n        if (gRecord) { gChosen = -1; for (size_t k = 0; k < gOpts.size(); k++) if (useGap ? (gOpts[k].g == bgs && gOpts[k].s == bs) : gOpts[k].g < 0) { gChosen = (int)k; break; } }', 1)

# клон с глубокой копией тени
s = s.replace('    ~Alloc() { delete shadow; }', '    ~Alloc() { delete shadow; }\n    Alloc* clone() const { Alloc* c = new Alloc(*this); c->shadow = shadow ? new Alloc(*shadow) : 0; return c; }', 1)

drv = r'''
#ifdef CFDUMP
// драйвер: для выборки укладок — все варианты, досчёт до конца теста по умолчанию; строка на вариант:
// R <тест> <ci> <idx> <nopt> <isdef> <Emax> <T> f0..f(NFT-1)
#include <cstring>
static std::vector<int> readCmds(const char* fn) { FILE* f = fopen(fn, "r"); double w; int n; std::vector<int> c; if (!f || fscanf(f, "%lf %d", &w, &n) != 2) return c; c.resize(n); for (int i = 0; i < n; i++) if (fscanf(f, "%d", &c[i]) != 1) { c.resize(i); break; } fclose(f); return c; }
static void stepRun(Alloc* A, const std::vector<int>& c, size_t from, ll& emax) {
    for (size_t i = from; i < c.size(); i++) { int x = c[i]; if (x == 0) break; if (x > 0) { ll N, D; A->place(x, N, D); } else A->remove(-x); ll e = A->tailStart(); if (e > emax) emax = e; }
}
int main(int argc, char** argv) {
    double prob = getenv("SAMPLE") ? atof(getenv("SAMPLE")) : 0.2; int maxopt = getenv("MAXOPT") ? atoi(getenv("MAXOPT")) : 6;
    unsigned long long rs = 88172645463325252ULL ^ (getenv("SEED") ? (unsigned long long)atoll(getenv("SEED")) : 1);
    for (int a = 1; a < argc; a++) {
        std::vector<int> c = readCmds(argv[a]);
        Alloc* A = new Alloc(); ll emax = 0;
        for (size_t i = 0; i < c.size(); i++) {
            int x = c[i]; if (x == 0) break;
            if (x > 0) {
                rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17;
                if ((double)(rs % 1000000) / 1e6 < prob && A->T >= 60) {
                    Alloc* P = A->clone(); gRecord = true; ll N, D; P->place(x, N, D); gRecord = false;
                    std::vector<CfOpt> opts = gOpts; int defc = gChosen; delete P;
                    int def = defc;
                    if (opts.size() >= 2) {
                        // отбор: maxopt дешёвых по базовой стоимости + вариант по умолчанию + хвост
                        std::vector<int> idx(opts.size()); for (size_t k = 0; k < opts.size(); k++) idx[k] = (int)k;
                        std::sort(idx.begin(), idx.end(), [&](int p, int q) { return opts[p].f[0] < opts[q].f[0]; });
                        std::vector<char> use(opts.size(), 0);
                        for (int k = 0; k < (int)idx.size() && k < maxopt; k++) use[idx[k]] = 1;
                        if (def >= 0) use[def] = 1; use[opts.size() - 1] = opts.back().g < 0 ? 1 : use[opts.size() - 1];
                        for (size_t k = 0; k < opts.size(); k++) {
                            if (!use[k]) continue;
                            Alloc* C = A->clone(); gForcePid = (int)C->pos.size(); gForceIdx = (int)k;
                            ll N2, D2; C->place(x, N2, D2); gForcePid = -1; gForceIdx = -1;
                            ll em = emax; { ll e = C->tailStart(); if (e > em) em = e; }
                            int K = getenv("CFK") ? atoi(getenv("CFK")) : 1; double jit = getenv("CFJ") ? atof(getenv("CFJ")) : 0;
                            ll emSum = 0, tSum = 0;
                            for (int kk = 0; kk < K; kk++) {
                                Alloc* E = C->clone(); ll e2 = em; gJitter = jit; gJs = 0x9E3779B97F4A7C15ULL * (kk + 1) + i;
                                stepRun(E, c, i + 1, e2); gJitter = 0;
                                emSum += e2; tSum += E->T; delete E;
                            }
                            em = emSum / K; C->T = tSum / K;
                            printf("R %s %zu %zu %zu %d %lld %lld", argv[a], i, k, opts.size(), (int)k == def, em, C->T);
                            for (int j = 0; j < NFT; j++) printf(" %.5g", opts[k].f[j]);
                            printf("\n");
                            delete C;
                        }
                    }
                }
                ll N, D; A->place(x, N, D);
            } else A->remove(-x);
            ll e = A->tailStart(); if (e > emax) emax = e;
        }
        fflush(stdout);
        delete A;
    }
    return 0;
}
#define main sf_main_unused_cf
#endif
'''
i = s.index('int main() {')
s = s[:i] + drv + s[i:]
s = s.replace('#include <cstdio>\n', '#include <cstdio>\n#include <algorithm>\n', 1)
open(out, 'w', encoding='utf-8').write(s)
print('written', out)
