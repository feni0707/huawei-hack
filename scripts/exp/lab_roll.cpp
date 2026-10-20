// scorefit_lab: ЭКСПЕРИМЕНТАЛЬНАЯ версия scorefit со всеми ручками P[] и оффлайн-симулятором (-DSIM).
// Посылать solutions/scorefit.cpp (чистая версия, ответы побитово совпадают с этой при параметрах по умолчанию).
// Сборка оффлайн-симулятора (без интерактора, для экспериментов):
//   c++ -std=c++17 -O2 -DSIM -o build/sim_sf solutions/scorefit_lab.cpp
//   build/sim_sf materials/tests/*.in      (параметры: env SF="i=v i=v ...")
//   SFCUTS=20 — KPI каждого теста усредняется по 20 концам в последних 10 % (убирает шум «где кончился тест»
//   относительно открытого участка потоков: на маленьких тестах он ≈ 400/√12 м)
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <unistd.h>

typedef long long ll;

// ---------------------------------------------------------------- параметры
// P0 затухание статистики длин, P1 порог доли для minRecent, P2 цена мёртвого метра остатка,
// P3 цена живого метра остатка, P4 сколько промежутков каждой длины смотреть, P5 бонус точного попадания,
// P6 цена метра остатка (best fit), P7 порог стоимости, выше которого кладём в хвост (0 — никогда),
// P8 выбор стороны (0 лево, 1 остаток к молодому соседу, 2 к старому, 3 право)
// P10 потоки (0 нет, 1 оракул, 2 по признакам), P11 горизонт оракула, P20 размер участка для потока B, м,
// P21 мин. число пар удалений для включения, P22 порог каппы, P23 затухание рождений, P24 затухание пар
// P15 мёртвый остаток >= P15 — в хвост (0 выкл), P16 штраф мёртвого остатка с учётом риска смерти соседей
//      0  1      2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
// P13 «крошечная траншея» (T < P13): точное попадание, иначе самая большая дыра; P32 — штраф точного
//      попадания у соседа с высокой долей смертей его длины (только в этом режиме)
// P33=1 — правило P15 только вне потоков, P33=2 — промежутки с мёртвым остатком ≥ P15 исключаются из кандидатов (вне потоков); P36 — мин. доля редкого значения признака для каппы;
// P39 — штраф замуровывания, только если s + len(соседа) встречается среди недавних укладок с долей ≥ P39
// P53 эндшпиль (последние P53 команд до лимита P54 или 200000 — чистый best fit; вредно, выкл.);
// P55 зона B при подтверждённых потоках (порог P9): 1 — B только в дыры между B, 2 — ещё и A не в дыры между двумя B;
// P80 зоны по длинам (1 — вкл.), P81 слотов в участке, P82 мин. доля длины, P83 макс. число частых длин, P85 макс. доля
//      самой частой, P86 мин. живых длины для участка, P87/P88 детектор смены режима (доля длин с зонами при затухании P88),
//      P89 прогрев признака LW25 (cntTot); P90 участок внутри незарезервированной дыры до хвоста; P91 мин. доля удалений
//      «по стеку длины» (среди 8 последних той же длины) для зон; LW24 хвост с докопкой, LW25/LW26 изменение избытка дыр над спросом
// P56 сторона остатка по порядку смертей (умирают старые — к старшему соседу), P57 затухание этой статистики
// P92 «гонка соседей» (какой из соседей промежутка с остатком умер первым) вместо P56, если доля дальше P95 от 0.5;
//      P93 затухание; P92 = 2..5 — статистика по ключам (возраст/длины соседей, P94 мин. наблюдений), не лучше 1
// P107 участки потоков по очереди зеркальны ([A→ ←B][B→ ←A]…): области одного потока соседних участков смежны;
// P99 тень без потоков (копия при первой возможной метке): участок ≤ max(P100, P106·L, P101·экономия), экономия =
//      свободное внутри траншеи у тени минус у основной (без хвоста и граничного промежутка). Отрицательные/выкл.:
//      P96/P97 участок ~ L, P98 ближний поток не копает участок, P102 участок кратен L, P103–P105 «помеченные не
//      умирают» → без участка, P108 нижняя граница участка ~ T; P109 — тень удаляется при T > P109 (память ×2)
#ifdef SIM
#include <vector>
static const std::vector<int>* gCmds = nullptr;
#include <unordered_map>
#include <map>
struct TokPred {
    static const int MAXO = 8;
    std::unordered_map<unsigned long long, std::map<int, int>> ctx[MAXO + 1];
    std::vector<int> hist;
    double minShare = 0.75; int minCnt = 2;
    static unsigned long long key(const std::vector<int>& h, int o) { unsigned long long k = 1469598103934665603ULL; for (int i = (int)h.size() - o; i < (int)h.size(); i++) k = (k ^ (unsigned long long)(h[i] + 7)) * 1099511628211ULL; return k; }
    void push(int t) {
        for (int o = 1; o <= MAXO && o <= (int)hist.size(); o++) ctx[o][key(hist, o)][t]++;
        hist.push_back(t);
    }
    int predict(const std::vector<int>& h) const {   // -1 — нет уверенного прогноза
        for (int o = MAXO; o >= 1; o--) {
            if ((int)h.size() < o) continue;
            auto it = ctx[o].find(key(h, o)); if (it == ctx[o].end()) continue;
            int tot = 0, bt = -1, bc = 0; for (auto& kv : it->second) { tot += kv.second; if (kv.second > bc) { bc = kv.second; bt = kv.first; } }
            if (tot >= minCnt && bc >= minShare * tot) return bt;
        }
        return -1;
    }
};
static TokPred* gTP = nullptr;
static TokPred* gTPv[3] = {nullptr, nullptr, nullptr};
static double gAcc[3] = {0, 0, 0}, gAccDecay = 0.98;
static int gView = 0;
static double gHacc = 0.0, gHdecay = 0.95;      // доля совпавших шагов прогноза (затухание)
static std::vector<int> gPend; static size_t gPendIdx = 0; static int gPendView = 0, gSkip = 0;
static void checkPend(int viewTok[3]) {          // viewTok[v] — токен реальной команды в представлении v
    if (gSkip > 0) { gSkip--; return; }
    if (gPendIdx < gPend.size()) {
        bool ok = gPend[gPendIdx] == viewTok[gPendView];
        gHacc = gHacc * gHdecay + (ok ? 1 - gHdecay : 0);
        gPendIdx++;
        if (!ok) gPend.clear();                   // после промаха остаток прогноза не проверяем
    }
}
static int encDel(int v, int L, int rank, int delta, int back) {
    if (v == 0) return 100 + L * 16 + (rank > 15 ? 15 : rank);
    if (v == 1) return 100000 + (delta + 5000);        // delta в [-5000, +5000]
    return 200000 + (back > 99999 ? 99999 : back);
}
static void tpPushP(int L) {
    { int vt[3] = {L, L, L}; checkPend(vt); }
    for (int v = 0; v < 3; v++) { if (!gTPv[v]) continue; std::vector<int> h(gTPv[v]->hist.end() - std::min<size_t>(gTPv[v]->hist.size(), TokPred::MAXO), gTPv[v]->hist.end()); int pr = gTPv[v]->predict(h); gAcc[v] = gAcc[v] * gAccDecay + (pr == L ? 1 - gAccDecay : 0); gTPv[v]->push(L); }
}
static void tpPushD(int L, int rank, int delta, int back) {
    { int vt[3] = {encDel(0, L, rank, delta, back), encDel(1, L, rank, delta, back), encDel(2, L, rank, delta, back)}; checkPend(vt); }
    for (int v = 0; v < 3; v++) { if (!gTPv[v]) continue; int tk = encDel(v, L, rank, delta, back); std::vector<int> h(gTPv[v]->hist.end() - std::min<size_t>(gTPv[v]->hist.size(), TokPred::MAXO), gTPv[v]->hist.end()); int pr = gTPv[v]->predict(h); gAcc[v] = gAcc[v] * gAccDecay + (pr == tk ? 1 - gAccDecay : 0); gTPv[v]->push(tk); }
}
#endif
static double P[128] = {0.999, 0.002, 0, 0, 1, 0, 1, 0, 1, 0.9, 2, 0, 3, 60, 6, 0, 0, 3, 4, 0, 400, 20, 0.6, 0.999, 0.998, 0, 0, 0.99, 0.5, 0, 4, 0.98, 2000, 2, 0, 0, 0.3, 0, 0.99, 0.05, 0, 0, 0, 0, 1, 0, 0, 0, 8, 0.999, 1.5, 0, 0, 0, 0, 2, 1, 0.99, 0, 1, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 16, 0.08, 8, 0, 0.25, 150, 0.7, 0.97, 300, 1, 0.7, 1, 0.99, 20, 0.2, 0, 0, 0, 1, 200, 1, 0, 0, 20, 0.99, 10, 1, 0, 3e6};

// P14=6 (по умолчанию): линейная оценка кандидатов вне потоков, веса LW (env SFW="w0 w1 ..." — переопределить);
// старое поведение (best fit + правила) — SF="14=0 15=4"
static double LW[32] = {1, -38.8, -22.4, 13.5, 47, -203, 107, -7.25, 2.78, 2.26, -8.96, 23.4, 1.05e+03, 69, 0, 46.7, 27.9, 8.7, 0, 0, 0, 0, 0, 0, 0, 20, 50};   // CMA-ES lin2 + LW25/26 (2026-10-19; LW24 «хвост с докопкой» ломал del_half_group/sawtooth)
// ---------------------------------------------------------------- аллокатор
// Свободные промежутки: граничные метки (склейка за O(1)) + корзины по точной длине (< NB) + set для больших.
static const int NB = 128;
static const int MAXL = 40;

struct Alloc {
    ll T = 0;                          // длина траншеи
    // tag[x] > 0 — в x начинается промежуток длины tag[x]; tag[x] < 0 — в x кончается промежуток,
    // начавшийся в -tag[x]-1; 0 — ни то ни другое (соседние промежутки всегда склеены, поэтому одно x
    // не бывает сразу началом и концом). owner[x] — PID трубы в клетке x (в свободных клетках — мусор).
    // Постраничное хранение: память растёт вместе с траншеей без копирования (худший случай T ≈ 8·10^6).
    struct Paged {
        static const int SH = 16, MASK = (1 << SH) - 1;
        std::vector<int*> pg; int fill;
        explicit Paged(int f) : fill(f) {}
        Paged(const Paged& o) : fill(o.fill) { for (int* p : o.pg) { int* q = new int[1 << SH]; memcpy(q, p, sizeof(int) << SH); pg.push_back(q); } }
        Paged& operator=(const Paged&) = delete;
        ~Paged() { for (int* p : pg) delete[] p; }
        void grow(ll n) { while ((ll)pg.size() << SH < n) { int* p = new int[1 << SH]; for (int i = 0; i <= MASK; i++) p[i] = fill; pg.push_back(p); } }
        int& operator[](ll i) { return pg[i >> SH][i & MASK]; }
        int operator[](ll i) const { return pg[i >> SH][i & MASK]; }
    };
    Paged tag{0}, owner{0}, zl{0};   // zl — длина-владелец зарезервированного участка (зоны, P80)
    std::set<int> bucket[NB];          // начала промежутков длины s (s < NB), по адресу
    std::set<std::pair<int,int>> big;  // (len, start) для len >= NB
    std::vector<int> pos, len;         // по PID
    int nPlaced = 0;

    // онлайн-статистика длин (с затуханием)
    double cnt[MAXL + 1] = {0}, cntTot = 0;

    void ensure(ll n) {
        tag.grow(n + 1); owner.grow(n + 1); zl.grow(n + 1);
    }
    void addGap(int s, int l) {
        tag[s] = l; tag[s + l] = -s - 1;
        if (l < NB) bucket[l].insert(s); else big.insert({l, s});
    }
    void delGap(int s) {
        int l = tag[s];
        if (l < NB) bucket[l].erase(s); else big.erase({l, s});
        tag[s] = 0; tag[s + l] = 0;
    }
    int freeRange(int s, int l) {   // возвращает конец склеенного промежутка
        int e = s + l;
        if (e < T && tag[e] > 0) { int r = tag[e]; delGap(e); e += r; }
        if (tag[s] < 0) { int ps = -tag[s] - 1; delGap(ps); s = ps; }
        addGap(s, e - s);
        return e;
    }
    void takeFrom(int gs, int d, int l) {
        int gl = tag[gs];
        delGap(gs);
        if (d > gs) addGap(gs, d - gs);
        if (d + l < gs + gl) addGap(d + l, gs + gl - (d + l));
    }
    void dig(ll n) {
        if (n <= 0) return;
        ll oldT = T;
        T += 20 * n;
        ensure(T);
        freeRange((int)oldT, (int)(T - oldT));
    }
    int tailStart() const { return (T > 0 && tag[T] < 0) ? -tag[T] - 1 : (int)T; }
    int pidLeftOf(int x) const { return x > 0 ? owner[x - 1] : 0; }      // труба, кончающаяся в x (x — начало дыры)
    int pidRightOf(int x) const { return x < T ? owner[x] : 0; }         // труба, начинающаяся в x (x — конец дыры)

    // ---- статистика
    void noteLifo(int pid, int L) { lastOfLen[L][lastPtr[L]] = pid; lastPtr[L] = (lastPtr[L] + 1) & 7; }
    void observe(int L) {
        for (int l = 1; l <= MAXL; l++) cntF[l] *= P[88];
        cntFTot = cntFTot * P[88] + 1; cntF[L] += 1;
        double dec = P[0];
        for (int l = 1; l <= MAXL; l++) cnt[l] *= dec;
        cntTot = cntTot * dec + 1; cnt[L] += 1;
    }
    // минимальная длина, которую «ещё ждём» (доля в недавних укладках > eps)
    int minRecent() const {
        double eps = P[1];
        for (int l = 1; l <= MAXL; l++) if (cnt[l] > eps * cntTot) return l;
        return 1;
    }

    // ---- заполнимость остатка: U[r] = u(r)·(средняя длина), u — вероятность обновления (какая-то
    // последовательность случайных труб из текущего распределения длин ляжет ровно в r). U ≈ 1 для «хороших» r,
    // 0 — для незаполнимых (короче всех длин или не представимых суммой частых длин, напр. 11 при длинах 5/10).
    static const int UR = 161;
    double U[UR] = {0};
    int uAge = 1 << 30;
    void rebuildU() {
        uAge = 0;
        if (cntTot <= 0) { for (int r = 0; r < UR; r++) U[r] = 1; return; }
        double p[MAXL + 1], mean = 0;
        for (int l = 1; l <= MAXL; l++) { p[l] = cnt[l] / cntTot; if (p[l] < P[71]) p[l] = 0; }
        double z = 0; for (int l = 1; l <= MAXL; l++) z += p[l];
        if (z <= 0) { for (int r = 0; r < UR; r++) U[r] = 1; return; }
        for (int l = 1; l <= MAXL; l++) { p[l] /= z; mean += l * p[l]; }
        double u[UR]; u[0] = 1;
        for (int r = 1; r < UR; r++) {
            double v = 0;
            for (int l = 1; l <= MAXL && l <= r; l++) v += p[l] * u[r - l];
            u[r] = v;
        }
        for (int r = 0; r < UR; r++) U[r] = std::min(1.0, u[r] * mean);
    }
    double fillU(int r) const { return r <= 0 ? 1 : r >= UR ? 1 : U[r]; }

    // ---- модель ожидаемых потерь дыры: W[r] — сколько метров дыры длины r останется пустыми,
    // если её заполнять случайными трубами из текущего распределения длин
    double W[NB] = {0};
    int wAge = 1 << 30;
    void rebuildW() {
        double pr[MAXL + 1], tot = 0;
        for (int l = 1; l <= MAXL; l++) { pr[l] = cnt[l] > P[1] * cntTot ? cnt[l] : 0; tot += pr[l]; }
        W[0] = 0;
        for (int r = 1; r < NB; r++) {
            double s = 0, w = 0;
            for (int l = 1; l <= MAXL && l <= r; l++) if (pr[l] > 0) { s += pr[l]; w += pr[l] * W[r - l]; }
            W[r] = s > 0 ? w / s : r;
        }
        wAge = 0;
    }
    double Wf(int r) const { return r < NB ? W[r] : 0; }

    // стоимость остатка r после укладки (r > 0)
    double fragCost(int r, int m) const {
        if (r < m) return r * P[2];        // мёртвый остаток
        return r * P[3];                                           // живой остаток (P3=0 — без разницы)
    }

    // ---- потоки (сегрегация по признаку, чьи носители умирают вместе)
    // stream[pid]: -1 — обычная, 0 — поток A (кладётся к левому краю граничного промежутка),
    // 1 — поток B (к правой границе роста F, растёт влево)
    std::vector<signed char> stream;
    int F = 0;                        // правая граница роста; F == T — обычный режим
    bool flip = false;                // P107: в текущем участке B растёт слева, A — справа
    double bVol = 0;                  // затухающий объём укладок потока B (для размера участка)
    double fbVol = 0;                 // то же, но только укладки B к границе роста
    std::vector<unsigned char> feat;  // биты признаков при рождении: bit0 — чётность индекса в серии равных длин, bit1 — чётность PID
    int runIdx = 0, lastL = 0, curRun = 0;
    std::vector<int> runOf;           // номер серии равных длин по PID
    static const int NF = 5;   // + bit2..4: индекс в серии mod 3 == 0 / 1 / 2
    double coSame[NF] = {0}, coTot[NF] = {0};   // соседние удаления в пачке: совпал ли признак
    double birthV[NF][2] = {{0}}, deathV[NF][2] = {{0}};
    int lastDelPid = 0; bool lastWasDel = false;
    int curFeat = -1;
    double tieAge = 0;
    double labSame = 0, labTot = 0;   // соседние удаления помеченных труб: совпал ли поток (подтверждение потоков)
    bool streamsVerified() const { return labTot >= 20 && labSame / labTot >= P[41]; }
    double bornLab[2] = {0, 0}, diedLab[2] = {0, 0};   // проверка прогноза: смерти среди помеченных A/B
    // ---- риск смерти соседа по возрасту
    // aliveBit — дерево Фенвика по PID (живые), deathAge[b] — число смертей в возрасте [2^b-1, 2^(b+1)-1)
    std::vector<int> aliveBit;
    double deathAge[24] = {0};
    void bitAdd(int i, int v) { for (; i < (int)aliveBit.size(); i += i & -i) aliveBit[i] += v; }
    int bitSum(int i) const { int s = 0; for (; i > 0; i -= i & -i) s += aliveBit[i]; return s; }
    static int ageBucket(int a) { int b = 0; while ((2 << b) - 1 <= a && b < 23) b++; return b; }
    // P(труба PID p когда-нибудь умрёт | дожила до текущего возраста)
    double hazard(int p) const {
        if (p <= 0) return 0;   // край траншеи
        int cur = (int)pos.size() - 1;
        int a = cur - p, b = ageBucket(a);
        double dGE = 0; for (int k = b; k < 24; k++) dGE += deathAge[k];
        int lo = cur - ((1 << b) - 1);   // живые с возрастом >= 2^b-1
        double aGE = lo > 0 ? bitSum(lo) : 0;
        return dGE + aGE > 0 ? dGE / (dGE + aGE) : 0.5;
    }
    // локальная частота недавних смертей по корзинам PID (окно удалений, ползущее по PID — как в 15)
    std::vector<double> bdv; std::vector<int> bdt, aliveB; int nDel = 0;
    void bEnsure(int b) { if ((int)bdv.size() <= b) { bdv.resize(b + 1, 0); bdt.resize(b + 1, 0); aliveB.resize(b + 1, 0); } }
    double rateB(int p) const {
        if (p <= 0) return 0;
        int b = p >> (int)P[48];
        if (b >= (int)bdv.size()) return 0;
        double v = bdv[b] * std::pow(P[49], nDel - bdt[b]);
        return v / std::max(1, aliveB[b]);
    }
    void bOnBirth(int p) { int b = p >> (int)P[48]; bEnsure(b); aliveB[b]++; }
    void bOnDeath(int p) {
        int b = p >> (int)P[48]; bEnsure(b); aliveB[b]--; nDel++;
        bdv[b] = bdv[b] * std::pow(P[49], nDel - bdt[b]) + 1; bdt[b] = nDel;
    }
    double birthL[MAXL + 1] = {0}, deathL[MAXL + 1] = {0};
    double pdieLen(int p) const { if (p <= 0) return 0; int l = len[p]; return (deathL[l] + 0.5) / (birthL[l] + 1); }
    // интенсивность смертей по длине: доля недавних удалений длины l / доля живых длины l
    double recDel[MAXL + 1] = {0}, recDelTot = 0; int aliveCntL[MAXL + 1] = {0}, aliveTot = 0;
    double hzRate(int p) const {
        if (p <= 0 || recDelTot < 1) return 0;
        int l = len[p];
        double rel = (recDel[l] / recDelTot) / ((double)std::max(1, aliveCntL[l]) / std::max(1, aliveTot));
        return std::max(std::min(1.0, rel / P[30]), modRate(p));
    }
    // P61 >= 2: интенсивность смертей по классу PID mod m (m = 2..P61): доля недавних удалений / доля живых, max по m
    static const int MODMAX = 16;
    double recDelMod[MODMAX + 1][MODMAX] = {}; int aliveMod[MODMAX + 1][MODMAX] = {};
    double modRate(int p) const {
        if (p <= 0 || recDelTot < 1 || P[61] < 2) return 0;
        double best = 0;
        for (int m = 2; m <= (int)P[61]; m++) {
            int r = p % m;
            double rel = (recDelMod[m][r] / recDelTot) / ((double)std::max(1, aliveMod[m][r]) / std::max(1, aliveTot));
            if (rel > best) best = rel;
        }
        return std::min(1.0, best / P[30]);
    }
    double hz2(int p) const { return P[18] == 4 ? hzRate(p) : P[18] == 1 ? pdieLen(p) : P[18] == 2 ? std::max(pdieLen(p), hazard(p)) : P[18] == 3 ? pdieLen(p) * hazard(p) : hazard(p); }
#ifdef SIM
    unsigned long long rs = 88172645463325252ULL; int epsRand = 0;
    int forcePid = -1, forceAlt = 0, listPid = -1, optBest = 0; std::vector<int> optList;
    unsigned long long rng() { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
#endif
    const int* oracleDeath = nullptr;   // только SIM-эксперименты: время смерти по PID (команда), -1 — не умрёт
    int curCmd = 0;

    double kappa(int f) const {
        if (coTot[f] < P[21]) return 0;
        double b0 = birthV[f][0], b1 = birthV[f][1], bt = b0 + b1;
        if (bt <= 0) return 0;
        if (std::min(b0, b1) < P[36] * bt) return 0;   // признак почти константа — каппа вырождена
        double e = (b0 * b0 + b1 * b1) / (bt * bt);   // ожидаемая доля совпадений при независимости
        double o = coSame[f] / coTot[f];
        return e >= 1 ? 0 : (o - e) / (1 - e);
    }
    double pdie(int f, int v) const { return (deathV[f][v] + 1) / (birthV[f][v] + 2); }
    // P99: тень без потоков (та же история команд, своя раскладка) — для выключения потоков, если они проигрывают
    bool noStreams = false, streamsOff = false, shadowDropped = false; Alloc* shadow = nullptr; ll occ = 0;
    ll freeIn() const {   // свободное внутри траншеи, кроме хвоста и граничного промежутка потоков
        ll f = T - occ - (T - tailStart());
        if (F < T && F > 0 && tag[F] < 0) f -= F - (-tag[F] - 1);
        return f;
    }
    void shadowCheck();
    // P103: доля помеченных (A/B) среди недавних удалений (считаются, пока живы помеченные) против их доли среди живых
    double delB = 0, delN = 0; int aliveB1 = 0;
    bool bFalsified() const { return P[103] > 0 && delN >= P[104] && aliveB1 > 0 && delB / delN < P[103] * aliveB1 / std::max(1, aliveTot); }
    int chooseStream(int pid, int L) {
        unsigned char fv = feat[pid];
        if (P[10] == 1 && oracleDeath) {
            if (L < P[62] || (P[63] > 0 && L >= P[63])) return -1;   // оракул только для длин [P62, P63)
            int dt = oracleDeath[pid];
            return dt >= 0 && dt - curCmd <= P[11] ? 1 : 0;
        }
        if (P[10] != 2 || noStreams || streamsOff) return -1;
        // признак: каппа «умирают вместе» выше порога; при равенстве — тот, где сильнее разница p_die
        // гистерезис: текущий признак меняем, только если другой лучше с запасом
        int bf = -1; double bs = -1, cs = -1;
        for (int f = 0; f < NF; f++) {
            double k = kappa(f);
            if (k < P[22]) continue;
            double s = k + std::fabs(pdie(f, 0) - pdie(f, 1));
            if (f == curFeat) cs = s;
            if (s > bs) { bs = s; bf = f; }
        }
        if (bf < 0) { curFeat = -1; return -1; }
        // P64 > 0: фальсификация — помечено >= P64 труб, а умерло из них < P65 доли: прогноз не сбывается
        if (P[64] > 0 && bornLab[0] + bornLab[1] >= P[64] && diedLab[0] + diedLab[1] < P[65] * (bornLab[0] + bornLab[1])) {
            curFeat = -1; return -1;
        }
        if (cs >= 0 && bs < cs + 0.2) bf = curFeat;
        curFeat = bf;
        int v = (fv >> bf) & 1;
        // тот, кто умирает чаще, — в поток B (правый край, у границы роста)
        double pv = pdie(bf, v), po = pdie(bf, v ^ 1);
        if (std::fabs(pv - po) < 0.1) return v;
        return pv > po ? 1 : 0;
    }

    bool kin(int p, int st) const { return p > 0 && stream[p] == st && runOf[p] == curRun; }

    // положить в промежуток gs длины s: выбрать сторону, вернуть D
    int chooseSide(int gs, int s, int L, int st) {
        int r = s - L;
        if (r == 0) return gs;
        int lp = pidLeftOf(gs), rp = pidRightOf(gs + s);
        if (st >= 0) {
            if (P[29] != 0) {
                bool kl = kin(lp, st), kr = kin(rp, st);
                if (kl && !kr) return gs;
                if (kr && !kl) return gs + r;
            }
            // прижимаем к соседу из того же потока
            int lc = lp ? stream[lp] : -2, rc = rp ? stream[rp] : -2;
            if (lc == st && rc != st) return gs;
            if (rc == st && lc != st) return gs + r;
            return st == 1 ? gs + r : gs;
        }
        if (P[8] == 4) { if (hz2(lp) > hz2(rp)) return gs + r; return gs; }
        if (P[8] == 5) {
            double a = rateB(lp), b = rateB(rp);
            if (a > P[50] * b && a > 0) return gs + r;   // остаток — к соседу с большей локальной частотой смертей
            if (b > P[50] * a && b > 0) return gs;
            if (lp > rp) return gs + r;
            return gs;
        }
        bool rOld;
        if (P[8] == 1 && P[92] >= 2 && lp > 0 && rp > 0 && raceSays(lp, rp, rOld)) return (rOld == (lp < rp)) ? gs + r : gs;
        if (P[8] == 1) {
            if (P[56] != 0 && oldDieFirst()) {   // P56: умирают скорее старые (очередь) — остаток к старшему соседу
                if (lp == 0) return gs;
                if (rp == 0) return gs + r;
                if (lp < rp) return gs + r;
                return gs;
            }
            if (lp > rp) return gs + r;
        }
        else if (P[8] == 3) return gs + r;
        return gs;
    }
    // средний нормированный ранг удаляемой трубы среди живых (0 — самая молодая, 1 — самая старая), с затуханием
    double qSum = 0, qTot = 0;
    // P92: «гонка соседей» — у дыр, куда клали с остатком, какой из двух соседей (старший/младший) умер первым
    // (затухание P93); если наблюдений >= 20 — сторона остатка по этой статистике вместо qSum
    double raceOld = 0, raceTot = 0;
    std::vector<int> raceHead, raceNext, raceOther;   // по PID — список записей; запись = пара (сосед, другой сосед)
    std::vector<char> raceDone;
    // ключ записи (P92): 1 — один на всех; 2 — корзина возраста старшего соседа; 3 — (корзина старшего, корзина младшего);
    // 4 — длины соседей (короткая/длинная относительно медианы недавних укладок) у старшего и младшего
    std::vector<int> raceKey;
    double raceOldK[600] = {0}, raceTotK[600] = {0};
    int raceKeyOf(int a, int b) const {   // a — старший (меньший PID)
        int cur = (int)pos.size() - 1;
        int ba = ageBucket(cur - a), bb = ageBucket(cur - b);
        if (P[92] == 2) return 1 + ba;
        if (P[92] == 3) return 1 + ba * 24 + bb;
        if (P[92] == 4) return 1 + (len[a] >= 8) * 2 + (len[b] >= 8);
        if (P[92] == 5) return 1 + (ba - bb > 2 ? 2 : ba - bb);
        return 0;
    }
    bool raceSays(int lp, int rp, bool& older) const {   // true — есть статистика; older — остаток к старшему
        int a = lp < rp ? lp : rp, b = lp < rp ? rp : lp;
        int k = raceKeyOf(a, b);
        if (k > 0 && k < 600 && raceTotK[k] >= P[94] && std::fabs(raceOldK[k] / raceTotK[k] - 0.5) > P[95]) { older = raceOldK[k] > 0.5 * raceTotK[k]; return true; }
        if (raceTot >= 20 && std::fabs(raceOld / raceTot - 0.5) > P[95]) { older = raceOld > 0.5 * raceTot; return true; }
        return false;
    }
    void raceAdd(int a, int b) {   // a, b — живые соседи дыры
        if ((int)raceHead.size() < (int)pos.size() + 1) raceHead.resize(pos.size() + 1024, -1);
        int id = (int)raceDone.size(); raceDone.push_back(0);
        raceKey.push_back(raceKeyOf(a < b ? a : b, a < b ? b : a));
        int pr[2] = {a, b};
        for (int k = 0; k < 2; k++) { raceNext.push_back(raceHead[pr[k]]); raceOther.push_back(pr[k ^ 1]); raceHead[pr[k]] = 2 * id + k; }
    }
    void raceOnDeath(int p) {
        if (p >= (int)raceHead.size()) return;
        for (int e = raceHead[p]; e >= 0; e = raceNext[e]) {
            int id = e >> 1;
            if (raceDone[id]) continue;
            raceDone[id] = 1;
            int o = p < raceOther[e] ? 1 : 0;
            raceOld = raceOld * P[93] + o; raceTot = raceTot * P[93] + 1;
            int k = raceKey[id];
            if (k > 0 && k < 600) { raceOldK[k] = raceOldK[k] * P[93] + o; raceTotK[k] = raceTotK[k] * P[93] + 1; }
        }
        raceHead[p] = -1;
    }
    bool oldDieFirst() const {
        if (P[92] != 0 && raceTot >= 20 && std::fabs(raceOld / raceTot - 0.5) > P[95]) return raceOld > 0.5 * raceTot;
        return qTot >= 20 && qSum / qTot > 0.5;
    }
    int sideNb(int lp, int rp) const {   // сосед, к которому ляжет остаток (вне потоков)
        bool rOld;
        if (P[92] >= 2 && lp > 0 && rp > 0 && raceSays(lp, rp, rOld)) return rOld ? (lp < rp ? lp : rp) : (lp > rp ? lp : rp);
        if (P[56] != 0 && oldDieFirst()) { if (lp == 0) return rp; if (rp == 0) return lp; return lp < rp ? lp : rp; }
        return lp > rp ? lp : rp;
    }

    std::pair<ll,ll> finishPlace(int pid, int L, int D, ll N, int st, bool useGap, int bs) {
        (void)bs;
#ifdef SIM
        if (getenv("SFPL") && pid >= atoi(getenv("SFPL")) && pid < atoi(getenv("SFPL")) + 30)
            fprintf(stderr, "place pid=%d L=%d st=%d D=%d N=%lld gap=%s bs=%d F=%d T=%lld\n", pid, L, st, D, N, useGap ? "hole" : "front", bs, F, T);
#endif
        for (int i = 0; i < L; i++) owner[D + i] = pid;
        fbVol = fbVol * P[38] + ((st == 1 && !useGap) ? L : 0);
        bOnBirth(pid);
        if (st >= 0) aliveB1++;
        pos.push_back(D); len.push_back(L); alive.push_back(1); birthL[L] += 1; aliveCntL[L]++; aliveTot++;
#ifdef SIM
        abl[L].push_back((int)pos.size() - 1);
        if (!inRoll && gTPv[0]) tpPushP(L);
#endif
        if (P[61] >= 2) for (int m = 2; m <= (int)P[61]; m++) aliveMod[m][pid % m]++;
        // после расширения rebuildBit уже учёл pid (alive[pid] = 1) — второй bitAdd давал «призрачных» живых
        if ((int)aliveBit.size() <= pid) { aliveBit.resize(std::max<size_t>(pid + 1, aliveBit.size() * 2 + 16), 0); rebuildBit(); }
        else bitAdd(pid, 1);
        noteLifo(pid, L); observe(L); nPlaced++;
        lastWasDel = false;
        return {N, D};
    }

    // ---- зоны по длинам (P80): трубы частой длины — в свои участки по порядку. Генератор открытых тестов
    // удаляет почти всегда одну из последних труб той же длины, поэтому в зоне дыры ровно под следующую трубу этой
    // длины, а умершие подряд сливаются в куски, кратные длине. Работает, пока активных длин мало.
    ll resEnd = 0;
    std::vector<int> slots[MAXL + 1];
    int nActive() const { int k = 0; for (int l = 1; l <= MAXL; l++) if (cntTot > 0 && cnt[l] > P[82] * cntTot) k++; return k; }
    bool zonesEnabled() const {
        if (!(P[80] > 0 && cntTot >= 50 && !(P[13] > 0 && T < P[13]))) return false;
        int k = nActive(); double mx = 0;
        for (int l = 1; l <= MAXL; l++) mx = std::max(mx, cnt[l] / cntTot);
        if (!(k >= 4 && k <= P[83] && (P[85] <= 0 || mx <= P[85]))) return false;
        if (P[91] > 0 && (lifoTot < 20 || lifoHit < P[91] * lifoTot)) return false;   // удаления не «по стеку длины»
        if (P[87] > 0) {   // смена режима: среди недавних укладок (быстрое затухание) мало длин с зонами — выключить
            double f = 0;
            for (int l = 1; l <= MAXL; l++) if (cnt[l] > P[82] * cntTot) f += cntF[l];
            if (f < P[87] * cntFTot) return false;
        }
        return true;
    }
    double cntF[MAXL + 1] = {0}, cntFTot = 0;   // статистика длин с быстрым затуханием (P88)
    // «по стеку своей длины»: удаляемая труба — среди 8 последних уложенных той же длины (доля, затухание 0.99)
    int lastOfLen[MAXL + 1][8] = {{0}}; int lastPtr[MAXL + 1] = {0};
    double lifoHit = 0, lifoTot = 0;
    bool zoneOn(int z, bool en) const { return en && z > 0 && cnt[z] > P[82] * cntTot; }
    bool occupied(ll x) const { int p = owner[x]; return p > 0 && alive[p] && pos[p] <= x && x < pos[p] + len[p]; }
    // свободна ли [d, d+L) от чужих зарезервированных клеток
    bool clearOf(int d, int L, int self, bool en) const {
        for (int x = d; x < d + L; x++) { int z = zl[x]; if (z && z != self && zoneOn(z, en)) return false; }
        return true;
    }
    bool zoneTry(int L, int st, int m, ll& N, int& D) {
        if (st >= 0) return false;
        bool en = zonesEnabled();
        if (!en) return false;
        int ts = tailStart();
        bool mine = zoneOn(L, en) && aliveCntL[L] >= P[86];
        N = 0; D = -1;
        // точное попадание: в свою зону или в незарезервированную дыру
        { int k = 0; for (int g : bucket[L]) { if (++k > 64) break; int z = zl[g]; if (z == L || !zoneOn(z, en) || P[84] > 0) { D = g; takeFrom(g, g, L); return true; } } }
        if (mine) {
            // свободный слот своей зоны (стек: последний освободившийся — первым)
            while (!slots[L].empty()) {
                int d = slots[L].back(); slots[L].pop_back();
                if (zl[d] != L || occupied(d) || occupied(d + L - 1)) continue;
                ll hs = d; int steps = 0;
                while (hs > 0 && !occupied(hs - 1) && steps < 4096) { hs--; steps++; }
                if (steps >= 4096 || tag[hs] <= 0 || hs + tag[hs] < d + L) continue;
                D = d; takeFrom((int)hs, D, L); return true;
            }
            // P90: сначала — участок внутри незарезервированной дыры (>= 4 слотов; участок ужимается под дыру),
            // иначе незарезервированные дыры по правилу редкой длины (ниже), и лишь потом новый участок в хвосте
            if (P[90] > 0) {
                int hg = -1, hl = 0;
                { int k = 0; for (int s2 = 4 * L; s2 < NB && hg < 0; s2++) { k = 0; for (int g : bucket[s2]) { if (++k > 64) break; if (g != ts && !zoneOn(zl[g], en) && clearOf(g, s2, L, en)) { hg = g; hl = s2; break; } } } }
                if (hg < 0) { int k = 0; for (auto& q : big) { if (++k > 64) break; if (q.second != ts && !zoneOn(zl[q.second], en) && clearOf(q.second, std::min<int>(q.first, (int)P[81] * L), L, en)) { hg = q.second; hl = q.first; break; } } }
                if (hg >= 0) {
                    ll C2 = std::min<ll>((ll)P[81] * L, (ll)(hl / L) * L);
                    for (ll x = hg; x < hg + C2; x++) zl[x] = L;
                    for (ll x = hg + C2 - L; x > hg; x -= L) slots[L].push_back((int)x);
                    D = hg; takeFrom(hg, D, L);
                    return true;
                }
                if (rarePlace(L, m, ts, en, D)) return true;
            }
            // новый участок из P81 слотов в конце траншеи, за всеми зарезервированными
            ll C = (ll)P[81] * L;
            ll stp = std::max<ll>(ts, resEnd);
            if (T < stp + C) { N = (stp + C - T + 19) / 20; dig(N); }
            for (ll x = stp; x < stp + C; x++) zl[x] = L;
            for (ll x = stp + C - L; x > stp; x -= L) slots[L].push_back((int)x);
#ifdef SIM
            if (getenv("ZDBG")) fprintf(stderr, "chunk L=%d at %lld T=%lld nact=%d cnt=%.3f alive=%d\n", L, stp, T, nActive(), cnt[L] / cntTot, aliveCntL[L]);
#endif
            resEnd = stp + C;
            D = (int)stp; takeFrom(tailStart(), D, L);
            return true;
        }
        if (rarePlace(L, m, ts, en, D)) return true;
        ll stp = std::max<ll>(ts, resEnd);
        if (T < stp + L) { N = (stp + L - T + 19) / 20; dig(N); }
        D = (int)stp; takeFrom(tailStart(), D, L);
        if (resEnd < stp + L && resEnd > ts) resEnd = stp + L;
        return true;
    }
    // редкая длина: остаток >= m самого редкого размера среди незарезервированных дыр (не хвост); false — нет такой
    bool rarePlace(int L, int m, int ts, bool en, int& D) {
        double bc = 1e18; int bg = -1;
        for (int s = L + m; s < NB; s++) { int k = 0; for (int g : bucket[s]) {
            if (++k > 64) break;
            if (g == ts || zoneOn(zl[g], en) || !clearOf(g, L, L, en)) continue;
            double c = (double)bucketSize(s - L) * 1000 + s;
            if (c < bc) { bc = c; bg = g; }
            break;
        } }
        { int k = 0; for (auto& q : big) { if (++k > 64) break; if (q.second == ts || zoneOn(zl[q.second], en) || !clearOf(q.second, L, L, en)) continue; double c = P[73] * 1000 + q.first; if (c < bc) { bc = c; bg = q.second; } break; } }
        if (bg >= 0) { D = bg; takeFrom(bg, bg, L); return true; }
        return false;
    }

#ifdef SIM
    int rollH = 0, rollK = 4, rollMode = 1, rollMin = 1; double accMin = 0.9, haccMin = 0.0; bool inRoll = false; double rollMargin = 0;
    std::vector<int> abl[41];   // живые PID каждой длины по возрастанию
    int lastVictim = 0;
    Alloc* cloneNoShadow() const { Alloc* c = new Alloc(*this); c->shadow = nullptr; c->shadowDropped = true; c->inRoll = true; c->forcePid = -1; c->listPid = -1; return c; }
    std::pair<ll,ll> place(int L) {
        if (rollH > 0 && !inRoll && T >= 60) {
            int pid = (int)pos.size();
            Alloc* pr = cloneNoShadow(); pr->listPid = pid; pr->placeCore(L);
            std::vector<int> opts = pr->optList; int best = pr->optBest; delete pr;
            if (opts.size() >= 2) {
                std::vector<int> cand; cand.push_back(best);
                for (int o : opts) if (o != best && (int)cand.size() < rollK) cand.push_back(o);
                std::vector<int> ptoks; int view = 0;
                if (rollMode == 2 && gTPv[0]) {
                    for (int v = 1; v < 3; v++) if (gAcc[v] > gAcc[view] + 1e-9) view = v;
                    if (gAcc[view] < accMin) return placeCore(L);
                    TokPred* tp = gTPv[view];
                    std::vector<int> h(tp->hist.end() - std::min<size_t>(tp->hist.size(), TokPred::MAXO), tp->hist.end());
                    h.push_back(L);
                    for (int t = 0; t < rollH; t++) { int p = tp->predict(h); if (p < 0) break; ptoks.push_back(p); h.push_back(p); if ((int)h.size() > TokPred::MAXO) h.erase(h.begin()); }
                    if ((int)ptoks.size() < rollMin) return placeCore(L);
                    if (getenv("HDBG") && (pid % 100) == 0) fprintf(stderr, "pid=%d view=%d acc=%.3f/%.3f/%.3f hacc=%.3f np=%zu\n", pid, view, gAcc[0], gAcc[1], gAcc[2], gHacc, ptoks.size());
                    if (gHacc < haccMin) { gPend.assign(ptoks.begin(), ptoks.begin() + std::min<size_t>(2, ptoks.size())); gPendIdx = 0; gPendView = view; gSkip = 1; ptoks.clear(); return placeCore(L); }
                    gPend.assign(ptoks.begin(), ptoks.begin() + std::min<size_t>(2, ptoks.size())); gPendIdx = 0; gPendView = view; gSkip = 1;
                }
                double bestScore = 1e300; int bestAlt = best;
                for (size_t ci2 = 0; ci2 < cand.size(); ci2++) {
                    int o = cand[ci2];
                    Alloc* c = cloneNoShadow();
                    if (o != best) { c->forcePid = pid; c->forceAlt = o; }
                    c->placeCore(L);
                    ll emax = c->tailStart();
                    if (rollMode == 1) {
                        for (int t = 1; t <= rollH; t++) {
                            long long idx = (long long)curCmd + t; if (idx >= (long long)gCmds->size()) break;
                            int cmd = (*gCmds)[idx]; if (cmd == 0) break;
                            c->curCmd = (int)idx;
                            if (cmd > 0) c->placeCore(cmd); else c->remove(-cmd);
                            ll e = c->tailStart(); if (e > emax) emax = e;
                        }
                    } else {
                        for (size_t t = 0; t < ptoks.size(); t++) {
                            int tk = ptoks[t]; c->curCmd = curCmd + 1 + (int)t;
                            if (tk < 100) c->placeCore(tk);
                            else {
                                int victim = 0;
                                if (view == 0) { int L2 = (tk - 100) / 16, rk = (tk - 100) % 16; std::vector<int>& v = c->abl[L2]; if ((int)v.size() <= rk) break; victim = v[v.size() - 1 - rk]; }
                                else if (view == 1) victim = c->lastVictim + (tk - 100000 - 5000);
                                else victim = (int)c->pos.size() - 1 - (tk - 200000);
                                if (victim <= 0 || victim >= (int)c->pos.size() || !c->alive[victim]) break;
                                c->remove(victim);
                            }
                            ll e = c->tailStart(); if (e > emax) emax = e;
                        }
                    }
                    double sc = (double)c->T * 1e9 + (double)emax * 1e3 + (double)std::min<ll>(c->freeIn(), 999);
                    if (ci2 == 0) sc -= rollMargin;   // небольшое предпочтение умолчанию
                    delete c;
                    if (sc < bestScore) { bestScore = sc; bestAlt = o; }
                }
                if (bestAlt != best) { forcePid = pid; forceAlt = bestAlt; }
            }
        }
        return placeCore(L);
    }
    std::pair<ll,ll> placeCore(int L) {
#else
    std::pair<ll,ll> place(int L) {
#endif
        if (P[99] > 0 && !noStreams) shadowCheck();   // до первой возможной метки потока тень совпала бы с основной
        if (shadow) { shadow->curCmd = curCmd; shadow->place(L); }
        occ += L;
        int m = minRecent();
        if (++uAge >= 64) rebuildU();
        int pid = (int)pos.size();
        if (L == lastL) runIdx++; else { runIdx = 0; curRun++; }
        runOf.push_back(curRun);
        lastL = L;
        feat.push_back((unsigned char)((runIdx & 1) | ((pid & 1) << 1) | (P[44] ? (1 << (2 + runIdx % 3)) : 0)));
        int st = chooseStream(pid, L);
        stream.push_back((signed char)st);
        bVol = bVol * P[27] + (st == 1 ? L : 0);
        if (st >= 0) bornLab[st] += 1;
#ifdef SIM
        if (getenv("SFDBG") && pid % atoi(getenv("SFDBG")) == 0)
            fprintf(stderr, "pid=%d L=%d feat=%d st=%d k0=%.3f k1=%.3f tot=%.1f,%.1f F=%d T=%lld\n", pid, L, feat[pid], st, kappa(0), kappa(1), coTot[0], coTot[1], F, T);
#endif
        for (int f = 0; f < NF; f++) {
            birthV[f][0] *= P[23]; birthV[f][1] *= P[23];
            deathV[f][0] *= P[23]; deathV[f][1] *= P[23];
            birthV[f][(feat[pid] >> f) & 1] += 1;
        }
        if (F > T) F = (int)T;
        int ts = tailStart();
        if (P[14] == 5 && ++wAge >= 64) rebuildW();
        double selfHz = 0;
        if (P[18] == 4 && recDelTot >= 1)
            selfHz = std::max(std::min(1.0, (recDel[L] / recDelTot) / ((double)std::max(1, aliveCntL[L]) / std::max(1, aliveTot)) / P[30]), modRate(pid));
        if (P[80] > 0) { ll zN = 0; int zD = -1; if (zoneTry(L, st, m, zN, zD)) return finishPlace(pid, L, zD, zN, st, false, 0); }
        double bestC = 1e18; int bgs = -1, bs = 0;
        int M = P[4] > 0 ? (int)P[4] : 1;  // сколько промежутков каждой длины смотреть
        // P53 > 0: последние P53 команд до лимита 200000 (P54 — лимит, если задан) — будущего почти нет
        bool endgame = P[53] > 0 && curCmd >= (P[54] > 0 ? P[54] : 200000) - P[53];
        bool verifiedB = P[12] == 3 && bornLab[1] >= 20 && bornLab[0] >= 20 &&
            diedLab[1] / bornLab[1] - diedLab[0] / bornLab[0] > P[9];
        bool bToFront = st == 1 && verifiedB;
        auto eval = [&](int gs, int s) {
            double cst;
            if (P[13] > 0 && T < P[13]) {
                // крошечная траншея: точное попадание (если сосед не «умирающий»), иначе самая большая дыра
                if (s == L) cst = -1000 + P[32] * std::max(pdieLen(pidLeftOf(gs)), pdieLen(pidRightOf(gs + s)));
                else cst = -s;
            } else if (P[14] == 1) {
                cst = gs;                                   // first fit по адресу
            } else if (P[14] == 2) {
                cst = -s;                                   // worst fit
            } else if (P[14] == 3) {
                cst = s == L ? -1e9 : -s;                   // точное, иначе worst
            } else if (P[14] == 5) {
                // изменение ожидаемых потерь + разрыв ничьих по best fit
                cst = Wf(s - L) - Wf(s) + 1e-3 * (s - L);
            } else if (P[14] == 4) {
                cst = -(double)gs;                          // last fit
            } else if (P[14] == 6 && st < 0) {
                int r = s - L;
                int lp = pidLeftOf(gs), rp = pidRightOf(gs + s);
                double hl = hazard(lp), hr = hazard(rp);
                double hside = sideNb(lp, rp) == lp ? hl : hr;   // остаток — к соседу, который скорее умрёт
                bool dead = r > 0 && r < m;
                double dr = (r >= 1 && r <= MAXL && cntTot > 0) ? cnt[r] / cntTot : 0;
                double ds = (s <= MAXL && cntTot > 0) ? cnt[s] / cntTot : 0;
                double hzn = std::max(hzRate(lp), hzRate(rp));
                cst = LW[0] * r + LW[1] * (r == 0) + LW[2] * dead + LW[3] * (dead ? r : 0) + LW[4] * (r >= m)
                    + LW[5] * dr + LW[6] * ds + LW[7] * (r > 0 ? hside : 0) + LW[8] * (dead ? r * (1 - std::max(hl, hr)) : 0)
                    + LW[9] * (r == 0 ? hzn : 0) + LW[10] * std::log(1.0 + (double)bucketSize(s)) + LW[11] * std::log((double)s / L)
                    + LW[12] * (dead && r >= 4) + LW[15] * (r == 0 ? std::min(hl, hr) : 0) + LW[16] * (r > 0 && r < 2 * m && !dead);
                if (r > 0 && LW[25] != 0 && cntTot >= P[89]) {   // изменение избытка дыр над ожидаемым спросом за LW[26] укладок: +r, −s
                    auto dem = [&](int z) { return z > MAXL ? 1e9 : (cntTot > 0 ? cnt[z] / cntTot * LW[26] : 0); };
                    double exr = (double)bucketSize(r) + 1 - dem(r), exs = (double)bucketSize(s) - dem(s);
                    cst += LW[25] * (std::log(1.0 + std::max(0.0, exr)) - std::log(1.0 + std::max(0.0, exs)));
                }
                if (r >= m && (LW[22] != 0 || LW[23] != 0)) {   // живой остаток: редкость его размера среди дыр
                    if (r <= MAXL) cst += LW[22] * std::log(1.0 + (double)bucketSize(r));
                    else cst += LW[23];
                }
                if (r > 0 && (LW[20] != 0 || LW[21] != 0)) {   // незаполнимый остаток
                    double nf = 1 - fillU(r);
                    cst += LW[20] * nf + LW[21] * nf * r;
                }
                if (LW[18] != 0 && r > 0) cst += LW[18] * ds / (1.0 + (double)bucketSize(s));   // ломаем дефицитный размер
                if (LW[19] != 0 && r > 0) cst += LW[19] * dr / (1.0 + (double)bucketSize(r));   // остаток — дефицитного размера
                if (r == 0 && LW[17] != 0) {   // исходный штраф «замуровывания» (как P17/P39)
                    double h = 0;
                    int nb[2] = {lp, rp};
                    for (int k = 0; k < 2; k++) {
                        if (nb[k] <= 0) continue;
                        int mg = s + len[nb[k]];
                        if (mg > MAXL || cnt[mg] < P[39] * cntTot) continue;
                        h = std::max(h, hzRate(nb[k]));
                    }
                    h -= selfHz;
                    if (h >= 0.5) cst += LW[17] * h;
                }
            } else {
                cst = s == L ? -P[5] : fragCost(s - L, m) + P[6] * (s - L);
                if (P[17] != 0 && s == L) {
                    // точное попадание «замуровывает» дыру у соседа, который скоро умрёт (дыра бы выросла)
                    double h;
                    if (P[39] > 0) {
                        // дыра вырастет после смерти соседа до s + len(соседа) — полезно, только если такая длина
                        // часто встречается среди недавних укладок
                        h = 0;
                        int nb[2] = {pidLeftOf(gs), pidRightOf(gs + s)};
                        for (int k = 0; k < 2; k++) {
                            if (nb[k] <= 0) continue;
                            int mg = s + len[nb[k]];
                            if (mg > MAXL || cnt[mg] < P[39] * cntTot) continue;
                            h = std::max(h, hz2(nb[k]));
                        }
                    } else h = std::max(hz2(pidLeftOf(gs)), hz2(pidRightOf(gs + s)));
                    if (P[18] == 4) h -= selfHz;   // вредно, только если новая труба живёт дольше соседа
                    if (h >= P[28]) cst += P[17] * h;
                }
                if (P[16] != 0 && s > L && s - L < m) {
                    // мёртвый остаток оживёт, только если умрёт сосед: штраф ~ r * P(оба соседа выживут)
                    double hl = hazard(pidLeftOf(gs)), hr = hazard(pidRightOf(gs + s));
                    cst += P[16] * (s - L) * (1 - std::max(hl, hr));
                }
            }
            if (P[51] != 0 && st < 0 && s > L) cst -= P[51] * std::max(rateB(pidLeftOf(gs)), rateB(pidRightOf(gs + s)));
            if (P[52] != 0 && st < 0 && s == L) cst += P[52] * std::max(rateB(pidLeftOf(gs)), rateB(pidRightOf(gs + s)));
            if (P[29] != 0 && st >= 0 && (kin(pidLeftOf(gs), st) || kin(pidRightOf(gs + s), st))) cst -= P[29];
            if (P[40 - 1 + 0] != 0 && false) {}
            if (P[15] > 0 && P[33] == 2 && st < 0 && s > L && s - L < m && s - L >= P[15]) cst = 1e17;   // мёртвый остаток — не кандидат
            if (P[43] != 0 && s != L && s <= MAXL && cntTot > 0) cst += P[43] * cnt[s] / cntTot;   // не ломать дыру «ходового» размера
            if (P[15] > 0 && P[33] == 3 && st < 0 && s > L && s - L < m && s - L >= P[15]) cst += 1e6;   // мёртвый остаток — в последнюю очередь
            if (tieAge != 0) {
                // разбивка ничьих по возрасту соседей (молодые = большой PID)
                int lp = pidLeftOf(gs), rp = pidRightOf(gs + s);
                double young = std::max(lp, rp), old = std::min(lp, rp);
                if (s == L) cst += tieAge * 1e-7 * old;        // точное — между старыми
                else cst -= tieAge * 1e-7 * young;             // остаток — к молодому
            }
            if (endgame && st < 0) cst = s - L;   // P53: эндшпиль — чистый best fit по дырам
            if (P[72] > 0 && st < 0 && !(P[13] > 0 && T < P[13])) {
                // P72: «точно, иначе остаток ≥ m самого редкого среди дыр размера, иначе хвост (с докопкой)»
                int r = s - L;
                if (r == 0) cst = -1e9;
                else if (r < m) cst = 1e17;
                else cst = (r > MAXL ? P[73] : (double)bucketSize(r)) * 1000 + s;
            }
            if (bToFront && P[55] != 0) {   // P55: подтверждённый B — только в дыры, где оба соседа из B (или край)
                int lp = pidLeftOf(gs), rp = pidRightOf(gs + s);
                if ((lp > 0 && stream[lp] != 1) || (rp > 0 && stream[rp] != 1)) cst = 1e17;
            }
            if (verifiedB && st == 0 && P[55] == 2) {   // P55=2: A не занимает дыры между двумя B — они для B
                int lp = pidLeftOf(gs), rp = pidRightOf(gs + s);
                if (lp > 0 && rp > 0 && stream[lp] == 1 && stream[rp] == 1) cst = 1e17;
            }
            if (cst < bestC) { bestC = cst; bgs = gs; bs = s; }
        };
        for (int s = L; s < NB; s++) {
            int k = 0;
            for (auto it = bucket[s].begin(); it != bucket[s].end() && k < M; ++it) {
                if (*it == ts || *it + s == F) continue;  // хвост и граничный промежуток — отдельно
                k++; eval(*it, s);
            }
        }
        {
            int k = 0;
            for (auto it = big.begin(); it != big.end() && k < M; ++it) {
                if (it->second == ts || it->second + it->first == F) continue;
                k++; eval(it->second, it->first);
            }
        }
#ifdef SIM
        if (forcePid == pid) {
            // принудительное отклонение: forceAlt = -1 — хвост, иначе дыра размера forceAlt
            if (forceAlt < 0) bgs = -1;
            else { bgs = -1; for (int g : bucket[forceAlt]) if (g != ts && g + forceAlt != F) { bgs = g; bs = forceAlt; bestC = 0; break; } }
        }
        if (listPid == pid) {
            optBest = bgs >= 0 ? bs : -1; optList.clear();
            for (int s = L; s < NB; s++) for (int g : bucket[s]) { if (g == ts || g + s == F) continue; optList.push_back(s); break; }
            optList.push_back(-1);
        }
        if (epsRand > 0 && (rng() % 1000) < epsRand) {
            // случайный выбор: любая подходящая дыра (по одной на размер) или хвост
            int opts[NB + 2], no = 0;
            for (int s = L; s < NB; s++) for (int g : bucket[s]) { if (g == ts || g + s == F) continue; opts[no++] = s; break; }
            int pick = (int)(rng() % (no + 1));
            if (getenv("SFEPSLOG")) fprintf(stderr, "dev t=%d pid=%d L=%d best=%s%d pick=%s%d opts=%d\n", curCmd, pid, L, bgs >= 0 ? "hole" : "tail", bgs >= 0 ? bs : 0, pick == no ? "tail" : "hole", pick == no ? 0 : opts[pick], no);
            if (pick == no) bgs = -1;
            else { int s = opts[pick]; for (int g : bucket[s]) if (g != ts && g + s != F) { bgs = g; bs = s; break; } bestC = 0; }
        }
#endif
        ll N = 0; int D;
        bool useGap = bgs >= 0 && bestC < P[7] + 1e17 * (P[7] == 0);
        if (st == 1 && P[12] == 2) useGap = false;
        if (bToFront && P[55] == 0) useGap = false;   // прогноз подтвердился: B — к границе роста
        if (useGap && P[14] == 5 && bestC > P[19]) useGap = false;   // дыра не стоит того — в хвост
        // P15: если лучший промежуток оставляет мёртвый остаток >= P15 — лучше в хвост
        if (bgs >= 0 && bestC >= 1e17) useGap = false;
        if (useGap && P[15] > 0 && P[33] != 2 && !(P[33] == 1 && st >= 0)) {
            int r = bs - L;
            if (r > 0 && r < m) {
                double keep = P[25] ? 1 - std::max(hz2(pidLeftOf(bgs)), hz2(pidRightOf(bgs + bs))) : 1;
                if (r * keep >= (P[34] > 0 ? P[34] * m : P[15])) useGap = false;
            }
        }
        if (useGap && P[14] == 6 && st < 0 && !endgame && !(P[72] > 0 && !(P[13] > 0 && T < P[13]))) {
            int tl = (int)(T - tailStart());
            if (tl >= L && LW[13] < bestC) useGap = false;   // хвост — только без докопки
            if (LW[24] != 0 && LW[24] < bestC) useGap = false;   // слишком дорогой кандидат — в хвост даже с докопкой
        }
        // P45: лучшая дыра оставит мёртвый остаток (>= P46), а в уже выкопанном хвосте труба помещается — в хвост
        if (useGap && P[45] != 0 && st < 0 && bs > L && bs - L < m && bs - L >= std::max(1.0, P[46]) &&
            T - tailStart() >= L) useGap = false;
        if (useGap) {
            if (P[92] != 0 && bs > L) { int lp = pidLeftOf(bgs), rp = pidRightOf(bgs + bs); if (lp > 0 && rp > 0) raceAdd(lp, rp); }
            D = chooseSide(bgs, bs, L, st);
            takeFrom(bgs, D, L);
        } else {
            if (P[40] == 1 && st == 1 && !streamsVerified()) st = 0;   // не подтверждено — B кладём как A
            // P58: «двусторонний хвост» — серия равных коротких длин (индекс в серии >= P59, L <= P60) до обнаружения
            // потоков: нечётные по индексу — к правой границе роста, чётные — к левой, без лишней докопки
            bool pseudo = P[58] != 0 && st < 0 && runIdx >= P[59] && L <= P[60];
            if (pseudo) st = runIdx & 1;
            bool smallChunk = P[40] == 2 && st >= 0 && !streamsVerified();
            int fgs = (F > 0 && tag[F] < 0) ? -tag[F] - 1 : F;
            int fl = F - fgs;
            int chunk = (int)P[20];
            if (P[26] > 0) chunk = std::min(chunk, std::max(L, (int)(P[26] * bVol)));
            if (P[35] > 0) chunk = std::min(chunk, std::max(2 * L, (int)(P[35] * T)));   // не больше доли траншеи
            if (P[37] > 0 && st == 1) chunk = std::min(chunk, std::max(L, (int)(P[37] * fbVol)));
            if (P[96] > 0) chunk = P[97] > 0 ? std::max(2 * L, (int)(P[96] * L)) : std::min(chunk, std::max(2 * L, (int)(P[96] * L)));   // P96: участок ~ длине трубы
            if (shadow && P[101] > 0) chunk = std::min<ll>(chunk, std::max<ll>(std::max<ll>(std::max<ll>((ll)P[100], (ll)(P[106] * L)), (ll)(P[108] * T)), (ll)(P[101] * (shadow->freeIn() - freeIn()))));   // P99: участок ~ доказанной экономии потоков
            if (st >= 0 && bFalsified()) chunk = 0;
#ifdef SIM
            if (getenv("SFCH") && st >= 0) fprintf(stderr, "chunk t=%d st=%d chunk=%d delB=%.2f delN=%.2f aliveLab=%d aliveTot=%d F=%d T=%lld\n", curCmd, st, chunk, delB, delN, aliveB1, aliveTot, F, T);
#endif
            if (smallChunk) chunk = std::min(chunk, (int)P[42]);
            if (pseudo) chunk = 0;
            bool far = st >= 0 && ((st == 1) != flip);   // этот поток растёт от правого края участка (F) влево; P107 — участки по очереди зеркальны
            if (fl >= L) {
                D = far ? F - L : fgs;
                takeFrom(fgs, D, L);
                if (far) F = D;
            } else if (F == T) {
                int need = L - fl;
                if (far) need = std::max(need, chunk - fl);
                N = (need + 19) / 20;
                dig(N);
                D = far ? (int)T - L : fgs;
                if (far && P[102] > 0) D = fgs + (int)((T - fgs) / L) * L - L;   // P102: участок кратен L — фронты A и B сойдутся без остатка
                takeFrom(fgs, D, L);
                F = far ? D : (int)T;
            } else {
                // граничный промежуток исчерпан — новый участок в конце траншеи
                if (P[107] > 0 && st >= 0) { flip = !flip; far = (st == 1) != flip; }
                ts = tailStart();
                int tl = (int)(T - ts);
                int need = L - tl;
                if (st >= 0 && (far || P[98] == 0)) need = std::max(need, chunk - tl);   // P98: ближний поток не копает участок — его копает дальний
                N = need > 0 ? (need + 19) / 20 : 0;
                dig(N);
                ts = tailStart();
                D = far ? (int)T - L : ts;
                if (far && P[102] > 0) D = ts + (int)((T - ts) / L) * L - L;
                takeFrom(ts, D, L);
                F = far ? D : (int)T;
            }
        }
return finishPlace(pid, L, D, N, st, useGap, bs);
    }
    void rebuildBit() {
        // после расширения массива Фенвика пересобрать по живым
        std::fill(aliveBit.begin(), aliveBit.end(), 0);
        for (int p = 1; p < (int)pos.size() && p < (int)aliveBit.size(); p++) if (alive[p]) bitAdd(p, 1);
    }
    std::vector<char> alive;
    ll remove(int pid) {
        int d = pos[pid], l = len[pid];
#ifdef SIM
        {   // токен удаления: длина и ранг от вершины стека длины
            std::vector<int>& v = abl[l]; int idx = (int)v.size() - 1; while (idx >= 0 && v[idx] != pid) idx--;
            int rank = (int)v.size() - 1 - idx; if (rank > 15) rank = 15;
            if (idx >= 0) v.erase(v.begin() + idx);
            if (!inRoll && gTPv[0]) tpPushD(l, rank, lastVictim ? pid - lastVictim : 0, (int)pos.size() - 1 - pid);
            lastVictim = pid;
        }
#endif
        occ -= l;
        if (shadow) shadow->remove(pid);
        if (aliveB1 > 0) { delB = delB * P[105] + (stream[pid] >= 0); delN = delN * P[105] + 1; }
        if (stream[pid] >= 0) aliveB1--;
        if (aliveTot > 1) {
            double q = (double)(aliveTot - bitSum(pid)) / (aliveTot - 1);   // доля живых моложе удаляемой
            qSum = qSum * P[57] + q; qTot = qTot * P[57] + 1;
        }
        alive[pid] = 0; bitAdd(pid, -1); deathL[l] += 1; bOnDeath(pid);
        if (P[92] != 0) raceOnDeath(pid);
        for (int k = 1; k <= MAXL; k++) recDel[k] *= P[31];
        recDelTot = recDelTot * P[31] + 1; recDel[l] += 1; aliveCntL[l]--; aliveTot--;
        if (P[61] >= 2) for (int m = 2; m <= (int)P[61]; m++) {
            for (int r = 0; r < m; r++) recDelMod[m][r] *= P[31];
            recDelMod[m][pid % m] += 1; aliveMod[m][pid % m]--;
        }
        deathAge[ageBucket((int)pos.size() - 1 - pid)] += 1;
        { bool hit = false; for (int k = 0; k < 8; k++) if (lastOfLen[l][k] == pid) hit = true;
          lifoHit = lifoHit * 0.99 + hit; lifoTot = lifoTot * 0.99 + 1; }
        if (P[80] > 0 && zl[d] == l) slots[l].push_back(d);
        int e = freeRange(d, l);
        if (d == F) F = e;
        unsigned char fv = feat[pid];
        if (stream[pid] >= 0) diedLab[stream[pid]] += 1;
        if (lastWasDel && stream[pid] >= 0 && stream[lastDelPid] >= 0) {
            labSame *= 0.998; labTot *= 0.998;
            labTot += 1; labSame += (stream[pid] == stream[lastDelPid]);
        }
        for (int f = 0; f < NF; f++) {
            int v = (fv >> f) & 1;
            deathV[f][v] += 1;
            if (lastWasDel) {
                coSame[f] *= P[24]; coTot[f] *= P[24];
                coTot[f] += 1; coSame[f] += (v == ((feat[lastDelPid] >> f) & 1));
            }
        }
        lastWasDel = true; lastDelPid = pid;
        return d;
    }
    size_t bucketSize(int s) const { return s < NB ? bucket[s].size() : 0; }
    size_t gapCount() const { size_t c = big.size(); for (int s = 0; s < NB; s++) c += bucket[s].size(); return c; }
    void holeStats(ll& holes, ll& small) const {
        holes = small = 0; int ts = tailStart();
        for (int s = 1; s < NB; s++) for (int g : bucket[s]) if (g != ts) { holes += s; if (s < 8) small += s; }
        for (auto& g : big) if (g.second != ts) holes += g.first;
    }
    ~Alloc() { delete shadow; }
    Alloc() { pos.push_back(-1); len.push_back(0); feat.push_back(0); stream.push_back(-1); alive.push_back(0); runOf.push_back(-1); ensure(1024); }
};

void Alloc::shadowCheck() {
    if (shadowDropped) return;
    if (shadow && P[109] > 0 && T > P[109]) { delete shadow; shadow = nullptr; shadowDropped = true; return; }   // P109: память; участок C/2 в такой траншее ничтожен
    if (!shadow) {
        bool any = false;
        for (int f = 0; f < NF; f++) if (kappa(f) >= P[22]) any = true;
        if (!any || P[10] != 2) return;
        shadow = new Alloc(*this); shadow->shadow = nullptr; shadow->noStreams = true;
    }
#ifdef SIM
    if (getenv("SFSH") && nPlaced % atoi(getenv("SFSH")) == 0)
        fprintf(stderr, "t=%d T=%lld in=%lld F=%d gap=%lld | sh T=%lld in=%lld | off=%d feat=%d\n", curCmd, T, freeIn(), F, (ll)((F < T && F > 0 && tag[F] < 0) ? F + tag[F] + 1 : 0), shadow->T, shadow->freeIn(), (int)streamsOff, curFeat);
#endif
}

#ifdef SIM
// ---------------------------------------------------------------- оффлайн-симулятор
// Читает тест целиком (это только для экспериментов: аллокатор видит команды строго по одной),
// проверяет корректность как интерактор и печатает KPI.
static std::vector<int> readTest(const char* fn, double& w) {
    FILE* f = fopen(fn, "r");
    if (!f) { perror(fn); exit(1); }
    int n; if (fscanf(f, "%lf %d", &w, &n) != 2) exit(1);
    std::vector<int> c(n);
    for (int i = 0; i < n; i++) if (fscanf(f, "%d", &c[i]) != 1) exit(1);
    fclose(f);
    return c;
}
// T в конце прогона (или -1 при ошибке); force/list — для анализа важности решений
static ll runOnce(const std::vector<int>& cmd, int forcePid, int forceAlt, int listPid, Alloc** keep) {
    Alloc* A = new Alloc();
    A->forcePid = forcePid; A->forceAlt = forceAlt; A->listPid = listPid;
    int ci = 0;
    for (int c : cmd) { A->curCmd = ci++; if (c > 0) A->place(c); else A->remove(-c); }
    ll T = A->T;
    if (keep) *keep = A; else delete A;
    return T;
}
int main(int argc, char** argv) {
    if (const char* s = getenv("SF")) {
        // формат: "i=v i=v ..."
        char* p = (char*)s;
        while (*p) {
            char* q; long i = strtol(p, &q, 10); if (q == p || *q != '=') break;
            p = q + 1; double v = strtod(p, &q); if (q == p) break; if (i >= 0 && i < 128) P[i] = v; p = q;
        }
    }
    if (const char* s = getenv("SFW")) {
        char* p = (char*)s; int i = 0;
        while (*p && i < 32) { char* q; double v = strtod(p, &q); if (q == p) break; LW[i++] = v; p = q; while (*p == ',' || *p == ' ') p++; }
    }
    bool quiet = getenv("SFQ") != nullptr;
    if (getenv("SFIMPORT")) {
        double w; std::vector<int> cmd = readTest(argv[1], w);
        ll T0 = runOnce(cmd, -1, 0, -1, nullptr);
        int np = 0; for (int c : cmd) if (c > 0) np++;
        printf("base T=%lld\n", T0);
        for (int p = 1; p <= np; p++) {
            Alloc* A; runOnce(cmd, -1, 0, p, &A);
            std::vector<int> opts = A->optList; int best = A->optBest; int L = A->len[p];
            delete A;
            for (int o : opts) {
                if (o == best) continue;
                ll T = runOnce(cmd, p, o, -1, nullptr);
                if (T != T0) printf("pid=%d L=%d best=%d alt=%d dT=%lld\n", p, L, best, o, T - T0);
            }
        }
        return 0;
    }
    double total = 0, wsum = 0;
    for (int a = 1; a < argc; a++) {
        double w; std::vector<int> cmd = readTest(argv[a], w);
        Alloc A;
        gCmds = &cmd; if (getenv("ROLLH")) { A.rollH = atoi(getenv("ROLLH")); A.rollK = getenv("ROLLK") ? atoi(getenv("ROLLK")) : 4; A.rollMargin = getenv("ROLLM") ? atof(getenv("ROLLM")) : 0; A.rollMode = getenv("ROLLMODE") ? atoi(getenv("ROLLMODE")) : 1; A.rollMin = getenv("ROLLMIN") ? atoi(getenv("ROLLMIN")) : 1; }
        for (int v = 0; v < 3; v++) { if (gTPv[v]) delete gTPv[v]; gTPv[v] = new TokPred(); gAcc[v] = 0; } if (getenv("ACCMIN")) A.accMin = atof(getenv("ACCMIN")); if (getenv("HACCMIN")) A.haccMin = atof(getenv("HACCMIN")); gHacc = 0; gPend.clear(); gPendIdx = 0; if (gTP) delete gTP; gTP = new TokPred(); for (int v = 0; v < 3; v++) { if (getenv("PSHARE")) gTPv[v]->minShare = atof(getenv("PSHARE")); if (getenv("PCNT")) gTPv[v]->minCnt = atoi(getenv("PCNT")); }
        std::vector<int> death(1, -1);
        { int p = 0; for (size_t i = 0; i < cmd.size(); i++) { if (cmd[i] > 0) { death.push_back(-1); p++; } else death[-cmd[i]] = (int)i; } }
        A.oracleDeath = death.data();
        if (getenv("SFTIE")) A.tieAge = atof(getenv("SFTIE"));
        if (getenv("SFFORCE")) { A.forcePid = atoi(getenv("SFFORCE")); A.forceAlt = atoi(getenv("SFALT")); }
        if (getenv("SFLIST")) A.listPid = atoi(getenv("SFLIST"));
        if (getenv("SFEPS")) { A.epsRand = atoi(getenv("SFEPS")); A.rs ^= (unsigned long long)atoll(getenv("SFSEED")) * 0x9E3779B97F4A7C15ULL + 1; }
        std::vector<int> cells;  // pid или 0
        std::vector<int> lastLen, freedAt;
        ll occ = 0; bool ok = true;
        long long step = getenv("SFTRACE") ? atoll(getenv("SFTRACE")) : 0; long long ci = 0;
        int ncut = getenv("SFCUTS") ? atoi(getenv("SFCUTS")) : 0; double cutSum = 0; int cutCnt = 0;   // средний KPI по концам в последних 10 %
        for (int c : cmd) {
            if (ncut > 1 && ci > 0) { long long n = (long long)cmd.size(), lo = n - n / 10; long long st = std::max(1LL, n / 10 / ncut); if (ci >= lo && (ci - lo) % st == 0 && !cells.empty()) { cutSum += (double)occ / cells.size(); cutCnt++; } }
            if (step && ci % step == 0) {
                // гистограмма дыр
                long long holes, small; A.holeStats(holes, small);
                printf("  t=%lld T=%lld occ=%lld waste=%lld holes=%lld small<8=%lld\n", ci, A.T, occ, A.T - occ, holes, small);
            }
            A.curCmd = (int)ci;
            ci++;
            if (c > 0) {
                auto r = A.place(c);
                ll T = (ll)cells.size() + 20 * r.first;
                if (r.first < 0 || T != A.T || r.second < 0 || r.second + c > T) { ok = false; break; }
                cells.resize(T, 0); lastLen.resize(T, 0); freedAt.resize(T, -1);
                for (int i = 0; i < c; i++) if (cells[r.second + i]) { ok = false; break; }
                if (!ok) break;
                int pid = (int)A.pos.size() - 1;
                if (getenv("SFDIG") && r.first > 0) {
                    ll holes, small; A.holeStats(holes, small);
                    int mx = 0; for (int s = NB - 1; s > 0; s--) if (!A.bucket[s].empty()) { for (int g : A.bucket[s]) if (g != A.tailStart() && g + s != A.F) { mx = s; break; } if (mx) break; }
                    printf("  dig t=%lld L=%d N=%lld T=%lld occ=%lld holes=%lld maxhole=%d\n", ci, c, r.first, A.T, occ, holes, mx);
                }
                for (int i = 0; i < c; i++) cells[r.second + i] = pid, lastLen[r.second + i] = c;
                occ += c;
            } else {
                int pid = -c;
                ll d = A.remove(pid);
                if (d < 0 || d >= (ll)cells.size() || cells[d] != pid) { ok = false; break; }
                int l = A.len[pid];
                for (int i = 0; i < l; i++) cells[d + i] = 0, freedAt[d + i] = (int)ci;
                occ -= l;
            }
        }
        if (getenv("SFORIG") && ok) {
            // происхождение свободных клеток в конце: 0 — не занималась, иначе длина последней трубы на этом месте
            long long byLen[41] = {0}, never = 0;
            for (size_t i = 0; i < cells.size(); i++) if (!cells[i]) { if (lastLen[i]) byLen[lastLen[i]]++; else never++; }
            printf("  free: never=%lld", never);
            long long sh = 0, lg = 0; for (int l = 1; l <= 40; l++) (l <= 7 ? sh : lg) += byLen[l];
            printf(" was_short=%lld was_long=%lld\n", sh, lg);
        }
        if (getenv("SFFREED") && ok) {
            int from = atoi(getenv("SFFREED")); long long n = (long long)cmd.size();
            long long b[8] = {0}; // освобождено: <100, <300, <1000, <3000, <10000, <30000, раньше, никогда
            int lim[6] = {100, 300, 1000, 3000, 10000, 30000};
            for (size_t i = from; i < cells.size(); i++) if (!cells[i]) {
                if (freedAt[i] < 0) { b[7]++; continue; }
                long long age = n - freedAt[i]; int k = 0; while (k < 6 && age >= lim[k]) k++; b[k]++;
            }
            printf("  freed-age from %d: <100:%lld <300:%lld <1k:%lld <3k:%lld <10k:%lld <30k:%lld older:%lld never:%lld\n", from, b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]);
        }
        if (getenv("SFMAP") && ok) {
            // карта траншеи: длины труб (с возрастом) и дыры [h]
            size_t i = 0; int np = (int)A.pos.size() - 1;
            while (i < cells.size()) {
                if (cells[i]) { int p = cells[i]; printf("%d", A.len[p]); if (getenv("SFMAPAGE")) printf("(%d)", np - p); printf(" "); i += A.len[p]; }
                else { size_t j = i; while (j < cells.size() && !cells[j]) j++; printf("[%zu] ", j - i); i = j; }
            }
            printf("\n");
        }
        if (getenv("SFHOLES") && ok) {
            int from = atoi(getenv("SFHOLES"));
            long long h[200] = {0}; long long tot = 0;
            for (size_t i = from; i < cells.size();) { if (cells[i]) { i++; continue; } size_t j = i; while (j < cells.size() && !cells[j]) j++; int s = (int)std::min<size_t>(j - i, 199); if (j < cells.size()) { h[s] += j - i; tot += j - i; } i = j; }
            printf("  holes from %d: total %lld; by size (meters):", from, tot);
            for (int s = 1; s < 200; s++) if (h[s]) printf(" %d:%lld", s, h[s]);
            printf("\n");
        }
        double kpi = ok && !cells.empty() ? (double)occ / cells.size() : 0;
        if (ncut > 1 && ok && !cells.empty()) kpi = (cutSum + kpi) / (cutCnt + 1);
        if (!ok) fprintf(stderr, "WA on %s\n", argv[a]);
        total += kpi * w; wsum += w;
        if (!quiet) printf("%-40s w=%-4g kpi=%.5f T=%zu occ=%lld\n", argv[a], w, kpi, cells.size(), occ);
    }
    printf("TOTAL %.5f / %.4g  avg %.5f\n", total, wsum, total / wsum);
}
#else
// ---------------------------------------------------------------- интерактивный режим
static char ibuf[1 << 16];
static int ipos = 0, ilen = 0;
static inline int gc() {
    if (ipos == ilen) { ilen = (int)read(0, ibuf, sizeof ibuf); ipos = 0; if (ilen <= 0) return -1; }
    return ibuf[ipos++];
}
static bool readInt(int& x) {
    int c = gc();
    while (c != '-' && (c < '0' || c > '9')) { if (c == -1) return false; c = gc(); }
    bool neg = false; if (c == '-') { neg = true; c = gc(); }
    x = 0;
    while (c >= '0' && c <= '9') { x = x * 10 + (c - '0'); c = gc(); }
    if (neg) x = -x;
    return true;
}
static void writeOut(const char* s, int n) {
    while (n > 0) { ssize_t k = write(1, s, n); if (k <= 0) return; s += k; n -= (int)k; }
}
static int putNum(char* b, ll v) {
    char t[24]; int n = 0;
    if (v == 0) t[n++] = '0';
    while (v > 0) { t[n++] = char('0' + v % 10); v /= 10; }
    for (int i = 0; i < n; i++) b[i] = t[n - 1 - i];
    return n;
}
int main() {
    if (const char* s = getenv("SF")) {   // те же ручки, что в SIM: "i=v i=v ..."
        char* p = (char*)s;
        while (*p) {
            char* q; long i = strtol(p, &q, 10); if (q == p || *q != '=') break;
            p = q + 1; double v = strtod(p, &q); if (q == p) break; if (i >= 0 && i < 128) P[i] = v; p = q;
        }
    }
    if (const char* s = getenv("SFW")) {
        char* p = (char*)s; int i = 0;
        while (*p && i < 32) { char* q; double v = strtod(p, &q); if (q == p) break; LW[i++] = v; p = q; while (*p == ',' || *p == ' ') p++; }
    }
    Alloc A;
    int c; char ob[64];
    while (readInt(c)) {
        if (c == 0) return 0;
        int n;
        if (c > 0) {
            auto r = A.place(c);
            n = putNum(ob, r.first); ob[n++] = ' '; n += putNum(ob + n, r.second);
        } else {
            n = putNum(ob, A.remove(-c));
        }
        ob[n++] = '\n';
        writeOut(ob, n);
    }
    return 0;
}
#endif
