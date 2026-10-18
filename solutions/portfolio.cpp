// Вариант B: портфель политик размещения + выбор политики по истории (теневые симуляции).
// Идея: набор «чистых» политик (best/first/last/worst/exact fit, арены для скоротечных труб по предсказателю времени жизни)
// прогоняется параллельно в теневых копиях на всей истории команд; реальная траншея ведётся текущим лидером.
// Сборка решения:  c++ -std=c++17 -O2 -o build/portfolio solutions/portfolio.cpp
// Стенд (внутрипроцессный прогон тестов): scripts/sim.cpp, сборка scripts/build_sim.bat; эксперименты: scripts/exp_*.cpp
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
#ifdef _MSC_VER
#include <intrin.h>
#endif

static inline int ctz64(uint64_t x) {
#ifdef _MSC_VER
    unsigned long r;
    _BitScanForward64(&r, x);
    return (int)r;
#else
    return __builtin_ctzll(x);
#endif
}
static inline int clz64idx(uint64_t x) {  // индекс старшего бита
#ifdef _MSC_VER
    unsigned long r;
    _BitScanReverse64(&r, x);
    return (int)r;
#else
    return 63 - __builtin_clzll(x);
#endif
}

// ---------------------------------------------------------------- политики

enum Sel { BEST, FIRST, LAST, WORST, EXACTWORST };

struct Pol {
    const char* name;
    Sel sel;          // как выбирать дыру для обычных (P) труб
    bool tieHigh;     // при равенстве размеров брать дыру с большим адресом
    bool right;       // класть к правому краю дыры
    bool tailFirst;   // хвост — обычный кандидат: если влезает, класть туда до поиска дыр
    int js;           // горизонт предсказателя для класса T (скоротечные); -1 — без арен, все трубы класса P
    int Cdiv, Cmin, Cmax;  // размер чанка арены: clamp(trench / Cdiv, Cmin, Cmax)
    int Tsmall;       // P-трубы могут занимать T-дыры размером <= Tsmall
};

static const Pol POLS[] = {
    {"bf-L", BEST, false, false, false, -1, 0, 0, 0, 0},
    {"bf-R", BEST, false, true, false, -1, 0, 0, 0, 0},
    {"bf-hiR", BEST, true, true, false, -1, 0, 0, 0, 0},
    {"ff-L", FIRST, false, false, false, -1, 0, 0, 0, 0},
    {"lf-L", LAST, false, false, false, -1, 0, 0, 0, 0},
    {"wf-L", WORST, false, false, false, -1, 0, 0, 0, 0},
    {"xw-L", EXACTWORST, false, false, false, -1, 0, 0, 0, 0},
    {"a2", BEST, false, false, false, 2, 100, 40, 400, 40},
    {"a3", BEST, false, false, false, 3, 100, 40, 400, 40},
    {"a4", BEST, false, false, false, 4, 100, 40, 400, 40},
    {"a5", BEST, false, false, false, 5, 100, 40, 400, 40},
    {"a4x", EXACTWORST, false, false, false, 4, 100, 40, 400, 40},
};
static const int NPOL = sizeof(POLS) / sizeof(POLS[0]);

// ---------------------------------------------------------------- множество свободных дыр одного тега

struct HoleSet {
    std::set<int> sz[40];                 // дыры длины 1..39: адреса
    std::set<std::pair<int, int>> big;    // дыры длины >= 40: (длина, адрес)
    std::set<int> bigAddr;                // они же по адресу
    uint64_t mask = 0;                    // бит s (1..39) — есть дыра длины s; бит 40 — есть большая

    void add(int a, int h) {
        if (h < 40) { sz[h].insert(a); mask |= 1ull << h; }
        else { big.insert({h, a}); bigAddr.insert(a); mask |= 1ull << 40; }
    }
    void del(int a, int h) {
        if (h < 40) { sz[h].erase(a); if (sz[h].empty()) mask &= ~(1ull << h); }
        else { big.erase({h, a}); bigAddr.erase(a); if (big.empty()) mask &= ~(1ull << 40); }
    }
    // кандидат по правилу; maxSize ограничивает размер дыры (для P-труб в T-дырах). false — нет подходящей.
    // h == -1 в ответе: большая дыра (>= 40), длину уточняет вызывающий.
    bool pick(Sel sel, bool tieHigh, int L, int maxSize, int& a, int& h) const {
        uint64_t mm = mask & (~0ull << L);
        if (maxSize < 40) mm &= (2ull << maxSize) - 1;
        else if ((mm >> 40 & 1) && maxSize < (1 << 30) && big.begin()->first > maxSize) mm &= ~(1ull << 40);
        if (!mm) return false;
        switch (sel) {
        case BEST: {
            int s = ctz64(mm);
            if (s < 40) { a = tieHigh ? *sz[s].rbegin() : *sz[s].begin(); h = s; }
            else {
                auto it = big.begin();
                if (tieHigh) { it = big.upper_bound({it->first, 0x7fffffff}); --it; }
                a = it->second; h = it->first;
            }
            return true;
        }
        case FIRST:
        case LAST: {
            bool first = sel == FIRST;
            int best = first ? 0x7fffffff : -1, bh = 0;
            for (uint64_t m = mm & ((1ull << 40) - 1); m; m &= m - 1) {
                int s = ctz64(m);
                int x = first ? *sz[s].begin() : *sz[s].rbegin();
                if (first ? x < best : x > best) { best = x; bh = s; }
            }
            if (mm >> 40 & 1) {
                int x = first ? *bigAddr.begin() : *bigAddr.rbegin();
                if (first ? x < best : x > best) { best = x; bh = -1; }
            }
            a = best; h = bh;
            return true;
        }
        case WORST:
            if (mm >> 40 & 1) { a = big.rbegin()->second; h = big.rbegin()->first; }
            else { int s = clz64idx(mm); a = *sz[s].begin(); h = s; }
            return true;
        case EXACTWORST:
            if (L < 40 && (mask >> L & 1) && L <= maxSize) { a = *sz[L].begin(); h = L; return true; }
            if (mm >> 40 & 1) { a = big.rbegin()->second; h = big.rbegin()->first; }
            else { int s = clz64idx(mm); a = *sz[s].begin(); h = s; }
            return true;
        }
        return false;
    }
};

// ---------------------------------------------------------------- состояние траншеи

struct Alloc {
    int trench = 0;      // длина траншеи
    int tailStart = 0;   // начало свободного хвоста (хвост = [tailStart, trench))
    struct HV { int len; int tag; };
    std::map<int, HV> byAddr;             // свободные дыры внутри траншеи (не касаются конца): начало -> (длина, тег)
    HoleSet hs[2];                        // 0 — обычные (U), 1 — арены/дыры скоротечных (T)
    std::vector<int> pos, len;            // по PID (с 1)
    std::vector<uint8_t> cls;             // класс трубы: 0 = P (обычная), 1 = T (скоротечная)
    long long live = 0;

    Alloc() { pos.assign(1, -1); len.assign(1, 0); cls.assign(1, 0); }

    void addHole(int a, int h, int tag) { byAddr[a] = {h, tag}; hs[tag].add(a, h); }
    void delHole(int a) {
        auto it = byAddr.find(a);
        hs[it->second.tag].del(a, it->second.len);
        byAddr.erase(it);
    }
    int holeLen(int a) const { return byAddr.find(a)->second.len; }

    // вырезать [D, D+L) из дыры (a, h, tag)
    void carve(int a, int h, int tag, int D, int L) {
        delHole(a);
        if (D > a) addHole(a, D - a, tag);
        if (D + L < a + h) addHole(D + L, a + h - (D + L), tag);
    }
    void commit(int D, int L, int c) { pos.push_back(D); len.push_back(L); cls.push_back((uint8_t)c); live += L; }

    // уложить трубу длины L класса c; возвращает N (докопка в 20-метровых сегментах) и D (позиция)
    void place(const Pol& p, int L, int c, int& N, int& D) {
        N = 0;
        int a, h;
        if (c == 1) {  // скоротечная: best fit по T-дырам, иначе новый чанк в хвосте
            if (hs[1].pick(BEST, false, L, 1 << 30, a, h)) {
                if (h < 0) h = holeLen(a);
                D = a; carve(a, h, 1, D, L); commit(D, L, 1); return;
            }
            int C = std::min(std::max(trench / p.Cdiv, p.Cmin), p.Cmax);
            int need = std::max(C, L);
            int tailLen = trench - tailStart;
            if (tailLen < need) { N = (need - tailLen + 19) / 20; trench += 20 * N; }
            D = tailStart;
            tailStart += need;
            if (need > L) addHole(D + L, need - L, 1);
            commit(D, L, 1);
            return;
        }
        int tailLen = trench - tailStart;
        if (p.tailFirst && tailLen >= L) {
            D = tailStart; tailStart += L; commit(D, L, 0); return;
        }
        // P: кандидаты из U-дыр и малых T-дыр
        int a0 = 0, h0 = 0, a1 = 0, h1 = 0;
        bool f0 = hs[0].pick(p.sel, p.tieHigh, L, 1 << 30, a0, h0);
        bool f1 = p.js >= 0 && hs[1].pick(p.sel, p.tieHigh, L, p.Tsmall, a1, h1);
        if (f0 && h0 < 0) h0 = holeLen(a0);
        if (f1 && h1 < 0) h1 = holeLen(a1);
        if (f0 || f1) {
            bool useT;
            if (!f1) useT = false;
            else if (!f0) useT = true;
            else {
                switch (p.sel) {
                case BEST: useT = h1 < h0 || (h1 == h0 && (p.tieHigh ? a1 > a0 : a1 < a0)); break;
                case FIRST: useT = a1 < a0; break;
                case LAST: useT = a1 > a0; break;
                default: useT = h1 > h0; break;
                }
            }
            int tag = useT ? 1 : 0;
            a = useT ? a1 : a0; h = useT ? h1 : h0;
            if (p.right) D = a + h - L; else D = a;
            carve(a, h, tag, D, L);
            commit(D, L, 0);
            return;
        }
        if (tailLen < L) { N = (L - tailLen + 19) / 20; trench += 20 * N; }
        D = tailStart;
        tailStart += L;
        commit(D, L, 0);
    }

    int remove(int pid) {
        int p = pos[pid], l = len[pid];
        int tag = cls[pid];
        pos[pid] = -1;
        live -= l;
        int s = p, e = p + l;
        auto it = byAddr.lower_bound(p);
        if (it != byAddr.begin()) {
            auto jt = std::prev(it);
            if (jt->first + jt->second.len == p) { s = jt->first; tag = std::max(tag, jt->second.tag); delHole(jt->first); }
        }
        if (e == tailStart) {
            tailStart = s;
        } else {
            auto r = byAddr.find(e);
            if (r != byAddr.end()) { tag = std::max(tag, r->second.tag); e += r->second.len; delHole(r->first); }
            addHole(s, e - s, tag);
        }
        return p;
    }
};

// ---------------------------------------------------------------- предсказатель времени жизни
// Для каждой трубы при укладке строятся ключи контекстов (семейства признаков). По завершённым наблюдениям
// («умерла за H команд» / «дожила до возраста H») каждое семейство ведёт таблицу P(смерть за H | ключ) по набору горизонтов.
// Качество семейства оценивается онлайн (префиксная кросс-энтропия с затуханием); для прогноза берётся лучшее.
struct Predictor {
    static const int NH = 7;                  // горизонты H_j = 8 * 4^j команд
    static const int NF = 15;                 // семейства признаков
    static int H(int j) { return 8 << (2 * j); }
    struct Cell { float a = 0, n = 0; };
    struct Entry { Cell c[NH]; };
    std::unordered_map<uint64_t, Entry> tab[NF];
    float loss[NF][NH];
    // по трубам (pid с 1)
    std::vector<int> birth;
    std::vector<uint8_t> resolved;            // маска разрешённых горизонтов
    std::vector<uint8_t> alive;
    std::vector<Entry*> ent;                  // NF указателей на трубу, подряд
    std::vector<int> queue[NH];               // очередь pid по рождению (для «дожила до H»)
    size_t qhead[NH];
    int t = 0;                                // время = число обработанных команд
    int lp1 = 0, lp2 = 0, rp = 0;             // контекст: две прошлые длины, позиция в серии одинаковых длин

    Predictor() {
        for (int f = 0; f < NF; f++) for (int j = 0; j < NH; j++) loss[f][j] = 0.0f;
        for (int j = 0; j < NH; j++) qhead[j] = 0;
        birth.assign(1, 0); resolved.assign(1, 0); alive.assign(1, 0); ent.assign(NF, nullptr);
    }

    static inline uint64_t mix(uint64_t a, uint64_t b) { return a * 1000003ull + b + 12345; }

    // вызывается при укладке очередной трубы (до размещения); возвращает вероятности p[j] смерти за H_j
    void onPlace(int L, float* p) {
        int pid = (int)birth.size();
        uint64_t keys[NF];
        int rpn = (L == lp1) ? rp + 1 : 0;
        keys[0] = 0;
        keys[1] = L;
        keys[2] = rpn % 2;
        keys[3] = mix(L, rpn % 2);
        keys[4] = rpn % 3;
        keys[5] = mix(L, rpn % 3);
        for (int m = 2; m <= 8; m++) keys[4 + m] = (uint64_t)(pid % m);  // 6..12
        keys[13] = mix(lp1, L);
        keys[14] = mix(mix(lp2, lp1), L);
        birth.push_back(t);
        resolved.push_back(0);
        alive.push_back(1);
        size_t base = ent.size();
        ent.resize(base + NF);
        for (int f = 0; f < NF; f++) ent[base + f] = &tab[f][keys[f]];
        for (int j = 0; j < NH; j++) {
            int bf = 0;
            float bl = 1e30f;
            for (int f = 0; f < NF; f++) {
                if (loss[f][j] < bl) { bl = loss[f][j]; bf = f; }
            }
            const Cell& c = ent[base + bf]->c[j];
            p[j] = (c.a + 0.5f) / (c.n + 1.0f);
        }
        for (int j = 0; j < NH; j++) queue[j].push_back(pid);
        rp = rpn; lp2 = lp1; lp1 = L;
    }

    void resolve(int pid, int j, int y) {
        resolved[pid] |= 1 << j;
        Entry** e = &ent[(size_t)pid * NF];
        for (int f = 0; f < NF; f++) {
            Cell& c = e[f]->c[j];
            float p = (c.a + 0.5f) / (c.n + 1.0f);
            float l = -std::log(y ? p : 1.0f - p);
            loss[f][j] = loss[f][j] * 0.995f + l;
            c.n += 1; c.a += y;
        }
    }

    void onDelete(int pid) {
        alive[pid] = 0;
        int age = t - birth[pid];
        for (int j = 0; j < NH; j++)
            if (!(resolved[pid] >> j & 1) && age <= H(j)) resolve(pid, j, 1);
    }

    // конец команды: продвинуть время, разрешить «дожила до H»
    void tick() {
        t++;
        for (int j = 0; j < NH; j++) {
            while (qhead[j] < queue[j].size()) {
                int pid = queue[j][qhead[j]];
                if (birth[pid] + H(j) > t) break;
                qhead[j]++;
                if (!(resolved[pid] >> j & 1) && alive[pid]) resolve(pid, j, 0);
            }
        }
    }
};

// ---------------------------------------------------------------- движок (портфель)

struct Engine {
    Alloc real;
    std::vector<Alloc> sh;   // теневые симуляции чистых политик (пусто в режиме одной политики)
    Predictor pred;
    int fixed;               // >= 0: одна чистая политика; -1: портфель
    int cur = 0;
    double hyst = 0.0;       // гистерезис: переключаться, только если лидер короче текущей на hyst*trench

    explicit Engine(int p = -1) : fixed(p) {
        if (p < 0) sh.resize(NPOL);
    }
    static int classOf(const Pol& p, const float* pr) { return p.js >= 0 && pr[p.js] >= 0.5f ? 1 : 0; }
    void place(int L, int& N, int& D) {
        float pr[Predictor::NH];
        pred.onPlace(L, pr);
        if (fixed >= 0) {
            real.place(POLS[fixed], L, classOf(POLS[fixed], pr), N, D);
        } else {
            for (int k = 0; k < NPOL; k++) { int n, d; sh[k].place(POLS[k], L, classOf(POLS[k], pr), n, d); }
            int best = cur;
            for (int k = 0; k < NPOL; k++) if (sh[k].trench < sh[best].trench) best = k;
            if (best != cur && sh[best].trench < sh[cur].trench - hyst * sh[cur].trench) cur = best;
            real.place(POLS[cur], L, classOf(POLS[cur], pr), N, D);
        }
        pred.tick();
    }
    int remove(int pid) {
        pred.onDelete(pid);
        pred.tick();
        for (auto& a : sh) a.remove(pid);
        return real.remove(pid);
    }
};

#ifndef SIM
// ---------------------------------------------------------------- ввод-вывод (боевой режим)
int main() {
    Engine eng(-1);
    int c;
    while (scanf("%d", &c) == 1) {
        if (c == 0) return 0;
        if (c > 0) {
            int N, D;
            eng.place(c, N, D);
            printf("%d %d\n", N, D);
        } else {
            printf("%d\n", eng.remove(-c));
        }
        fflush(stdout);
    }
}
#endif
