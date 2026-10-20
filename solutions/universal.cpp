// universal.cpp — онлайн-упаковка труб в траншею («Gas pipeline construction», NSUts).
//
// ПОЧЕМУ ТАК. KPI = (занятые клетки в конце) / (длина траншеи в конце). Живые трубы в конце заданы входом,
// траншея не укорачивается => максимизировать KPI = МИНИМИЗИРОВАТЬ ЧИСЛО ДОКОПОК по 20 м. Это online
// dynamic storage allocation с монотонным ростом кучи; качество = число «лишних» блоков
// (T_final − 20·ceil(peak_live/20)).
//
// КОМПОНЕНТЫ.
//   GapSet — свободные промежутки: set (len,start) для best fit + map начало→длина и конец→начало для
//            слияния за O(1); промежутки арен помечены владельцем и не видны best fit.
//   Арены  — если подряд кладут одну и ту же длину L (серия без удалений), труба ранга r идёт в поток
//            r mod 2, каждый поток растёт в своих блоках фиксированной ёмкости. В тестах «серия одинаковых
//            труб, затем удаляют каждую вторую» жертвы оказываются в одном потоке, и их смерть освобождает
//            ОДИН непрерывный кусок, пригодный для длинных труб. Неиспользованный остаток блока при конце
//            серии возвращается в общий фонд.
//   Engine — стратегия: best fit, сторона остатка по «молодости» соседа, ровно столько докопок, сколько нужно.
//
// Ручки экспериментов: при сборке с -DSIM массив K[] задаётся снаружи (стенд scripts/sim2.cpp, env K).
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

typedef long long ll;

static const ll MAX_T = 100000000;

// Значения по умолчанию:
//   K[0]=RUN_MIN (серия для арен), K[1]=резерв, K[2]=включить арены, K[3]=политика малой траншеи,
//   K[4]=порог малой траншеи (клетки), K[5]=ёмкость блока арены (клетки), K[6]=макс. слотов в блоке,
//   K[7]=гейт арен по направлению длин (0=всегда, 1=только растущие серии, 2=растущие или первая).
static const double KD[64] = { 8, 0, 1, 2, 100, 96, 16, 1 };
#ifdef SIM
extern double K[64];
#else
static double K[64];
#endif

// ------------------------------------------------------------------ свободные промежутки
struct GapSet {
    std::set<std::pair<ll, ll>> byLen;          // (len, start)
    std::unordered_map<ll, ll> lenAt;           // start -> len
    std::unordered_map<ll, ll> startAtEnd;      // end -> start
    std::unordered_map<ll, int> owner;          // start -> тег блока арены (0 = свободный)

    int getLen(ll s) const {
        auto it = lenAt.find(s);
        return it == lenAt.end() ? -1 : (int)it->second;
    }
    int getOwner(ll s) const {
        auto it = owner.find(s);
        return it == owner.end() ? 0 : it->second;
    }

    // вставка; приклеиваем только соседей с тем же владельцем (теги блоков уникальны => блоки не слипаются)
    void addGap(ll s, ll len, int own = 0) {
        if (len <= 0) return;
        auto it = startAtEnd.find(s);
        if (it != startAtEnd.end()) {
            ll ls = it->second;
            if (getOwner(ls) == own) {
                ll l2 = lenAt[ls];
                delGap(ls);
                s = ls;
                len += l2;
            }
        }
        int rl = getLen(s + len);
        if (rl >= 0 && getOwner(s + len) == own) {
            ll rs = s + len;
            delGap(rs);
            len += rl;
        }
        lenAt[s] = len;
        startAtEnd[s + len] = s;
        byLen.insert({ len, s });
        if (own) owner[s] = own;
    }

    void delGap(ll s) {
        auto it = lenAt.find(s);
        if (it == lenAt.end()) return;
        startAtEnd.erase(s + it->second);
        byLen.erase({ it->second, s });
        owner.erase(s);
        lenAt.erase(it);
    }

    // вырезать кусок [p, p+n) из промежутка [s, s+gl); остатки сохраняют владельца
    void cut(ll s, ll gl, ll p, ll n, int own) {
        delGap(s);
        if (p > s) addGap(s, p - s, own);
        if (s + gl > p + n) addGap(p + n, s + gl - (p + n), own);
    }

    ll bestFit(int L) {
        auto it = byLen.lower_bound({ L, -1 });
        while (it != byLen.end()) {
            ll s = it->second;
            int l = getLen(s);
            if (l != (int)it->first) { it = byLen.erase(it); continue; }
            if (getOwner(s) == 0) return s;
            ++it;
        }
        return -1;
    }

    // точное попадание (промежуток ровно длины L)
    ll exactFit(int L) {
        auto it = byLen.lower_bound({ L, -1 });
        while (it != byLen.end() && it->first == L) {
            ll s = it->second;
            if (getLen(s) != L) { it = byLen.erase(it); continue; }
            if (getOwner(s) == 0) return s;
            ++it;
        }
        return -1;
    }

    // самый большой промежуток длины >= L (worst fit)
    ll worstFit(int L) {
        for (auto it = byLen.rbegin(); it != byLen.rend(); ++it) {
            if (it->first < L) break;
            ll s = it->second;
            if (getLen(s) != (int)it->first) continue;
            if (getOwner(s) == 0) return s;
        }
        return -1;
    }
};

// ------------------------------------------------------------------ движок
struct Engine {
    static const int MAXL = 40;
    static const int NCLASS = 2;
    static const int NARENA = (MAXL + 1) * NCLASS;

    struct Arena {
        ll start = -1, len = 0, used = 0;   // блок: начало, ёмкость, занято трубами
        int slots = 0, usedSlots = 0, blocks = 0;
        int tag = 0;
    };

    GapSet g;
    std::map<ll, ll> pipeAt;
    std::vector<ll> pos;
    std::vector<int> plen;
    std::vector<int> parena;                 // pid -> id арены (L*2+c) + 1, 0 = обычная труба
    ll T = 0;

    int runLen = 0, runL = 0;
    int prevRunL = 0;                        // длина предыдущей завершённой серии
    bool ascRun = false;                     // текущая серия длиннее предыдущей (растущие группы)
    bool arenaOn[MAXL + 1] = {};
    Arena ar[MAXL + 1][NCLASS];
    int nextTag = 1;
    ll dugBlocks = 0;                        // докопки, сделанные открытием блока

    Engine() {
        pos.push_back(-1);
        plen.push_back(0);
        parena.push_back(0);
    }

    int pidAt(ll cell) const {
        if (cell < 0 || cell >= T) return 0;
        auto it = pipeAt.upper_bound(cell);
        if (it == pipeAt.begin()) return 0;
        --it;
        ll s = it->first;
        return (cell < s + plen[it->second]) ? (int)it->second : 0;
    }

    static int arenaId(int L, int c) { return L * NCLASS + c; }

    // ------------------------------------------------------------------ укладка
    void place(int L, ll& N, ll& D) {
        if (L != runL) {
            if (runL >= 1 && runL <= MAXL && arenaOn[runL]) closeArenasOf(runL);
            if (runL >= 1) prevRunL = runL;      // длина предыдущей серии (runLen мог обнулиться удалением)
            runL = L;
            runLen = 0;
            ascRun = prevRunL > 0 && L > prevRunL;
        }
        runLen++;
        ll dug = 0;                              // все докопки этого ответа

        // Арены нужны там, где освободившиеся дыры НЕ годятся для будущих труб (группы длин растут),
        // и вредны там, где годятся (группы убывают: best fit и так переиспользует дыры).
        bool arenaGate = K[7] < 0.5 || (K[7] < 1.5 ? ascRun : (ascRun || prevRunL == 0));
        if (K[2] > 0 && arenaGate && L <= MAXL && runLen >= (int)K[0]) {
            int c = runLen & 1;
            Arena& A = ar[L][c];
            if (A.slots == 0 || A.usedSlots >= A.slots) { dugBlocks = 0; blockOpen(L, c); dug += dugBlocks; }
            if (A.slots > 0 && A.usedSlots < A.slots) {
                ll p = A.start + A.used;
                int gl = g.getLen(p);
                if (gl >= L && g.getOwner(p) == A.tag) {
                    A.used += L;
                    A.usedSlots++;
                    g.cut(p, gl, p, L, A.tag);
                    commit(p, L, arenaId(L, c) + 1);
                    N = dug; D = p;
                    return;
                }
                A.slots = 0;                     // страховка: блок недоступен
            }
        }

        ll s = -1;
        int tp = (int)K[3];
        if (tp > 0 && T <= (ll)K[4]) {           // режим малой траншеи: головоломки
            if (tp == 2 || tp == 3) s = g.exactFit(L);
            if (s < 0) s = (tp == 1 || tp == 2) ? g.worstFit(L) : g.bestFit(L);
        } else {
            s = g.bestFit(L);
        }
        ll tailStart = -1;
        auto itEnd = g.startAtEnd.find(T);
        if (itEnd != g.startAtEnd.end()) tailStart = itEnd->second;
        int tailLen = tailStart >= 0 ? g.getLen(tailStart) : -1;
        if (tailLen >= L && s >= 0 && tailLen == g.getLen(s)) s = tailStart;

        if (s >= 0) {
            int gl = g.getLen(s);
            bool toLeft = sideChoice(s, gl, L);
            ll p = toLeft ? s : (s + gl - L);
            g.cut(s, gl, p, L, 0);
            commit(p, L, 0);
            N = dug; D = p;
            return;
        }
        ll need = L - std::max(0, tailLen);
        ll nb = need <= 0 ? 0 : (need + 19) / 20;
        if (nb > 0) {
            if (tailLen > 0) { ll ts = tailStart; int tl = tailLen; g.delGap(ts); g.addGap(ts, tl + nb * 20, 0); }
            else { tailStart = T; g.addGap(tailStart, nb * 20, 0); }
            T += nb * 20;
            dug += nb;
        }
        int gl = g.getLen(tailStart);
        g.cut(tailStart, gl, tailStart, L, 0);
        commit(tailStart, L, 0);
        N = dug; D = tailStart;
    }

    void commit(ll p, int L, int arenaTag) {
        int pid = (int)pos.size();
        pos.push_back(p);
        plen.push_back(L);
        parena.push_back(arenaTag);
        pipeAt[p] = pid;
    }

    bool sideChoice(ll s, int gl, int L) {
        if (s == 0) return true;
        int leftPid = pidAt(s - 1);
        int rightPid = pidAt(s + gl);
        if (leftPid == 0) return false;
        if (rightPid == 0) return true;
        return rightPid > leftPid;
    }

    // ------------------------------------------------------------------ арены
    void closeArena(int L, int c) {
        Arena& A = ar[L][c];
        if (A.slots > 0 && A.usedSlots < A.slots) {
            ll s = A.start + A.used;
            if (g.getLen(s) >= 0 && g.getOwner(s) == A.tag) {
                int len = g.getLen(s);
                g.delGap(s);
                g.addGap(s, len, 0);          // возврат остатка блока в общий фонд
            }
        }
        A = Arena();
    }

    void closeArenasOf(int L) {
        for (int c = 0; c < NCLASS; c++) closeArena(L, c);
        arenaOn[L] = false;
    }

    // Открыть блок арены.
    //   Ёмкость ограничена по КЛЕТКАМ (K[5]): большие блоки простаивают, пока зарезервированы, и заставляют
    //   копать в другом месте (проверено на desc_groups: 16 слотов × 40 м → KPI 0.71 вместо 0.97).
    //   Проверено и отвергнуто: подстройка блока под найденную дыру (slots = gl/L) — тест 5 падает
    //   0.748 → 0.679: классы размазываются по мелким дырам и теряют непрерывность.
    bool blockOpen(int L, int c) {
        Arena& A = ar[L][c];
        int capCells = std::max(20, (int)K[5]);
        int maxSlots = std::max(1, (int)K[6]);
        int slots = std::max(1, std::min(maxSlots, capCells / L));
        int tag = nextTag++;
        ll s = g.bestFit(L * slots);
        if (s < 0 || g.getLen(s) < L * slots) {          // нет готового промежутка — копаем в хвосте
            int fit = std::max(1, (int)K[6]);
            ll want = (ll)slots * L;
            ll tailStart = -1;
            auto itEnd = g.startAtEnd.find(T);
            if (itEnd != g.startAtEnd.end()) tailStart = itEnd->second;
            int tailLen = tailStart >= 0 ? g.getLen(tailStart) : -1;
            ll need = want - std::max(0, tailLen);
            if (need > 0) {
                ll nb = (need + 19) / 20;
                if (tailLen > 0) { ll ts = tailStart; int tl = tailLen; g.delGap(ts); g.addGap(ts, tl + nb * 20, 0); }
                else { tailStart = T; g.addGap(tailStart, nb * 20, 0); }
                T += nb * 20;
                dugBlocks += nb;
            }
            s = tailStart;
            (void)fit;
        }
        ll want = (ll)slots * L;
        int gl = g.getLen(s);
        if (gl < want) { A.slots = 0; return false; }
        g.delGap(s);
        g.addGap(s, want, tag);
        if (gl > want) g.addGap(s + want, gl - want, 0);
        A.start = s;
        A.len = want;
        A.used = 0;
        A.slots = slots;
        A.usedSlots = 0;
        A.blocks++;
        A.tag = tag;
        arenaOn[L] = true;
        return true;
    }

    // ------------------------------------------------------------------ удаление
    ll remove(int pid) {
        ll p = pos[pid];
        ll len = plen[pid];
        int aid = parena[pid];
        if (aid) {                              // труба из арены: освобождаем остаток блока
            int id = aid - 1;
            int L = id / NCLASS, c = id % NCLASS;
            if (L >= 1 && L <= MAXL) closeArenasOf(L);
        }
        ll s = p, e = p + len;
        auto itEnd = g.startAtEnd.find(s);
        if (itEnd != g.startAtEnd.end()) { s = itEnd->second; g.delGap(s); }
        int rl = g.getLen(e);
        if (rl >= 0) { g.delGap(e); e += rl; }
        g.addGap(s, e - s, 0);
        pipeAt.erase(p);
        runLen = 0;
        return p;
    }
};

// ------------------------------------------------------------------ ввод-вывод
// ВАЖНО: ввод — через read(), а не fread(stdin): fread с большим буфером НА ПАЙПЕ блокируется до
// заполнения буфера целиком (интерактор ждёт ответа, решение ждёт 64 КБ → TL). fread здесь допустим
// только как «читаю столько, сколько дали», поэтому используем системный вызов.
#ifdef _WIN32
#include <io.h>
static inline int fdRead(int fd, void* buf, unsigned n) { return _read(fd, buf, n); }
static inline int fdWrite(int fd, const void* buf, unsigned n) { return _write(fd, buf, n); }
#else
#include <unistd.h>
static inline int fdRead(int fd, void* buf, unsigned n) { return (int)read(fd, buf, n); }
static inline int fdWrite(int fd, const void* buf, unsigned n) { return (int)write(fd, buf, n); }
#endif

static char ibuf[1 << 16];
static int ipos = 0, ilen = 0;
static inline int gc() {
    if (ipos == ilen) {
        ilen = fdRead(0, ibuf, sizeof(ibuf));
        ipos = 0;
        if (ilen <= 0) return -1;
    }
    return ibuf[ipos++];
}
static bool readInt(int& x) {
    int c = gc();
    while (c == ' ' || c == '\n' || c == '\r' || c == '\t') c = gc();
    if (c < 0) return false;
    bool neg = false;
    if (c == '-') { neg = true; c = gc(); }
    ll v = 0;
    while (c >= '0' && c <= '9') { v = v * 10 + (c - '0'); c = gc(); }
    x = (int)(neg ? -v : v);
    return true;
}
static inline int fmtNum(ll v, char* b) {
    char t[24];
    int m = 0;
    if (v == 0) t[m++] = '0';
    while (v > 0) { t[m++] = (char)('0' + v % 10); v /= 10; }
    int n = 0;
    while (m) b[n++] = t[--m];
    return n;
}
static void outOne(ll a) {
    char b[32];
    int n = fmtNum(a, b);
    b[n++] = '\n';
    fdWrite(1, b, n);            // один системный вызов на ответ = и есть требуемый flush
}
static void outTwo(ll a, ll b) {
    char buf[48];
    int n = fmtNum(a, buf);
    buf[n++] = ' ';
    n += fmtNum(b, buf + n);
    buf[n++] = '\n';
    fdWrite(1, buf, n);
}

#ifndef SIM
int main() {
    for (int i = 0; i < 64; i++) K[i] = KD[i];
    Engine E;
    int c;
    while (readInt(c)) {
        if (c == 0) return 0;
        if (c > 0) {
            ll N, D;
            E.place(c, N, D);
            outTwo(N, D);
        } else {
            outOne(E.remove(-c));
        }
    }
    return 0;
}
#endif
