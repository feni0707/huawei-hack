// Вариант A: «чистые» аллокаторные политики + мягкая сегрегация по размерам.
//
// Состояние: множество свободных промежутков [s, s+len) внутри траншеи [0, T). Соседние промежутки
// всегда склеены. На укладку трубы длины L выбирается промежуток по политике (зависит от того,
// «короткая» это труба (L <= thr) или «длинная»); если подходящего нет — траншея докапывается
// минимально (с использованием свободного хвоста) и труба кладётся в начало хвоста.
//
// Политики выбора промежутка (pol):
//   0 best   — наименьший промежуток длины >= L (среди равных — самый левый / правый, tie)
//   1 first  — самый левый подходящий
//   2 last   — самый правый подходящий
//   3 worst  — самый большой
//   4 bwl    — самый левый среди промежутков с запасом (len - L) <= w; если таких нет — best
//   5 bwr    — самый правый среди промежутков с запасом <= w; если таких нет — best
//   6 sliver — best fit со штрафом: стоимость = (len - L) + p, если остаток 1..m-1 («щепка»)
// Сторона в промежутке (side): 0 — к левому краю, 1 — к правому.
//
// Сборка решения:  c++ -std=c++17 -O2 -o build/variant_a solutions/variant_a.cpp
// Сборка лаборатории (оффлайн-симулятор, без интерактора):
//   c++ -std=c++17 -O2 -DLAB -o build/lab_a solutions/variant_a.cpp
//   build/lab_a [key=val ...] test1.in test2.in ...    → строка на тест: имя kpi T occ peak ub
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <vector>
using namespace std;
typedef long long ll;

struct Params {
    int polS = 0, polL = 0;      // политика для коротких (L <= thr) и длинных труб
    int thr = 40;                // порог: L <= thr — «короткая»
    int tieS = 0, tieL = 0;      // best fit: среди равных промежутков 0 — самый левый, 1 — самый правый
    int sideS = 0, sideL = 0;    // сторона внутри найденного промежутка
    int wS = 0, wL = 0;          // допустимый запас для политик 4/5
    int ext = 0;                 // докопка: 0 — труба в начало хвоста, 1 — труба прижата к новому концу
    int lanes = 1;               // A+: 1 — две полосы по чётности индекса в серии равных длин
    int laneCells = 320;         // A+: размер резервируемого участка полосы в клетках (постоянный)
    int laneR0 = 2;              // A+: полосы включаются, когда серия равных длин подряд длиннее R0
    int learn = 1;               // A+: 1 — полосы только после подтверждённой закономерности удалений в прошлой серии
    int laneFrac = 25;            // A+: если нет дыры >= R, взять самую большую, если она >= laneFrac % от R (0 — не брать)
    int laneMinRun = 8;          // A+: минимальная длина серии, по которой изучаем закономерность
    int mS = 0, mL = 0;          // политика 6: остаток 1..m-1 считается «щепкой» и штрафуется
    int pS = 0, pL = 0;          // политика 6: штраф за щепку (в клетках)
};

class Allocator {
    static const int BIG = 40;   // промежутки длины >= BIG вмещают любую трубу (L <= 40)
    Params P;
    ll T = 0;                    // длина траншеи
    map<ll, ll> gaps;            // start -> len, все свободные промежутки (максимальные)
    set<ll> cls[BIG];            // cls[c] — начала промежутков длины ровно c (1..39)
    set<ll> bigStart;            // начала промежутков длины >= 40
    set<pair<ll, ll>> bigLen;    // (len, start) для промежутков длины >= 40
    // A+ (полосы): состояние серии равных длин и зарезервированные участки двух полос
    bool lastPlace = false;
    ll lastLen = 0;
    int runLen = 0;
    ll laneCur[2] = {0, 0}, laneEnd[2] = {0, 0};
    // закономерность удалений внутри серии: индекс i (с нуля) в серии умирает, если i % lm == lr
    bool lvalid = false;
    int lm = 2, lr = 0;
    ll runStart = 0;              // PID первой трубы текущей серии
    ll prevStart = 0;             // PID первой трубы последней завершённой серии
    int prevN = 0;                // её длина (0 — нет серии, ожидающей оценки)
    vector<char> died;            // died[i] — i-я труба последней серии удалена
    int diedCnt = 0;
    bool prevEndedByRemoval = false;
public:
    ll live = 0, peak = 0;
    vector<ll> pos, len;
    vector<char> alive;

    explicit Allocator(const Params& p) : P(p), pos(1, -1), len(1, 0), alive(1, 0) {}
    ll trench() const { return T; }

private:
    void insertGap(ll s, ll l) {
        gaps[s] = l;
        if (l < BIG) cls[l].insert(s);
        else { bigStart.insert(s); bigLen.insert({l, s}); }
    }
    void eraseGap(ll s, ll l) {
        gaps.erase(s);
        if (l < BIG) cls[l].erase(s);
        else { bigStart.erase(s); bigLen.erase({l, s}); }
    }
    // Добавить свободный отрезок с немедленной склейкой соседей.
    void addFree(ll s, ll l) {
        if (l <= 0) return;
        auto it = gaps.lower_bound(s);
        if (it != gaps.end() && it->first == s + l) {
            ll ns = it->first, nl = it->second;
            eraseGap(ns, nl);
            l += nl;
        }
        it = gaps.lower_bound(s);
        if (it != gaps.begin()) {
            auto pv = prev(it);
            if (pv->first + pv->second == s) {
                ll ps = pv->first, pl = pv->second;
                eraseGap(ps, pl);
                s = ps;
                l += pl;
            }
        }
        insertGap(s, l);
    }

    static ll better(bool right, ll a, ll b) {   // a — текущий (может быть -1), b — кандидат
        if (a < 0) return b;
        return right ? max(a, b) : min(a, b);
    }

    bool selectBest(ll L, int tie, ll& gs, ll& gl) {
        for (int c = (int)L; c < BIG; c++)
            if (!cls[c].empty()) {
                gs = tie ? *cls[c].rbegin() : *cls[c].begin();
                gl = c;
                return true;
            }
        if (!bigLen.empty()) {
            auto it = bigLen.begin();
            ll bl = it->first;
            if (!tie) gs = it->second;
            else gs = prev(bigLen.lower_bound({bl + 1, -1}))->second;
            gl = bl;
            return true;
        }
        return false;
    }
    // first (right=false) / last (right=true): крайний по адресу из всех подходящих
    bool selectByAddress(ll L, bool right, ll& gs, ll& gl) {
        ll best = -1, bl = 0;
        for (int c = (int)L; c < BIG; c++)
            if (!cls[c].empty()) {
                ll s = right ? *cls[c].rbegin() : *cls[c].begin();
                ll nb = better(right, best, s);
                if (nb != best) { best = nb; bl = c; }
            }
        if (!bigStart.empty()) {
            ll s = right ? *bigStart.rbegin() : *bigStart.begin();
            ll nb = better(right, best, s);
            if (nb != best) { best = nb; bl = gaps[s]; }
        }
        if (best < 0) return false;
        gs = best; gl = bl;
        return true;
    }
    bool selectWorst(ll L, ll& gs, ll& gl) {
        if (!bigLen.empty()) {
            auto it = prev(bigLen.end());
            gl = it->first; gs = it->second;
            return true;
        }
        for (int c = BIG - 1; c >= (int)L; c--)
            if (!cls[c].empty()) { gl = c; gs = *cls[c].begin(); return true; }
        return false;
    }
    // крайний по адресу среди промежутков с запасом len - L <= w; иначе best fit
    bool selectBounded(ll L, int w, bool right, int tie, ll& gs, ll& gl) {
        ll best = -1, bl = 0;
        ll hi = min<ll>(BIG - 1, L + w);
        for (ll c = L; c <= hi; c++)
            if (!cls[c].empty()) {
                ll s = right ? *cls[c].rbegin() : *cls[c].begin();
                ll nb = better(right, best, s);
                if (nb != best) { best = nb; bl = c; }
            }
        if (L + w >= BIG && !bigLen.empty()) {   // большие промежутки с запасом <= w (ограниченный просмотр)
            auto it = bigLen.lower_bound({max<ll>(BIG, L), -1});
            for (int k = 0; k < 48 && it != bigLen.end() && it->first <= L + w; ++k, ++it) {
                ll nb = better(right, best, it->second);
                if (nb != best) { best = nb; bl = it->first; }
            }
        }
        if (best >= 0) { gs = best; gl = bl; return true; }
        return selectBest(L, tie, gs, gl);
    }
    // best fit со штрафом за «щепки»: стоимость = остаток + p, если 0 < остаток < m
    bool selectSliver(ll L, int m, int pen, int tie, ll& gs, ll& gl) {
        ll bestCost = (ll)4e18, bs = -1, bl = 0;
        for (int c = (int)L; c < BIG; c++) {
            if (cls[c].empty()) continue;
            ll r = c - L;
            ll cost = r + ((r > 0 && r < m) ? pen : 0);
            if (cost < bestCost) {
                bestCost = cost;
                bs = tie ? *cls[c].rbegin() : *cls[c].begin();
                bl = c;
            }
        }
        if (!bigLen.empty()) {
            auto it = bigLen.begin();
            ll r = it->first - L;
            ll cost = r + ((r > 0 && r < m) ? pen : 0);
            if (cost < bestCost) {
                bestCost = cost;
                bl = it->first;
                bs = tie ? prev(bigLen.lower_bound({bl + 1, -1}))->second : it->second;
            }
        }
        if (bs < 0) return false;
        gs = bs; gl = bl;
        return true;
    }
    bool select(ll L, ll& gs, ll& gl) {
        bool small = L <= P.thr;
        int pol = small ? P.polS : P.polL;
        int tie = small ? P.tieS : P.tieL;
        int w = small ? P.wS : P.wL;
        switch (pol) {
            case 6: return selectSliver(L, small ? P.mS : P.mL, small ? P.pS : P.pL, tie, gs, gl);
            case 0: return selectBest(L, tie, gs, gl);
            case 1: return selectByAddress(L, false, gs, gl);
            case 2: return selectByAddress(L, true, gs, gl);
            case 3: return selectWorst(L, gs, gl);
            case 4: return selectBounded(L, w, false, tie, gs, gl);
            default: return selectBounded(L, w, true, tie, gs, gl);
        }
    }

    // Серия завершилась (пришла команда не-«та же длина»): запомнить её для оценки закономерности удалений.
    void endRun(bool byRemoval) {
        releaseLanes();
        if (runLen >= P.laneMinRun) {
            prevStart = runStart;
            prevN = runLen;
            died.assign(runLen, 0);
            diedCnt = 0;
            prevEndedByRemoval = byRemoval;
        } else if (!byRemoval) {
            // короткая серия: прошлую оценку не трогаем
        }
        runLen = 0;
    }
    // Вызывается при первой укладке после пачки удалений: оценить, как умирали трубы последней серии.
    void evaluatePattern() {
        if (prevN == 0 || !prevEndedByRemoval) { prevN = 0; return; }
        int n = prevN, c = diedCnt;
        prevN = 0;
        lvalid = false;
        if (c < max(4, n / 4)) return;
        for (int m = 2; m <= 4; m++)
            for (int r = 0; r < m; r++) {
                int in = 0, tot = 0;
                for (int i = 0; i < n; i++)
                    if (i % m == r) { tot++; in += died[i]; }
                if (tot >= 4 && in * 10 >= tot * 9 && (c - in) * 10 <= c) {
                    lvalid = true; lm = m; lr = r;
                    return;
                }
            }
    }
    void releaseLanes() {
        for (int c = 0; c < 2; c++) {
            if (laneEnd[c] > laneCur[c]) addFree(laneCur[c], laneEnd[c] - laneCur[c]);
            laneCur[c] = laneEnd[c] = 0;
        }
    }
    // Зарезервировать непрерывный участок размера R: из наименьшего достаточного большого промежутка либо докопкой.
    // Возвращает начало участка; N — сколько 20-метровых блоков докопано.
    ll reserveRegion(ll R, ll& N, ll& Reff) {
        N = 0;
        Reff = R;
        auto it = bigLen.lower_bound({R, -1});
        if (it != bigLen.end()) {
            ll gl = it->first, gs = it->second;
            eraseGap(gs, gl);
            addFree(gs + R, gl - R);
            return gs;
        }
        if (P.laneFrac > 0 && !bigLen.empty()) {   // нет достаточной дыры: взять самую большую, если она не слишком мала
            auto jt = prev(bigLen.end());
            if (jt->first * 100 >= R * P.laneFrac) {
                ll gl = jt->first, gs = jt->second;
                eraseGap(gs, gl);
                Reff = gl;
                return gs;
            }
        }
        ll ts = T, tl = 0;
        if (!gaps.empty()) {
            auto last = prev(gaps.end());
            if (last->first + last->second == T) { ts = last->first; tl = last->second; }
        }
        N = (R - tl + 19) / 20;
        if (tl > 0) eraseGap(ts, tl);
        T += N * 20;
        addFree(ts + R, T - (ts + R));
        return ts;
    }

public:
    // Укладка трубы длины L: возвращает (N, D).
    pair<ll, ll> place(ll L) {
        if (P.lanes) {
            if (lastPlace && L == lastLen) {
                runLen++;
            } else {
                if (!lastPlace) evaluatePattern();   // пачка удалений закончилась
                if (lastPlace) endRun(false);        // серия оборвалась сменой длины
                runLen = 1;
                runStart = (ll)pos.size();            // PID этой трубы
            }
            lastLen = L;
            lastPlace = true;
            int idx = runLen - 1;
            int c = -1;
            if (P.learn) {
                if (lvalid && runLen > P.laneR0) c = (idx % lm == lr) ? 0 : 1;
            } else if (runLen > P.laneR0) {
                c = idx & 1;
            }
            if (c >= 0) {
                ll N = 0;
                if (laneCur[c] + L > laneEnd[c]) {
                    if (laneEnd[c] > laneCur[c]) addFree(laneCur[c], laneEnd[c] - laneCur[c]);
                    ll R = P.laneCells, Reff = R;
                    laneCur[c] = reserveRegion(R, N, Reff);
                    laneEnd[c] = laneCur[c] + Reff;
                }
                ll D = laneCur[c];
                laneCur[c] += L;
                pos.push_back(D);
                len.push_back(L);
                alive.push_back(1);
                live += L;
                if (live > peak) peak = live;
                return {N, D};
            }
        }
        bool small = L <= P.thr;
        int side = small ? P.sideS : P.sideL;
        ll gs, gl, N = 0, D;
        if (select(L, gs, gl)) {
            eraseGap(gs, gl);
            if (side == 0) { D = gs; addFree(gs + L, gl - L); }
            else { D = gs + gl - L; addFree(gs, gl - L); }
        } else {
            ll ts = T, tl = 0;
            if (!gaps.empty()) {
                auto last = prev(gaps.end());
                if (last->first + last->second == T) { ts = last->first; tl = last->second; }
            }
            N = (L - tl + 19) / 20;
            if (tl > 0) eraseGap(ts, tl);
            T += N * 20;
            if (P.ext == 0) { D = ts; addFree(D + L, T - (D + L)); }
            else { D = T - L; addFree(ts, D - ts); }
        }
        pos.push_back(D);
        len.push_back(L);
        alive.push_back(1);
        live += L;
        if (live > peak) peak = live;
        return {N, D};
    }
    // Извлечение трубы PID: возвращает её позицию.
    ll remove(ll pid) {
        if (P.lanes) {
            if (lastPlace) { endRun(true); lastPlace = false; }
            if (prevN > 0 && pid >= prevStart && pid < prevStart + prevN) {
                died[pid - prevStart] = 1;
                diedCnt++;
            }
        }
        ll s = pos[pid], l = len[pid];
        alive[pid] = 0;
        live -= l;
        addFree(s, l);
        return s;
    }
};

#ifdef LAB
// ---------------------------------------------------------------- лаборатория
// Картинка раскладки для маленьких тестов: '.' — свободно, буква — труба (a + PID mod 26).
static string renderLayout(const Allocator& A) {
    string r(A.trench(), '.');
    for (size_t pid = 1; pid < A.alive.size(); pid++)
        if (A.alive[pid])
            for (ll i = 0; i < A.len[pid]; i++) r[A.pos[pid] + i] = 'a' + (pid % 26);
    return r;
}
static int gVis = 0;
static bool parseParam(Params& P, const char* arg) {
    const char* eq = strchr(arg, '=');
    if (!eq) return false;
    string k(arg, eq - arg);
    int v = atoi(eq + 1);
    if (k == "polS") P.polS = v;
    else if (k == "polL") P.polL = v;
    else if (k == "pol") P.polS = P.polL = v;
    else if (k == "thr") P.thr = v;
    else if (k == "tieS") P.tieS = v;
    else if (k == "tieL") P.tieL = v;
    else if (k == "tie") P.tieS = P.tieL = v;
    else if (k == "sideS") P.sideS = v;
    else if (k == "sideL") P.sideL = v;
    else if (k == "side") P.sideS = P.sideL = v;
    else if (k == "wS") P.wS = v;
    else if (k == "wL") P.wL = v;
    else if (k == "w") P.wS = P.wL = v;
    else if (k == "ext") P.ext = v;
    else if (k == "vis") gVis = v;
    else if (k == "lanes") P.lanes = v;
    else if (k == "laneCells") P.laneCells = v;
    else if (k == "laneR0") P.laneR0 = v;
    else if (k == "learn") P.learn = v;
    else if (k == "laneFrac") P.laneFrac = v;
    else if (k == "laneMinRun") P.laneMinRun = v;
    else if (k == "mS") P.mS = v;
    else if (k == "mL") P.mL = v;
    else if (k == "m") P.mS = P.mL = v;
    else if (k == "pS") P.pS = v;
    else if (k == "pL") P.pL = v;
    else if (k == "p") P.pS = P.pL = v;
    else { fprintf(stderr, "unknown param %s\n", k.c_str()); exit(2); }
    return true;
}

int main(int argc, char** argv) {
    Params P;
    vector<const char*> files;
    for (int i = 1; i < argc; i++)
        if (!parseParam(P, argv[i])) files.push_back(argv[i]);
    for (const char* fn : files) {
        FILE* f = fopen(fn, "r");
        if (!f) { fprintf(stderr, "cannot open %s\n", fn); return 2; }
        double w; ll n;
        if (fscanf(f, "%lf %lld", &w, &n) != 2) return 2;
        Allocator A(P);
        ll c;
        for (ll i = 0; i < n && fscanf(f, "%lld", &c) == 1; i++) {
            if (c > 0) A.place(c);
            else A.remove(-c);
        }
        fclose(f);
        if (gVis) printf("%s\n", renderLayout(A).c_str());
        ll T = A.trench();
        ll ubT = 20 * ((A.peak + 19) / 20);
        printf("%s %.6f %lld %lld %lld %.6f\n", fn, T ? (double)A.live / T : 1.0, T, A.live, A.peak,
               ubT ? (double)A.live / ubT : 1.0);
    }
    return 0;
}
#else
// ---------------------------------------------------------------- решение (интерактивный протокол)
int main() {
    Params P;   // значения по умолчанию = выбранная конфигурация (см. EXPERIMENTS.md)
    Allocator A(P);
    ll x;
    while (scanf("%lld", &x) == 1) {
        if (x == 0) break;
        if (x > 0) {
            pair<ll, ll> r = A.place(x);
            printf("%lld %lld\n", r.first, r.second);
        } else {
            printf("%lld\n", A.remove(-x));
        }
        fflush(stdout);
    }
    return 0;
}
#endif
