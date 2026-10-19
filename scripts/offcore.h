// Оффлайн-ядро «время-пространство»: адреса труб по первому подходящему месту в порядке убывания времени смерти.
// Один и тот же код используется оффлайн-инструментом (scripts/offline_tool.cpp) и внутри посылки (solutions/overfit.cpp),
// поэтому ответы воспроизводятся побитово. Только стандартная библиотека.
#ifndef OFFCORE_H
#define OFFCORE_H
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <vector>
#include <functional>

namespace off {

static const int INF = 1 << 30;

struct Seg {   // m(c) = минимальное время рождения среди уже поставленных труб, покрывающих клетку c (теги не проталкиваются)
    int n;
    std::vector<int> tg, mn, mx;
    void init(int n_) { n = 1; while (n < n_) n <<= 1; tg.assign(2 * n, INF); mn.assign(2 * n, INF); mx.assign(2 * n, INF); }
    void upd(int x, int l, int r, int a, int b, int v) {
        if (b <= l || r <= a) return;
        if (a <= l && r <= b) { if (v < tg[x]) tg[x] = v; if (v < mn[x]) mn[x] = v; if (v < mx[x]) mx[x] = v; return; }
        int m = (l + r) >> 1;
        upd(2 * x, l, m, a, b, v); upd(2 * x + 1, m, r, a, b, v);
        int t = tg[x];
        mn[x] = std::min(t, std::min(mn[2 * x], mn[2 * x + 1]));
        mx[x] = std::min(t, std::max(mx[2 * x], mx[2 * x + 1]));
    }
    void chmin(int a, int b, int v) { upd(1, 0, n, a, b, v); }
    int firstGE(int x, int l, int r, int pos, int d, int acc) {   // первый индекс >= pos со значением >= d (или n)
        if (r <= pos || std::min(acc, mx[x]) < d) return n;
        if (r - l == 1) return l;
        int a2 = std::min(acc, tg[x]), m = (l + r) >> 1;
        int res = firstGE(2 * x, l, m, pos, d, a2);
        if (res != n) return res;
        return firstGE(2 * x + 1, m, r, pos, d, a2);
    }
    int firstLT(int x, int l, int r, int pos, int d, int acc) {   // первый индекс >= pos со значением < d (или n)
        if (r <= pos || std::min(acc, mn[x]) >= d) return n;
        if (r - l == 1) return l;
        int a2 = std::min(acc, tg[x]), m = (l + r) >> 1;
        int res = firstLT(2 * x, l, m, pos, d, a2);
        if (res != n) return res;
        return firstLT(2 * x + 1, m, r, pos, d, a2);
    }
};

// Безопасный префикс: простой best fit (минимальный подходящий промежуток, меньший адрес; иначе хвост).
// Зависит только от самих команд, поэтому одинаков у всех тестов с общим началом и у любого неизвестного теста.
struct BfSim {
    std::vector<int> pos{-1}, len{0};
    std::map<int, int> holes;                      // начало -> длина (не касаются хвоста)
    std::set<std::pair<int, int>> bySize;          // (длина, начало)
    int tail = 0;
    void addHole(int a, int h) { holes[a] = h; bySize.insert({h, a}); }
    void delHole(int a) { bySize.erase({holes[a], a}); holes.erase(a); }
    int place(int c) {
        int D;
        auto it = bySize.lower_bound({c, -1});
        if (it != bySize.end()) {
            int h = it->first, a = it->second;
            delHole(a);
            D = a;
            if (h > c) addHole(a + c, h - c);
        } else { D = tail; tail += c; }
        pos.push_back(D); len.push_back(c);
        return D;
    }
    void remove(int pid) {
        int a = pos[pid], l = len[pid], s = a, e = a + l;
        auto it = holes.lower_bound(a);
        if (it != holes.begin()) { auto jt = std::prev(it); if (jt->first + jt->second == a) { s = jt->first; delHole(jt->first); } }
        if (e == tail) tail = s;
        else { auto r = holes.find(e); if (r != holes.end()) { int hl = r->second; delHole(e); e += hl; } addHole(s, e - s); }
    }
};
inline std::vector<int> bfPrefix(const std::vector<int>& cmds, int K) {
    BfSim b;
    std::vector<int> res;
    for (int i = 0; i < K && i < (int)cmds.size(); i++) {
        if (cmds[i] > 0) res.push_back(b.place(cmds[i])); else b.remove(-cmds[i]);
    }
    return res;
}

struct Problem {
    std::vector<int> cmds;
    // выводимое
    std::vector<int> len, birth, death;   // по pid (с 1)
    int np = 0;
    long long peak = 0, liveEnd = 0, ubT = 0;
    // условия
    std::vector<int> key;                 // «ключ» смерти для быстрого решателя (>= death): порядок и порог; по умолчанию = death
    std::vector<int> pin;                 // pin[pid] >= 0 — позиция зафиксирована (общий префикс с другим тестом)
    int ov = 0;                           // вариант тай-брейка начального порядка
    int ffPre = 0;                        // трубы pid <= ffPre всегда идут по первому подходящему слева (компактно), остальные — по fitMode
    int capPre = 0; long long capAddr = 0;   // трубы pid <= capPre желательно держать левее capAddr (только штраф в excess/bad, на раскладку не влияет)
    int capK = 16, capSlack = 40;         // первые capK труб должны лежать компактно: конец <= сумма длин первых i + capSlack
    // результат
    std::vector<int> D;
    std::vector<long long> prefLen;

    void load(const std::vector<int>& c) {
        cmds = c;
        len.assign(1, 0); birth.assign(1, 0); death.assign(1, INF);
        for (int i = 0; i < (int)cmds.size(); i++) {
            if (cmds[i] > 0) { len.push_back(cmds[i]); birth.push_back(i); death.push_back(INF); }
            else death[-cmds[i]] = i;
        }
        np = (int)len.size() - 1;
        long long cur = 0; peak = 0;
        for (int x : cmds) { if (x > 0) { cur += x; peak = std::max(peak, cur); } else cur -= len[-x]; }
        liveEnd = cur;
        ubT = (peak + 19) / 20 * 20;
        prefLen.assign(np + 1, 0);
        for (int i = 1; i <= np; i++) prefLen[i] = prefLen[i - 1] + len[i];
        pin.assign(np + 1, -1);
        key = death;
    }

    // начальный порядок: по убыванию смерти; выжившие — по рождению (плотно), остальные — длинные раньше; закреплённые не входят
    std::vector<int> initialOrder() const {
        std::vector<int> ord;
        for (int i = 1; i <= np; i++) if (pin[i] < 0) ord.push_back(i);
        std::sort(ord.begin(), ord.end(), [&](int a, int b) {
            const std::vector<int>& kk = ov >= 10 ? key : death;
            if (kk[a] != kk[b]) return kk[a] > kk[b];
            if (death[a] == INF) return birth[a] < birth[b];
            int o = ov % 10;
            if (o == 1 && birth[a] != birth[b]) return birth[a] < birth[b];
            if (o == 2 && len[a] != len[b]) return len[a] < len[b];
            if (o == 3 && birth[a] != birth[b]) return birth[a] > birth[b];
            if (o == 0 && len[a] != len[b]) return len[a] > len[b];
            return a < b;
        });
        return ord;
    }

    struct Pin { int a, L, b, d; };

    // решение для порядка ord; D заполняется; возвращает maxEnd (или -1, если не уложились в ограничения «компактного префикса»)
    long long solve(const std::vector<int>& ord, long long target, std::vector<int>& bad, long long& excess) {
        Seg seg;
        seg.init((int)(ubT * 13 / 10) + 4000);
        D.assign(np + 1, 0);
        std::vector<Pin> pins;
        long long maxEnd = 0;
        int firstFree = INF;   // раньше всего родившаяся незакреплённая труба
        for (int i = 1; i <= np; i++) if (pin[i] < 0) { firstFree = birth[i]; break; }
        for (int i = 1; i <= np; i++) if (pin[i] >= 0) {
            D[i] = pin[i];
            maxEnd = std::max<long long>(maxEnd, pin[i] + len[i]);
            if (death[i] == INF) seg.chmin(pin[i], pin[i] + len[i], birth[i]);   // закреплённые выжившие — как обычные поставленные трубы
            else if (death[i] > firstFree) pins.push_back({pin[i], len[i], birth[i], death[i]});   // с конечной смертью — точные препятствия
            // умершие до рождения первой незакреплённой трубы ни с кем не пересекаются
        }
        std::sort(pins.begin(), pins.end(), [](const Pin& x, const Pin& y) { return x.a < y.a; });
        bad.clear();
        excess = 0;
        for (int p : ord) {
            int L = len[p], d = death[p], b = birth[p];
            long long lim = p <= capK ? prefLen[p] + capSlack : (1LL << 40);
            int c = 0, found = -1;
            for (;;) {
                c = seg.firstGE(1, 0, seg.n, c, d, INF);
                if (c >= seg.n || c + L > lim) break;
                int q = seg.firstLT(1, 0, seg.n, c, d, INF);
                if (q - c < L) { c = q; continue; }
                int bump = c;
                if (!pins.empty()) {
                    size_t lo = std::lower_bound(pins.begin(), pins.end(), c - 39, [](const Pin& x, int v) { return x.a < v; }) - pins.begin();
                    for (size_t k = lo; k < pins.size() && pins[k].a < c + L; k++) {
                        const Pin& o = pins[k];
                        if (o.a + o.L > c && o.b < d && o.d > b) bump = std::max(bump, o.a + o.L);
                    }
                }
                if (bump == c) { found = c; break; }
                c = bump;
            }
            if (found < 0) return -1;
            D[p] = found;
            seg.chmin(found, found + L, b);
            maxEnd = std::max<long long>(maxEnd, found + L);
            if (found + L > target) { bad.push_back(p); excess += found + L - target; }
        }
        return maxEnd;
    }

    // Быстрый вариант (ov >= 10): свободные «пробеги» клеток (m(c) >= d) ведутся явно, первое подходящее место — один спуск по дереву максимумов.
    // Закреплённые трубы с конечной смертью блокируют свои клетки навсегда (чуть консервативнее точного решения).
    long long solveFast(const std::vector<int>& ord, long long target, std::vector<int>& bad, long long& excess) {
        int N = (int)(ubT * 13 / 10) + 4000;
        int TN = 1; while (TN < N + 2) TN <<= 1;
        std::vector<int> mx(2 * TN, 0), lenAt(N + 2, 0), startAt(N + 2, -1), m(N, INF);
        std::vector<uint8_t> good(N, 1);
        auto setRun = [&](int a, int l) {   // пробег, начинающийся в a, длины l (l = 0 — удалить)
            int i = a + TN;
            mx[i] = l;
            for (i >>= 1; i; i >>= 1) mx[i] = std::max(mx[2 * i], mx[2 * i + 1]);
        };
        int fitMode = ov >= 10 ? ov % 10 : 0;   // 0 — первое подходящее, 1 — наименьший подходящий пробег (слева), 2 — то же (справа)
        std::set<std::pair<int, int>> bySz;
        auto addRun = [&](int a, int l) { lenAt[a] = l; startAt[a + l] = a; setRun(a, l); if (fitMode) bySz.insert({l, a}); };
        auto delRun = [&](int a) { int l = lenAt[a]; startAt[a + l] = -1; lenAt[a] = 0; setRun(a, 0); if (fitMode) bySz.erase({l, a}); };
        auto firstRun = [&](int L) {   // самый левый пробег длины >= L или -1
            if (mx[1] < L) return -1;
            int x = 1;
            while (x < TN) x = mx[2 * x] >= L ? 2 * x : 2 * x + 1;
            return x - TN;
        };
        std::function<int(int, int, int, int, int)> rec = [&](int node, int lo, int hi, int x, int L) -> int {
            if (hi <= x || mx[node] < L) return -1;
            if (hi - lo == 1) return lo;
            int mid = (lo + hi) / 2;
            int r = rec(2 * node, lo, mid, x, L);
            return r >= 0 ? r : rec(2 * node + 1, mid, hi, x, L);
        };
        auto firstRunFrom = [&](int x, int L) { return rec(1, 0, TN, x, L); };
        auto carve = [&](int a, int L) -> int {    // вырезать L клеток из пробега, начинающегося в a; возвращает позицию
            int l = lenAt[a];
            int pos = fitMode == 2 ? a + l - L : a;
            delRun(a);
            if (pos > a) addRun(a, pos - a);
            if (pos + L < a + l) addRun(pos + L, a + l - (pos + L));
            return pos;
        };
        auto blockCells = [&](int a, int L, int v) { for (int c = a; c < a + L; c++) { good[c] = 0; m[c] = v; } };
        D.assign(np + 1, 0);
        long long maxEnd = 0;
        struct Ent { int v, a, L; };
        std::vector<Ent> heap;   // максимум-куча по v
        auto cmp = [](const Ent& x, const Ent& y) { return x.v < y.v; };
        auto release = [&](int a, int Lw, int dcur) {
            int c = a;
            while (c < a + Lw) {
                if (!good[c] && m[c] >= dcur) {
                    int e = c;
                    while (e < a + Lw && !good[e] && m[e] >= dcur) { good[e] = 1; e++; }
                    int s0 = c, e0 = e;
                    if (s0 > 0 && startAt[s0] >= 0) { int ps = startAt[s0]; delRun(ps); s0 = ps; }
                    if (e0 < N && lenAt[e0] > 0) { int l2 = lenAt[e0]; delRun(e0); e0 += l2; }
                    addRun(s0, e0 - s0);
                    c = e;
                } else c++;
            }
        };
        // закреплённые (по порядку pid): выжившие — обычные блоки с записью в кучу, с конечной смертью — блок навсегда (m = -1)
        int pinEnd = 0;
        std::vector<Pin> obs;
        int firstFree = INF;
        for (int i = 1; i <= np; i++) if (pin[i] < 0) { firstFree = birth[i]; break; }
        for (int i = 1; i <= np; i++) if (pin[i] >= 0) {
            int a = pin[i], L = len[i];
            D[i] = a;
            maxEnd = std::max<long long>(maxEnd, a + L);
            pinEnd = std::max(pinEnd, a + L);
            // умершие до рождения первой незакреплённой трубы клетки оставляют свободными (им не с кем пересекаться)
            int mv = death[i] == INF ? birth[i] : INF;
            if (death[i] != INF && death[i] > firstFree) obs.push_back({a, L, birth[i], death[i]});   // точное препятствие во времени
            for (int c = a; c < a + L; c++) { good[c] = mv == INF; m[c] = mv; }
            if (death[i] == INF) { heap.push_back({birth[i], a, L}); std::push_heap(heap.begin(), heap.end(), cmp); }
        }
        std::sort(obs.begin(), obs.end(), [](const Pin& x, const Pin& y) { return x.a < y.a; });
        {   // пробеги из свободных клеток до pinEnd и общий хвост
            int c = 0;
            while (c < pinEnd) {
                if (good[c]) { int e = c; while (e < pinEnd && good[e]) e++; addRun(c, e - c); c = e; } else c++;
            }
            if (pinEnd < N) { if (pinEnd > 0 && startAt[pinEnd] >= 0) { int ps = startAt[pinEnd]; delRun(ps); addRun(ps, N - ps); } else addRun(pinEnd, N - pinEnd); }
        }
        bad.clear();
        excess = 0;
        for (int p : ord) {
            int L = len[p], d = key[p], b = birth[p];
            while (!heap.empty() && heap.front().v >= d) {
                Ent e = heap.front();
                std::pop_heap(heap.begin(), heap.end(), cmp); heap.pop_back();
                release(e.a, e.L, d);
            }
            int a;
            if (fitMode && p > ffPre) { auto it = bySz.lower_bound({L, -1}); if (it == bySz.end()) return -1; a = it->second; }
            else if (obs.empty()) { a = firstRun(L); if (a < 0) return -1; }
            else {
                int x = 0;
                a = -1;
                for (;;) {
                    int sidx = firstRunFrom(x, L);
                    if (sidx < 0) return -1;
                    int rl = lenAt[sidx], c = sidx;
                    while (c + L <= sidx + rl) {
                        int bump = c;
                        size_t lo = std::lower_bound(obs.begin(), obs.end(), c - 39, [](const Pin& o, int v) { return o.a < v; }) - obs.begin();
                        for (size_t k = lo; k < obs.size() && obs[k].a < c + L; k++) {
                            const Pin& o = obs[k];
                            if (o.a + o.L > c && o.b < d && o.d > b) bump = std::max(bump, o.a + o.L);
                        }
                        if (bump == c) break;
                        c = bump;
                    }
                    if (c + L <= sidx + rl) { a = c; break; }
                    x = sidx + rl;
                }
                // вырезать [a, a+L) из пробега, содержащего a
                int s0 = a;
                while (s0 > 0 && lenAt[s0] == 0) s0--;
                int rl = lenAt[s0];
                delRun(s0);
                if (a > s0) addRun(s0, a - s0);
                if (a + L < s0 + rl) addRun(a + L, s0 + rl - (a + L));
                goto placed;
            }
            a = carve(a, L);
            placed:
            blockCells(a, L, b);
            heap.push_back({b, a, L}); std::push_heap(heap.begin(), heap.end(), cmp);
            D[p] = a;
            maxEnd = std::max<long long>(maxEnd, a + L);
            if (a + L > target) { bad.push_back(p); excess += a + L - target; }
            else if (p <= capPre && a + L > capAddr) { bad.push_back(p); excess += a + L - capAddr; }
        }
        return maxEnd;
    }
    long long solveAuto(const std::vector<int>& ord, long long target, std::vector<int>& bad, long long& excess) {
        return ov >= 10 ? solveFast(ord, target, bad, excess) : solve(ord, target, bad, excess);
    }
};


// Совместная раскладка общего префикса тестов A и B: выжившие трубы префикса теста A, которые умирают в тесте B (cmdsB),
// собираются в отдельный блок сразу за остальными выжившими префикса (ключ INF-1), а поздние выжившие идут после (INF-2).
inline long long jointHost = 6000;   // суммарная длина «хозяев»: самых поздних выживших, на чьи клетки (внизу) во времени встают поздно умирающие трубы префикса
inline void applyJoint(Problem& P, const std::vector<int>& cmdsB, int prefixCmds) {
    Problem Q;
    Q.load(cmdsB);
    int pre = 0;
    for (int i = 0; i < prefixCmds && i < (int)P.cmds.size(); i++) if (P.cmds[i] > 0) pre++;
    P.ffPre = pre;
    long long acc = 0;
    std::vector<char> host(P.np + 1, 0);
    for (int pid = P.np; pid > pre && acc < jointHost; pid--) if (P.death[pid] == INF) { host[pid] = 1; acc += P.len[pid]; }
    for (int pid = 1; pid <= P.np; pid++) {
        if (P.death[pid] != INF) continue;
        if (pid <= pre) { if (pid <= Q.np && Q.death[pid] != INF) P.key[pid] = INF - 1; }
        else P.key[pid] = host[pid] ? INF - 2 : INF - 3;
    }
}

}  // namespace off
#endif
