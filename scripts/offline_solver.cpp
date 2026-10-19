// Оффлайн-решатель (знает всё будущее теста): адреса труб по «первому подходящему по времени-пространству»
// в порядке убывания времени смерти. Нужен для эксперимента «максимум на открытых тестах» (overfit.cpp).
//   offline_solver <test.in> [out.pos]     печатает KPI; out.pos — по одному D на укладку
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <cstdint>
using namespace std;

static const int INF = 1 << 30;

struct Seg {
    int n;
    vector<int> mn, mx, lz;
    void init(int n_) { n = 1; while (n < n_) n <<= 1; mn.assign(2 * n, INF); mx.assign(2 * n, INF); lz.assign(2 * n, INF); }
    inline void apply(int x, int v) { if (v < mn[x]) mn[x] = v; if (v < mx[x]) mx[x] = v; if (v < lz[x]) lz[x] = v; }
    inline void push(int x) { if (lz[x] != INF) { apply(2 * x, lz[x]); apply(2 * x + 1, lz[x]); lz[x] = INF; } }
    void upd(int x, int l, int r, int a, int b, int v) {
        if (b <= l || r <= a) return;
        if (a <= l && r <= b) { apply(x, v); return; }
        push(x);
        int m = (l + r) / 2;
        upd(2 * x, l, m, a, b, v); upd(2 * x + 1, m, r, a, b, v);
        mn[x] = min(mn[2 * x], mn[2 * x + 1]); mx[x] = max(mx[2 * x], mx[2 * x + 1]);
    }
    void chmin(int a, int b, int v) { upd(1, 0, n, a, b, v); }
    // первый индекс >= pos со значением >= d (или n)
    int firstGE(int x, int l, int r, int pos, int d) {
        if (r <= pos || mx[x] < d) return n;
        if (r - l == 1) return l;
        push(x);
        int m = (l + r) / 2;
        int res = firstGE(2 * x, l, m, pos, d);
        if (res != n) return res;
        return firstGE(2 * x + 1, m, r, pos, d);
    }
    // первый индекс >= pos со значением < d (или n)
    int firstLT(int x, int l, int r, int pos, int d) {
        if (r <= pos || mn[x] >= d) return n;
        if (r - l == 1) return l;
        push(x);
        int m = (l + r) / 2;
        int res = firstLT(2 * x, l, m, pos, d);
        if (res != n) return res;
        return firstLT(2 * x + 1, m, r, pos, d);
    }
    int fit(int L, int d) {
        int pos = 0;
        for (;;) {
            int p = firstGE(1, 0, n, pos, d);
            if (p >= n) return -1;
            int q = firstLT(1, 0, n, p, d);
            if (q - p >= L) return p;
            pos = q;
        }
    }
};

int main(int argc, char** argv) {
    ifstream f(argv[1]);
    double w; int n;
    f >> w >> n;
    vector<int> cmds(n);
    for (auto& c : cmds) f >> c;
    vector<int> len(1), birth(1), death(1);
    for (int i = 0; i < n; i++) {
        if (cmds[i] > 0) { len.push_back(cmds[i]); birth.push_back(i); death.push_back(INF); }
        else death[-cmds[i]] = i;
    }
    int np = (int)len.size() - 1;
    // пик живых и UB
    long long cur = 0, peak = 0;
    for (int c : cmds) { if (c > 0) { cur += c; peak = max(peak, cur); } else cur -= len[-c]; }
    long long liveEnd = cur;
    long long ubT = (peak + 19) / 20 * 20;
    vector<int> order(np);
    for (int i = 0; i < np; i++) order[i] = i + 1;
    sort(order.begin(), order.end(), [&](int a, int b) {
        if (death[a] != death[b]) return death[a] > death[b];
        if (death[a] == INF) return birth[a] < birth[b];   // выжившие — по порядку рождения, плотно
        return len[a] > len[b];
    });
    // решение для заданного порядка; возвращает maxEnd, заполняет D, перечень «вылезших» труб (конец > target)
    long long excess = 0;
    auto solve = [&](const vector<int>& ord, vector<int>& Dout, long long target, vector<int>& bad) {
        excess = 0;
        Seg seg;
        seg.init((int)(ubT * 1.3) + 4000);
        long long maxEnd = 0;
        bad.clear();
        for (int p : ord) {
            int a = seg.fit(len[p], death[p]);
            Dout[p] = a;
            seg.chmin(a, a + len[p], birth[p]);
            maxEnd = max<long long>(maxEnd, a + len[p]);
            if (a + len[p] > target) { bad.push_back(p); excess += a + len[p] - target; }
        }
        return maxEnd;
    };
    vector<int> D(np + 1), bad, bestD;
    long long target = ubT;
    long long bestEnd = solve(order, D, target, bad);
    bestD = D;
    unsigned long long rs = 88172645463325252ull;
    auto rnd = [&]() { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; };
    int iters = argc > 3 ? atoi(argv[3]) : 0;
    long long curEx = excess, bestEx = excess;
    vector<int> curOrder = order, bestOrder = order;
    vector<int> curBad = bad;
    for (int it = 0; it < iters && bestEx > 0; it++) {
        order = curOrder;
        vector<int> pos(np + 1);
        for (int i = 0; i < np; i++) pos[order[i]] = i;
        int moves = 1 + (int)(rnd() % 3);
        for (int mv = 0; mv < moves && !curBad.empty(); mv++) {
            int bp = curBad[rnd() % curBad.size()];
            int p = bp;
            if (rnd() % 3 == 0) {   // вместо самой трубы двигаем случайную из пересекающихся с ней по времени
                for (int tries = 0; tries < 50; tries++) {
                    int q = 1 + (int)(rnd() % np);
                    if (birth[q] < death[bp] && death[q] > birth[bp]) { p = q; break; }
                }
            }
            int cur = pos[p];
            int span = 20 + (int)(rnd() % 3000);
            int to;
            if (rnd() % 4 == 0) to = min(np - 1, cur + 1 + (int)(rnd() % (unsigned)(1 + min(np - 1 - cur, span))));
            else to = max(0, cur - 1 - (int)(rnd() % (unsigned)(1 + min(cur, span))));
            int q = order[cur];
            order.erase(order.begin() + cur);
            order.insert(order.begin() + to, q);
            for (int i = min(to, cur); i <= max(to, cur); i++) pos[order[i]] = i;
        }
        long long e = solve(order, D, target, bad);
        if (excess <= curEx) {
            curEx = excess; curOrder = order; curBad = bad;
            if (excess < bestEx || (excess == bestEx && e < bestEnd)) { bestEx = excess; bestEnd = e; bestD = D; fprintf(stderr, "it %d maxEnd=%lld excess=%lld bad=%zu\n", it, e, excess, bad.size()); }
        }
    }
    if (getenv("DBG")) {
        solve(bestOrder.empty() ? order : curOrder, D, target, bad);
        for (int p : bad) fprintf(stderr, "BAD pid=%d birth=%d death=%d len=%d D=%d\n", p, birth[p], death[p] == INF ? -1 : death[p], len[p], D[p]);
    }
    D = bestD;
    long long maxEnd = bestEnd;
    long long trench = (maxEnd + 19) / 20 * 20;
    printf("%s n=%d peak=%lld ubT=%lld maxEnd=%lld trench=%lld kpi=%.5f ub=%.5f\n", argv[1], n, peak, ubT, maxEnd, trench, (double)liveEnd / trench, (double)liveEnd / ubT);
    if (argc > 2) {
        FILE* o = fopen(argv[2], "w");
        for (int p = 1; p <= np; p++) fprintf(o, "%d\n", D[p]);
        fclose(o);
    }
}
