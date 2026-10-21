// bestfit: классический best fit (эталон для сравнения со scorefit и portfolio).
// Труба длины L кладётся в самый короткий свободный промежуток внутри траншеи, куда она влезает (среди равных — левый,
// труба прижата к левому краю промежутка). Если такого нет — в хвост (свободная часть траншеи за последней трубой),
// докопка минимальная: ceil((L - свободный хвост) / 20) сегментов. При удалении промежутки склеиваются с соседями
// (и с хвостом, если труба была последней). Память: два массива int на длину траншеи, до 8·10^6 клеток (≈ 64 МБ).
// Открытые тесты: Σ ≈ 28.37 из 32.1 (scorefit 30.15). Сборка: c++ -std=c++17 -O2 -o bestfit solutions/bestfit.cpp
#include <cstdio>
#include <set>
#include <utility>
#include <vector>
#ifdef _WIN32
#include <io.h>
static int rd(void* b, unsigned n) { return _read(0, b, n); }
static int wr(const void* b, unsigned n) { return _write(1, b, n); }
#else
#include <unistd.h>
static int rd(void* b, unsigned n) { return (int)read(0, b, n); }
static int wr(const void* b, unsigned n) { return (int)write(1, b, n); }
#endif

typedef long long ll;

struct BestFit {
    ll T = 0;                          // длина траншеи
    ll tail = 0;                       // свободен хвост [tail, T); перед tail — труба или граница (промежутки с хвостом склеены)
    std::set<std::pair<int, int>> gaps;   // свободные промежутки внутри траншеи: (длина, начало)
    std::vector<int> gapLen;           // gapLen[x] — длина промежутка, начинающегося в x (0 — нет)
    std::vector<int> gapBegin;         // gapBegin[x] — (начало + 1) промежутка, кончающегося в x (0 — нет)
    std::vector<int> pos, len;         // по PID (с 1)
    BestFit() { pos.push_back(-1); len.push_back(0); gapLen.assign(21, 0); gapBegin.assign(21, 0); }

    void addGap(int s, int l) { gaps.insert({l, s}); gapLen[s] = l; gapBegin[s + l] = s + 1; }
    void delGap(int s) { int l = gapLen[s]; gaps.erase({l, s}); gapLen[s] = 0; gapBegin[s + l] = 0; }

    void place(int L, ll& N, ll& D) {
        N = 0;
        auto it = gaps.lower_bound({L, -1});
        if (it != gaps.end()) {
            int l = it->first, s = it->second;
            delGap(s);
            if (l > L) addGap(s + L, l - L);
            D = s;
        } else {
            ll freeTail = T - tail;
            if (freeTail < L) {
                N = (L - freeTail + 19) / 20;
                T += 20 * N;
                gapLen.resize(T + 2, 0); gapBegin.resize(T + 2, 0);
            }
            D = tail;
            tail += L;
        }
        pos.push_back((int)D); len.push_back(L);
    }

    ll remove(int pid) {
        int d = pos[pid], l = len[pid];
        int s = d, e = d + l;
        if (gapBegin[s]) { int ps = gapBegin[s] - 1; delGap(ps); s = ps; }      // левый сосед — промежуток
        if (e == tail) { tail = s; return d; }                                   // последняя труба: промежуток уходит в хвост
        if (gapLen[e]) { int r = gapLen[e]; delGap(e); e += r; }                 // правый сосед — промежуток
        addGap(s, e - s);
        return d;
    }
};

static char ibuf[1 << 16];
static int ipos = 0, ilen = 0;
static int gc() {
    if (ipos == ilen) { ilen = rd(ibuf, sizeof ibuf); ipos = 0; if (ilen <= 0) return -1; }
    return ibuf[ipos++];
}
static bool readInt(int& x) {
    int c = gc();
    while (c != '-' && (c < '0' || c > '9')) { if (c == -1) return false; c = gc(); }
    bool neg = false;
    if (c == '-') { neg = true; c = gc(); }
    x = 0;
    while (c >= '0' && c <= '9') { x = x * 10 + (c - '0'); c = gc(); }
    if (neg) x = -x;
    return true;
}
static int putNum(char* b, ll v) {
    char t[24]; int n = 0;
    if (v == 0) t[n++] = '0';
    while (v > 0) { t[n++] = (char)('0' + v % 10); v /= 10; }
    for (int i = 0; i < n; i++) b[i] = t[n - 1 - i];
    return n;
}

int main() {
    BestFit* A = new BestFit();
    int c; char ob[64];
    while (readInt(c)) {
        if (c == 0) return 0;
        int n;
        if (c > 0) {
            ll N, D;
            A->place(c, N, D);
            n = putNum(ob, N); ob[n++] = ' '; n += putNum(ob + n, D);
        } else {
            n = putNum(ob, A->remove(-c));
        }
        ob[n++] = '\n';
        int off = 0;
        while (off < n) { int k = wr(ob + off, n - off); if (k <= 0) return 0; off += k; }   // без буфера — ответ уходит сразу
    }
    return 0;
}
