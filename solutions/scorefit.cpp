// scorefit — решение задачи «Gas pipeline construction» (вариант feni, см. docs/variants/scorefit.md).
//
// Сборка: c++ -std=c++17 -O2 -o build/bin/scorefit solutions/scorefit.cpp
// Эксперименты (все ручки и оффлайн-симулятор): solutions/scorefit_lab.cpp — ответы этой версии с ним совпадают.
//
// Стратегия:
//   1. Best fit: из свободных промежутков, куда помещается труба, берём самый короткий (среди равных — левый).
//      Хвост траншеи (промежуток, упирающийся в конец) — в последнюю очередь; докопка минимальная.
//   2. Остаток в промежутке оставляем рядом с более молодым соседом (он скорее умрёт, и остаток сольётся).
//   3. Потоки. Онлайн ищем признак, носители одного значения которого умирают вместе (каппа совпадения
//      признака у соседних удалений пачки). Признаки: чётность индекса в серии равных длин, чётность PID.
//      Если нашёлся — поток A кладём к левому краю «граничного» промежутка, поток B растёт влево от правой
//      границы роста F (участками по 400 м). Какая бы половина ни умерла, освобождается сплошной кусок.
//      Если прогноз «B умирает чаще» подтвердился на помеченных трубах, B идёт сразу к границе роста.
//   4. Мелкие правила: не «замуровывать» точным попаданием дыру рядом с трубой, чья длина сейчас часто
//      удаляется (если слитая дыра — частая длина); мёртвый остаток (< минимальной недавней длины) ≥ 5 м — лучше в хвост (вне потоков);
//      в крошечной траншее (< 60 м) — точное попадание, иначе самая большая дыра.
//
// Код написан с оглядкой на перенос на Cangjie: без std::set/map/priority_queue и лямбд — кучи и
// постраничные массивы сделаны вручную, все структуры — массивы.
#include <cstdio>
#include <vector>
#include <unistd.h>

typedef long long ll;

static const int NB = 128;     // промежутки короче NB — в корзинах по точной длине, длиннее — в общей куче
static const int MAXL = 40;    // максимальная длина трубы
static const int CHUNK = 400;  // участок, который докапывается под поток B, м

static inline double dmax(double a, double b) { return a > b ? a : b; }
static inline double dmin(double a, double b) { return a < b ? a : b; }
static inline double dabs(double a) { return a < 0 ? -a : a; }
static inline int imax(int a, int b) { return a > b ? a : b; }

// ---------------------------------------------------------------- постраничный массив int
// Память растёт вместе с траншеей без копирования (худший случай — длина 8·10^6 м).
struct Paged {
    static const int SH = 16, MASK = (1 << SH) - 1;
    std::vector<int*> pg;
    ~Paged() { for (size_t i = 0; i < pg.size(); i++) delete[] pg[i]; }
    void grow(ll n) {
        while ((ll)pg.size() << SH < n) {
            int* p = new int[1 << SH];
            for (int i = 0; i <= MASK; i++) p[i] = 0;
            pg.push_back(p);
        }
    }
    int get(ll i) const { return pg[i >> SH][i & MASK]; }
    void set(ll i, int v) { pg[i >> SH][i & MASK] = v; }
};

// ---------------------------------------------------------------- двоичная куча (минимум) на массиве
struct Heap {
    std::vector<ll> a;
    bool empty() const { return a.empty(); }
    ll top() const { return a[0]; }
    void push(ll v) {
        a.push_back(v);
        size_t i = a.size() - 1;
        while (i > 0) {
            size_t p = (i - 1) / 2;
            if (a[p] <= a[i]) break;
            ll t = a[p]; a[p] = a[i]; a[i] = t;
            i = p;
        }
    }
    void pop() {
        a[0] = a.back();
        a.pop_back();
        size_t n = a.size(), i = 0;
        while (true) {
            size_t l = 2 * i + 1, r = l + 1, m = i;
            if (l < n && a[l] < a[m]) m = l;
            if (r < n && a[r] < a[m]) m = r;
            if (m == i) break;
            ll t = a[m]; a[m] = a[i]; a[i] = t;
            i = m;
        }
    }
};

// ---------------------------------------------------------------- аллокатор
struct Alloc {
    ll T = 0;                 // длина траншеи
    // tag[x] > 0 — в x начинается свободный промежуток длины tag[x]; tag[x] < 0 — в x кончается промежуток,
    // начавшийся в -tag[x]-1; 0 — ни то ни другое (соседние промежутки всегда склеены). owner[x] — PID трубы
    // в клетке x (в свободных клетках — устаревшее значение, туда не смотрим).
    Paged tag, owner;
    // Корзины: bucket[s] — начала промежутков длины s (куча по адресу); big — ключ len·2^32 + start.
    // Удаление ленивое: элемент действителен, только если tag[start] всё ещё равен его длине.
    Heap bucket[NB], big;

    std::vector<int> pos, len;           // по PID
    std::vector<unsigned char> feat;     // bit0 — чётность индекса в серии равных длин, bit1 — чётность PID
    std::vector<signed char> stream;     // -1 — вне потоков, 0 — A, 1 — B

    // статистика длин укладок (с затуханием) — для «минимальной ожидаемой длины»
    double cnt[MAXL + 1], cntTot = 0;
    // доля смертей по длине за всё время (для крошечной траншеи)
    double birthL[MAXL + 1], deathL[MAXL + 1];
    // текущая интенсивность смертей по длине: недавние удаления (с затуханием) и число живых
    double recDel[MAXL + 1], recDelTot = 0;
    int aliveCntL[MAXL + 1], aliveTot = 0;

    // потоки
    int F = 0;                           // правая граница роста потока B; F == T — обычный режим
    int runIdx = 0, lastL = 0;
    double coSame[2], coTot[2];          // соседние удаления в пачке: совпал ли признак (с затуханием)
    double birthV[2][2], deathV[2][2];   // рождения и смерти по значению признака (с затуханием)
    int curFeat = -1;
    double bornLab[2], diedLab[2];       // проверка прогноза: рождения и смерти среди помеченных A/B
    int lastDelPid = 0; bool lastWasDel = false;

    Alloc() {
        pos.push_back(-1); len.push_back(0); feat.push_back(0); stream.push_back(-1);
        for (int l = 0; l <= MAXL; l++) { cnt[l] = birthL[l] = deathL[l] = recDel[l] = 0; aliveCntL[l] = 0; }
        for (int f = 0; f < 2; f++) {
            coSame[f] = coTot[f] = 0; bornLab[f] = diedLab[f] = 0;
            for (int v = 0; v < 2; v++) birthV[f][v] = deathV[f][v] = 0;
        }
        ensure(1024);
    }

    // ---- промежутки
    void ensure(ll n) { tag.grow(n + 1); owner.grow(n + 1); }
    void addGap(int s, int l) {
        tag.set(s, l); tag.set(s + l, -s - 1);
        if (l < NB) bucket[l].push(s); else big.push(((ll)l << 32) | s);
    }
    void delGap(int s) { int l = tag.get(s); tag.set(s, 0); tag.set(s + l, 0); }
    int freeRange(int s, int l) {   // освободить [s, s+l) со склейкой соседей; вернуть конец склеенного
        int e = s + l;
        if (e < T && tag.get(e) > 0) { int r = tag.get(e); delGap(e); e += r; }
        if (tag.get(s) < 0) { int ps = -tag.get(s) - 1; delGap(ps); s = ps; }
        addGap(s, e - s);
        return e;
    }
    void takeFrom(int gs, int d, int l) {   // занять [d, d+l) внутри промежутка, начинающегося в gs
        int gl = tag.get(gs);
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
    int tailStart() const { return (T > 0 && tag.get(T) < 0) ? -tag.get(T) - 1 : (int)T; }
    int pidLeftOf(int x) const { return x > 0 ? owner.get(x - 1) : 0; }   // x — начало промежутка
    int pidRightOf(int x) const { return x < T ? owner.get(x) : 0; }      // x — конец промежутка

    // Первый действительный промежуток длины s (по адресу), кроме хвоста ts и граничного (кончается в F).
    // -1 — нет. Устаревшие элементы и дубликаты выбрасываются.
    int firstGap(int s, int ts) {
        Heap& h = bucket[s];
        int keep[3], nk = 0, res = -1;
        while (!h.empty() && nk < 3) {
            int g = (int)h.top();
            if (tag.get(g) != s || (nk > 0 && keep[nk - 1] == g)) { h.pop(); continue; }
            h.pop();
            keep[nk++] = g;
            if (g != ts && g + s != F) { res = g; break; }
        }
        for (int i = 0; i < nk; i++) h.push(keep[i]);
        return res;
    }
    // То же для больших промежутков: самый короткий (среди равных — левый). Возвращает start, длину — в outLen.
    int firstBig(int ts, int& outLen) {
        ll keep[3]; int nk = 0, res = -1;
        while (!big.empty() && nk < 3) {
            ll k = big.top();
            int l = (int)(k >> 32), g = (int)(k & 0xffffffffLL);
            if (tag.get(g) != l || (nk > 0 && keep[nk - 1] == k)) { big.pop(); continue; }
            big.pop();
            keep[nk++] = k;
            if (g != ts && g + l != F) { res = g; outLen = l; break; }
        }
        for (int i = 0; i < nk; i++) big.push(keep[i]);
        return res;
    }

    // ---- статистика
    int minRecent() const {   // минимальная длина, доля которой в недавних укладках > 0.2 %
        for (int l = 1; l <= MAXL; l++) if (cnt[l] > 0.002 * cntTot) return l;
        return 1;
    }
    double pdieLen(int p) const {
        if (p <= 0) return 0;
        int l = len[p];
        return (deathL[l] + 0.5) / (birthL[l] + 1);
    }
    double rateOfLen(int l) const {   // (доля недавних удалений длины l) / (доля живых длины l), /4, ≤ 1
        double rel = (recDel[l] / recDelTot) / ((double)imax(1, aliveCntL[l]) / imax(1, aliveTot));
        return dmin(1.0, rel / 4);
    }
    double hzRate(int p) const { return (p <= 0 || recDelTot < 1) ? 0 : rateOfLen(len[p]); }

    // ---- потоки
    double kappa(int f) const {
        if (coTot[f] < 20) return 0;
        double b0 = birthV[f][0], b1 = birthV[f][1], bt = b0 + b1;
        if (bt <= 0) return 0;
        if ((b0 < b1 ? b0 : b1) < 0.3 * bt) return 0;   // признак почти константа — каппа вырождена
        double e = (b0 * b0 + b1 * b1) / (bt * bt);     // ожидаемая доля совпадений при независимости
        double o = coSame[f] / coTot[f];
        return e >= 1 ? 0 : (o - e) / (1 - e);
    }
    double pdie(int f, int v) const { return (deathV[f][v] + 1) / (birthV[f][v] + 2); }
    int chooseStream(int pid) {
        // признак с каппой ≥ 0.6; при нескольких — где сильнее разница p_die; текущий меняем только с запасом
        int bf = -1; double bs = -1, cs = -1;
        for (int f = 0; f < 2; f++) {
            double k = kappa(f);
            if (k < 0.6) continue;
            double s = k + dabs(pdie(f, 0) - pdie(f, 1));
            if (f == curFeat) cs = s;
            if (s > bs) { bs = s; bf = f; }
        }
        if (bf < 0) { curFeat = -1; return -1; }
        if (cs >= 0 && bs < cs + 0.2) bf = curFeat;
        curFeat = bf;
        int v = (feat[pid] >> bf) & 1;
        double pv = pdie(bf, v), po = pdie(bf, v ^ 1);
        if (dabs(pv - po) < 0.1) return v;
        return pv > po ? 1 : 0;           // кто умирает чаще — в поток B
    }

    int chooseSide(int gs, int s, int L, int st) {   // где в промежутке gs длины s положить трубу
        int r = s - L;
        if (r == 0) return gs;
        int lp = pidLeftOf(gs), rp = pidRightOf(gs + s);
        if (st >= 0) {   // к соседу из того же потока; иначе A — слева, B — справа
            int lc = lp ? stream[lp] : -2, rc = rp ? stream[rp] : -2;
            if (lc == st && rc != st) return gs;
            if (rc == st && lc != st) return gs + r;
            return st == 1 ? gs + r : gs;
        }
        return lp > rp ? gs + r : gs;   // остаток — к более молодому соседу
    }

    // стоимость укладки трубы L в промежуток gs длины s (меньше — лучше)
    double gapCost(int gs, int s, int L, double selfHz) const {
        if (T < 60) {   // крошечная траншея
            if (s == L) return -1000 + 2000 * dmax(pdieLen(pidLeftOf(gs)), pdieLen(pidRightOf(gs + s)));
            return -s;
        }
        if (s != L) return s - L;   // best fit
        // точное попадание «замуровывает» дыру. Плохо, если сосед скоро умрёт (длина его класса сейчас часто
        // удаляется), новая труба живёт дольше, а слитая дыра s + len(соседа) — частая длина среди укладок.
        double h = 0;
        int nb[2] = {pidLeftOf(gs), pidRightOf(gs + s)};
        for (int k = 0; k < 2; k++) {
            if (nb[k] <= 0) continue;
            int mg = s + len[nb[k]];
            if (mg > MAXL || cnt[mg] < 0.05 * cntTot) continue;
            h = dmax(h, hzRate(nb[k]));
        }
        h -= selfHz;
        return h >= 0.5 ? 3 * h : 0;
    }

    // ---- команды
    void place(int L, ll& outN, ll& outD) {
        int m = minRecent();
        int pid = (int)pos.size();
        if (L == lastL) runIdx++; else runIdx = 0;
        lastL = L;
        feat.push_back((unsigned char)((runIdx & 1) | ((pid & 1) << 1)));
        int st = chooseStream(pid);
        stream.push_back((signed char)st);
        if (st >= 0) bornLab[st] += 1;
        for (int f = 0; f < 2; f++) {
            birthV[f][0] *= 0.999; birthV[f][1] *= 0.999;
            deathV[f][0] *= 0.999; deathV[f][1] *= 0.999;
            birthV[f][(feat[pid] >> f) & 1] += 1;
        }
        if (F > T) F = (int)T;
        int ts = tailStart();
        double selfHz = recDelTot >= 1 ? rateOfLen(L) : 0;

        double bestC = 1e18; int bgs = -1, bs = 0;
        for (int s = L; s < NB; s++) {
            int g = firstGap(s, ts);
            if (g < 0) continue;
            double c = gapCost(g, s, L, selfHz);
            if (c < bestC) { bestC = c; bgs = g; bs = s; }
        }
        {
            int bl = 0, g = firstBig(ts, bl);
            if (g >= 0) {
                double c = gapCost(g, bl, L, selfHz);
                if (c < bestC) { bestC = c; bgs = g; bs = bl; }
            }
        }

        bool useGap = bgs >= 0;
        // прогноз «B умирает» подтвердился — B сразу к границе роста
        if (st == 1 && bornLab[1] >= 20 && bornLab[0] >= 20 &&
            diedLab[1] / bornLab[1] - diedLab[0] / bornLab[0] > 0.5) useGap = false;
        // мёртвый остаток (короче всех ожидаемых труб) ≥ 5 м — лучше в хвост (вне потоков)
        if (useGap && st < 0) {
            int r = bs - L;
            if (r > 0 && r < m && r >= 5) useGap = false;
        }

        ll N = 0; int D;
        if (useGap) {
            D = chooseSide(bgs, bs, L, st);
            takeFrom(bgs, D, L);
        } else {
            // граничный промежуток — кончается в F (в обычном режиме F == T, это хвост)
            int fgs = (F > 0 && tag.get(F) < 0) ? -tag.get(F) - 1 : F;
            int fl = F - fgs;
            if (fl >= L) {
                D = st == 1 ? F - L : fgs;
                takeFrom(fgs, D, L);
                if (st == 1) F = D;
            } else if (F == T) {
                int need = L - fl;
                if (st == 1 && CHUNK - fl > need) need = CHUNK - fl;
                N = (need + 19) / 20;
                dig(N);
                D = st == 1 ? (int)T - L : fgs;
                takeFrom(fgs, D, L);
                F = st == 1 ? D : (int)T;
            } else {
                // граничный промежуток исчерпан — новый участок в конце траншеи
                ts = tailStart();
                int tl = (int)(T - ts);
                int need = L - tl;
                if (st >= 0 && CHUNK - tl > need) need = CHUNK - tl;
                N = need > 0 ? (need + 19) / 20 : 0;
                dig(N);
                ts = tailStart();
                D = st == 1 ? (int)T - L : ts;
                takeFrom(ts, D, L);
                F = st == 1 ? D : (int)T;
            }
        }
        for (int i = 0; i < L; i++) owner.set(D + i, pid);
        pos.push_back(D); len.push_back(L);
        birthL[L] += 1; aliveCntL[L]++; aliveTot++;
        for (int l = 1; l <= MAXL; l++) cnt[l] *= 0.999;
        cntTot = cntTot * 0.999 + 1; cnt[L] += 1;
        lastWasDel = false;
        outN = N; outD = D;
    }

    ll remove(int pid) {
        int d = pos[pid], l = len[pid];
        deathL[l] += 1;
        for (int k = 1; k <= MAXL; k++) recDel[k] *= 0.98;
        recDelTot = recDelTot * 0.98 + 1; recDel[l] += 1; aliveCntL[l]--; aliveTot--;
        int e = freeRange(d, l);
        if (d == F) F = e;
        if (stream[pid] >= 0) diedLab[stream[pid]] += 1;
        for (int f = 0; f < 2; f++) {
            int v = (feat[pid] >> f) & 1;
            deathV[f][v] += 1;
            if (lastWasDel) {
                coSame[f] *= 0.998; coTot[f] *= 0.998;
                coTot[f] += 1; coSame[f] += (v == ((feat[lastDelPid] >> f) & 1));
            }
        }
        lastWasDel = true; lastDelPid = pid;
        return d;
    }
};

// ---------------------------------------------------------------- ввод-вывод
// Читаем через read(): он возвращает то, что уже пришло, и не ждёт заполнения буфера (интерактив).
static char ibuf[1 << 16];
static int ipos = 0, ilen = 0;
static int gc() {
    if (ipos == ilen) { ilen = (int)read(0, ibuf, sizeof ibuf); ipos = 0; if (ilen <= 0) return -1; }
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
static void writeOut(const char* s, int n) {   // write() без буфера — ответ уходит сразу (аналог flush)
    while (n > 0) { ssize_t k = write(1, s, n); if (k <= 0) return; s += k; n -= (int)k; }
}

int main() {
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
