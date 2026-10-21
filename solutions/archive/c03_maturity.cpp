// scorefit_guards — scorefit.cpp + обе страховки от «раздувания» траншеи (2026-10-21):
//   (1) предохранитель зон из scorefit_zguard.cpp — ZGUARD_FREE, inflated(), две строки в zoneTry;
//   (2) страховка через тень best fit из scorefit_guard.cpp — GUARD_*, описание ниже.
// Сборка: c++ -std=c++17 -O2 -o build/bin/scorefit_guards solutions/scorefit_guards.cpp
// -DZGUARD_ON=0 / -DGUARD_ON=0 выключают каждую часть (обе выключены — ответы побитово = scorefit.cpp).
// Описание страховки (2) из scorefit_guard.cpp:
//
// ПОЧЕМУ ОТДЕЛЬНЫЙ ФАЙЛ. Исходный solutions/scorefit.cpp оставлен без изменений как эталон: страховка проверена только
// на оффлайн-наборах (DEV 744 / HOLD 375 тестов) и на открытых тестах не меняет НИ ОДНОГО ответа, т. е. её польза на
// закрытых тестах — предположение. Если она там вредит (ложное срабатывание: наш алгоритм временно отстаёт от best fit,
// а потом обгоняет — так было на two_phase_real при порогах 0.5–2 %), вернуться к scorefit.cpp. Подробности и цифры —
// docs/variants/scorefit.md, раздел «Страховка через тень best fit», и build/hypB/01/REPORT.md.
//
// Что делает. Параллельно ведётся тень — копия траншеи, где каждая труба кладётся чистым best fit (самый короткий
// подходящий промежуток, иначе хвост с минимальной докопкой); тень получает те же укладки и удаления. Если наша
// траншея длиннее тени больше чем на 10 % И больше чем на 1000 м (абсолютный порог — чтобы в крошечной траншее разница
// в 2 блока по 20 м не включала страховку), дальше решаем как best fit на своём состоянии, пока разрыв не станет < 5 %.
// Срабатывает только на «катастрофах», где эвристики scorefit (зоны, потоки, линейная оценка) проигрывают best fit
// в разы: size_shift/08 0.526 → 0.911, slivers/02 0.905 → 0.958. На остальных тестах ответы побитово = scorefit.cpp.
// Цена: память тени (пик 50 МБ против 25 МБ на stress/rand), CPU +~5–15 %. Тень выключается при T > SHADOW_MAX_T.
// Порт на Cangjie (solutions/cangjie/scorefit.cj) страховку НЕ содержит.
//
// scorefit — решение задачи «Gas pipeline construction» (вариант feni, см. docs/variants/scorefit.md).
//
// Сборка: c++ -std=c++17 -O2 -o build/bin/scorefit solutions/scorefit.cpp
// Эксперименты (все ручки и оффлайн-симулятор): solutions/scorefit_lab.cpp — ответы этой версии с ним совпадают.
//
// Стратегия:
//   1. Кандидаты: по одному промежутку каждой длины (самый левый) + самый короткий из больших; хвост траншеи
//      (промежуток, упирающийся в конец) — отдельно; докопка минимальная и только если ничто не подходит.
//   2. Вне потоков кандидат оценивается линейной функцией признаков (gapCost, веса LW): остаток r, точное
//      попадание, «мёртвый» остаток (короче минимальной недавней длины), доля длин r и s среди недавних укладок
//      (не ломать дыру ходового размера, оставлять остаток ходового размера), риск соседей по возрасту и др.
//      Веса подобраны CMA-ES на сгенерированных классах и проверены на независимых наборах
//      (docs/variants/scorefit.md). Если лучший кандидат дороже LW[13], а в хвосте место есть без докопки — хвост.
//   3. Остаток в промежутке оставляем рядом с соседом, который скорее умрёт (и остаток сольётся): обычно с более
//      молодым, а если по статистике удалений умирают скорее старые (очередь) — с более старым.
//   4. Потоки. Онлайн ищем признак, носители одного значения которого умирают вместе (каппа совпадения
//      признака у соседних удалений пачки). Признаки: чётность индекса в серии равных длин, чётность PID,
//      «индекс в серии mod 3 = k», после длинных серий — ещё mod 4 и mod 5 («удаляется каждая k-я»).
//      Если нашёлся — поток A кладём к левому краю «граничного» промежутка, поток B растёт влево от правой
//      границы роста F (участками по 400 м). Какая бы половина ни умерла, освобождается сплошной кусок.
//      Участки по очереди зеркальны ([A→ ←B][B→ ←A]…): области одного потока соседних участков смежны и при
//      вымирании сливаются. Недозаполненный участок в конце теста — чистая потеря (≈ C/2), поэтому участок не
//      больше max(200, 10·L, экономия) (в длинной серии равных длин — до +80 м, чтобы щель между потоками была
//      минимальной: sorted_groups +0.4…0.7 %, 9 +0.9 %), где экономия — насколько меньше свободного внутри траншеи у нас, чем у
//      «тени» без потоков (её копия создаётся при первой возможной метке потока): pairs/periodic +0.5/+0.2 %.
//      Если подтвердилось, что B вымирает почти целиком, а A — нет (разница долей смертей > 0.9), зона B
//      отделяется: B — только в дыры между B (иначе к границе роста), A не занимает дыры между двумя B.
//      В потоках — best fit и штраф «замуровывания» дыры рядом с трубой, чья длина сейчас часто удаляется.
//   5. В крошечной траншее (< 60 м) — точное попадание, иначе самая большая дыра.
//   6. Зоны по длинам (zoneTry). Генератор открытых тестов удаляет почти всегда одну из последних труб ТОЙ ЖЕ
//      длины, а старые короткие во 2-й фазе 15 — «ходоком» вниз по стеку длины. Пока 4–8 длин примерно одинаково
//      часты, у каждой ≥ 150 живых и удаления действительно идут «по стеку длины», трубы частой длины кладутся
//      в свои участки по порядку (стек свободных слотов): умершие подряд сливаются в сплошные куски для длинных.
//      15: 0.953 → 0.961. Редкие длины в это время — «точно / остаток ≥ m самого редкого размера / хвост».
//   7. Признак LW25: изменение избытка дыр над ожидаемым спросом (держать разнообразие размеров дыр, чтобы
//      следующие трубы находили точные дыры).
//   8. Сторона остатка — по «гонке соседей» (какой из соседей промежутка умирает первым), если статистика уверенная.
//   9. «Мёртвый» остаток — короче минимальной ожидаемой длины m; m не меньше самой короткой длины за последние
//      1000 укладок (после смены режима медленная статистика отставала: two_phase_real 10 тыс. команд +0.7 %).
//
// Код написан с оглядкой на перенос на Cangjie: без std::set/map/priority_queue и лямбд — кучи и
// постраничные массивы сделаны вручную, все структуры — массивы.
#include <cstdio>
#include <cmath>
#include <vector>
#include <unistd.h>

typedef long long ll;

static const int NB = 128;     // промежутки короче NB — в корзинах по точной длине, длиннее — в общей куче
static const int MAXL = 40;    // максимальная длина трубы
static const int CHUNK = 400, CHUNK_SAME = 800;   // участок под поток B, м (CHUNK_SAME — при группах одного размера)
static const int CHUNK_MIN = 200, CHUNK_PER_L = 10;   // нижняя граница участка, пока потоки не доказали экономию
static const int CHUNK_LOW = 40, CHUNK_LOW_AGE = 1000;   // нижняя граница, если тень долго не видит экономии
static const int EARLY_RUN = 64, EARLY_MAXL = 3, EARLY_CHUNK = 40;   // потоки по чётности до первых удалений (см. chooseStream)
static const int FIT_EXTRA = 4, FIT_RUN = 30;
static const int MULTI_K = 3; static const ll MULTI_T = 50000;   // кандидатов каждой длины в большой траншее   // подгонка участка под длину серии (см. fitChunk)
// Страховка (см. шапку): тень чистого best fit; если наш T > T_тени·(1 + GUARD_PCT/100) и T − T_тени ≥ GUARD_ABS
// (не раньше GUARD_MIN укладок), решаем как best fit, пока разрыв не станет < GUARD_PCT/200 (гистерезис).
// GUARD_ON=0 — ответы побитово как у scorefit.cpp.
#ifndef GUARD_ON
#define GUARD_ON 1
#endif
#ifndef GUARD_PCT
#define GUARD_PCT 10.0
#endif
#ifndef GUARD_ABS
#define GUARD_ABS 1000
#endif
#ifndef GUARD_MIN
#define GUARD_MIN 2000
#endif
static const ll SHADOW_MAX_T = 3000000;   // тень удваивает память; в траншее длиннее участок C/2 ничтожен — тень удаляем
// Зоны по длинам (см. zoneTry): слотов в участке, мин. доля длины среди недавних укладок, макс. число таких длин,
// макс. доля самой частой длины, мин. число живых труб длины для своего участка.
static const int ZSLOTS = 16;
static const double ZSHARE = 0.08;
static const int ZMAXACT = 8;
static const double ZMAXSHARE = 0.25;
static const int ZMINALIVE = 150;
static const double ZFAST = 0.7;     // смена режима: доля длин с зонами среди недавних укладок (затухание 0.97) ниже — выкл.
static const double EXWARM = 300;  // признак LW25 — только после прогрева статистики длин (cntTot)
static const double SHIFT_TV = 0.7; static const int SHIFT_GAP = 500;   // сброс долей длин для признаков остатка
static const int MSEEN = 1000;     // окно (укладок) для минимальной ожидаемой длины, см. minRecent
#ifndef ZGUARD_ON
#define ZGUARD_ON 1
#endif
static const double ZGUARD_FREE = 0.15;   // предохранитель зон: доля свободного в траншее, выше которой зоны не копают
static const double ZLIFO = 0.7;     // зоны — только если удаления «по стеку длины»: удаляемая среди 8 последних той же длины

// Веса линейной оценки кандидата (вне потоков). Подобраны CMA-ES на сгенерированных классах (two_phase,
// alternating, lifo_mix, phases, sorted_groups, periodic; scripts/cmaes_weights.py) и проверены на независимых наборах.
// Признаки — см. gapCost: 0 остаток r, 1 точное попадание, 2 мёртвый остаток, 3 r·[мёртвый], 4 живой остаток,
// 5 доля длины r среди недавних укладок, 6 то же для длины s, 7 риск соседа со стороны остатка, 8 r·[мёртвый]·
// P(оба соседа выживут), 9 точное попадание: max интенсивность смертей соседей, 10 log(1+число дыр длины s),
// 11 log(s/L), 12 мёртвый остаток ≥ 4, 13 цена хвоста без докопки, 15 точное: min риск соседей,
// 16 живой остаток < 2m, 17 штраф «замуровывания», 18/19 дефицит размера s / r, 25 изменение избытка дыр над
// ожидаемым спросом за LW[26] укладок (+r, −s). 20–24 не используются (24 «хвост с докопкой» ломал del_half_group).
// LW17 поднят вручную 8.7 → 50: в обучении CMA не было сценария 10 («удалили 3 рядом с дырой 1 → кладут 4»),
// на 10 +1.2 %, на остальных 1700 тестах ±0.02 %.
static const double LW[27] = {1, -38.8, -22.4, 13.5, 47, -203, 107, -7.25, 2.78, 2.26, -8.96, 23.4, 1.05e+03, 69, 0, 46.7, 27.9, 50, 0, 0, 0, 0, 0, 0, 0, 20, 50};

static inline double dmax(double a, double b) { return a > b ? a : b; }
static inline double dmin(double a, double b) { return a < b ? a : b; }
static inline double dabs(double a) { return a < 0 ? -a : a; }
static inline int imax(int a, int b) { return a > b ? a : b; }

// ---------------------------------------------------------------- постраничный массив int
// Память растёт вместе с траншеей без копирования (худший случай — длина 8·10^6 м).
struct Paged {
    static const int SH = 16, MASK = (1 << SH) - 1;
    std::vector<int*> pg;
    std::vector<int> nz;                 // ленивый режим: число ненулевых значений на странице
    // Ленивый режим: страница создаётся при первой ненулевой записи и освобождается, когда снова вся нулевая
    // (у tag ненулевые только концы промежутков, у zl — участки зон): плотно занятая траншея их не держит.
    bool lazy;
    explicit Paged(bool lz = false) : lazy(lz) {}
    Paged(const Paged& o) : nz(o.nz), lazy(o.lazy) {   // глубокая копия (тень аллокатора)
        for (size_t k = 0; k < o.pg.size(); k++) {
            int* p = 0;
            if (o.pg[k]) { p = new int[1 << SH]; for (int i = 0; i <= MASK; i++) p[i] = o.pg[k][i]; }
            pg.push_back(p);
        }
    }
    Paged& operator=(const Paged&) = delete;
    ~Paged() { for (size_t i = 0; i < pg.size(); i++) delete[] pg[i]; for (size_t i = 0; i < spare.size(); i++) delete[] spare[i]; }
    std::vector<int*> spare;             // освобождённые страницы (уже нулевые) — без повторного обнуления
    int* newPage() {
        if (!spare.empty()) { int* p = spare.back(); spare.pop_back(); return p; }
        int* p = new int[1 << SH]; for (int i = 0; i <= MASK; i++) p[i] = 0; return p;
    }
    void grow(ll n) {
        while ((ll)pg.size() << SH < n) { pg.push_back(lazy ? 0 : newPage()); if (lazy) nz.push_back(0); }
    }
    int get(ll i) const { const int* p = pg[i >> SH]; return p ? p[i & MASK] : 0; }
    void set(ll i, int v) {
        size_t k = (size_t)(i >> SH);
        int*& p = pg[k];
        if (!lazy) { p[i & MASK] = v; return; }
        if (!p) { if (v == 0) return; p = newPage(); }
        int& c = p[i & MASK];
        nz[k] += (v != 0) - (c != 0);
        c = v;
        if (nz[k] == 0) { if (spare.size() < 4) spare.push_back(p); else delete[] p; p = 0; }
    }
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
#ifndef C3_ORACLE
#define C3_ORACLE 0
#endif
#if C3_ORACLE
extern std::vector<int> ORD; extern long long CUR; extern double K[64];
#endif
struct Alloc {
    ll T = 0;                 // длина траншеи
    // tag[x] > 0 — в x начинается свободный промежуток длины tag[x]; tag[x] < 0 — в x кончается промежуток,
    // начавшийся в -tag[x]-1; 0 — ни то ни другое (соседние промежутки всегда склеены). owner[x] — PID трубы
    // в клетке x (в свободных клетках — устаревшее значение, туда не смотрим).
    Paged tag{true}, owner, zl{true};               // zl[x] — длина-владелец зарезервированного участка (зоны), 0 — нет
    // Корзины: bucket[s] — начала промежутков длины s (куча по адресу); big — ключ len·2^32 + start.
    // Удаление ленивое: элемент действителен, только если tag[start] всё ещё равен его длине.
    Heap bucket[NB], big;
    int gapCnt[NB];                      // точное число промежутков длины s < NB

    std::vector<int> pos, len;           // по PID
    std::vector<char> alive;
    // признаки при рождении: bit0 — чётность индекса в серии равных длин, bit1 — чётность PID,
    // bit2/3/4 — индекс в серии mod 3 равен 0/1/2, bit5..8 — mod 4 равен 0..3, bit9..13 — mod 5 равен 0..4
    // (mod 4/5 — только после хотя бы одной серии ≥ FIT_RUN: в коротких сериях periodic совпадения случайны)
    static const int NF = 14;
    std::vector<unsigned short> feat;
    std::vector<signed char> stream;     // -1 — вне потоков, 0 — A, 1 — B

    // статистика длин укладок (с затуханием) — для «минимальной ожидаемой длины»
    double cnt[MAXL + 1], cntTot = 0;
    int lastSeen[MAXL + 1];              // PID последней трубы каждой длины (0 — не было)
    double cntR[MAXL + 1]; int lastShift = 0;   // доли длин для признаков остатка (см. finishPlace)
    // доля смертей по длине за всё время (для крошечной траншеи)
    double birthL[MAXL + 1], deathL[MAXL + 1];
    // текущая интенсивность смертей по длине: недавние удаления (с затуханием) и число живых
    double recDel[MAXL + 1], recDelTot = 0;
    int aliveCntL[MAXL + 1], aliveTot = 0;

    // потоки
    int F = 0;                           // правая граница роста потока B; F == T — обычный режим
    bool flip = false;                   // в текущем участке B растёт слева, A — справа (участки зеркальны)
    // тень без потоков: та же история команд, своя раскладка (до первой возможной метки совпадает с нами)
    bool noStreams = false, shadowDropped = false;
    Alloc* shadow = 0; ll shadowBorn = 0;
    // страховка (GUARD_ON): pureBF — сама является тенью чистого best fit (без потоков и зон)
    bool pureBF = false;
    Alloc* gsh = 0;                      // тень чистого best fit на той же истории команд
    bool gDropped = false, gOn = false;  // gOn — основной сейчас решает как best fit
    ll occ = 0;                          // занято метров
    int runIdx = 0, lastL = 0;
    int runPrev = 0, runLast = 0;        // длины двух последних серий равных длин (≥ FIT_RUN)
    double coSame[NF], coTot[NF];        // соседние удаления в пачке: совпал ли признак (с затуханием)
    double birthV[NF][2], deathV[NF][2]; // рождения и смерти по значению признака (с затуханием)
    int curFeat = -1;
    double bornLab[2], diedLab[2];       // проверка прогноза: рождения и смерти среди помеченных A/B
    // порядок смертей: средняя доля живых, которые моложе удаляемой трубы (0 — умирают свежие, 1 — старые), с затуханием
    double qSum = 0, qTot = 0;
    // «гонка соседей»: у промежутков, куда клали с остатком, какой из двух соседей (старший/младший) умер первым
    // (затухание 0.99). Записи — списки по PID: raceHead[p] -> элемент e (запись e/2), raceNext[e], raceOther[e].
    double raceOld = 0, raceTot = 0;
    std::vector<int> raceHead, raceNext, raceOther;
    std::vector<char> raceDone;
    int lastDelPid = 0; bool lastWasDel = false;

    Alloc() {
        if (logN[1] == 0) initLogs();
        pos.push_back(-1); len.push_back(0); feat.push_back(0); stream.push_back(-1); alive.push_back(0);
        for (int l = 0; l <= MAXL; l++) {
            cnt[l] = cntF[l] = birthL[l] = deathL[l] = recDel[l] = 0; aliveCntL[l] = 0; lastPtr[l] = 0; lastSeen[l] = 0; cntR[l] = 0;
            for (int k = 0; k < 8; k++) lastOfLen[l][k] = 0;
        }
        bornLab[0] = bornLab[1] = diedLab[0] = diedLab[1] = 0;
        for (int i = 0; i < NB; i++) gapCnt[i] = 0;
        for (int i = 0; i < 24; i++) deathAge[i] = deathSuf[i] = 0;
        aliveBit.assign(BITN, 0);
        raceHead.assign(1024, -1);
        for (int f = 0; f < NF; f++) {
            coSame[f] = coTot[f] = 0;
            for (int v = 0; v < 2; v++) birthV[f][v] = deathV[f][v] = 0;
        }
        ensure(1024);
    }
    ~Alloc() { delete shadow; delete gsh; }

    // ---- промежутки
    void ensure(ll n) { tag.grow(n + 1); owner.grow(n + 1); zl.grow(n + 1); }
    void addGap(int s, int l) {
        tag.set(s, l); tag.set(s + l, -s - 1);
        if (l < NB) { bucket[l].push(s); gapCnt[l]++; } else big.push(((ll)l << 32) | s);
    }
    void delGap(int s) { int l = tag.get(s); if (l < NB) gapCnt[l]--; tag.set(s, 0); tag.set(s + l, 0); }
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
        if (gapCnt[s] == 0) return -1;   // счётчик точный: действительных промежутков нет
        Heap& h = bucket[s];
        int keep[3], nk = 0, res = -1;
        while (!h.empty() && nk < 3) {
            int g = (int)h.top();
            if (tag.get(g) != s || (nk > 0 && keep[nk - 1] == g)) { h.pop(); continue; }
            if (g != ts && g + s != F) { res = g; break; }   // найденный остаётся в куче
            h.pop();
            keep[nk++] = g;
        }
        for (int i = 0; i < nk; i++) h.push(keep[i]);
        return res;
    }
    // До k действительных промежутков длины s (по адресу), кроме хвоста и граничного; элементы кучи — на место.
    int someGaps(int s, int ts, int k, int* out) {
        if (gapCnt[s] == 0) return 0;
        Heap& h = bucket[s];
        int keep[16], nk = 0, no = 0;
        while (!h.empty() && nk < 16 && no < k) {
            int g = (int)h.top(); h.pop();
            if (tag.get(g) != s || (nk > 0 && keep[nk - 1] == g)) continue;
            keep[nk++] = g;
            if (g != ts && g + s != F) out[no++] = g;
        }
        for (int i = 0; i < nk; i++) h.push(keep[i]);
        return no;
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
    // Минимальная ожидаемая длина: самая короткая, доля которой в недавних укладках > 0.2 %, но не меньше самой
    // короткой длины среди последних MSEEN укладок. Медленная статистика после смены режима (2-я фаза 11/15:
    // только длины ≥ 8) отстаёт на ~4000 укладок — всё это время мёртвые остатки 1–7 м считались живыми.
    int minRecent() const {
        int a = 1;
        for (int l = 1; l <= MAXL; l++) if (cnt[l] > 0.002 * cntTot) { a = l; break; }
        int np = (int)pos.size() - 1;
        for (int l = 1; l <= MAXL; l++) if (lastSeen[l] > 0 && np - lastSeen[l] < MSEEN) { if (l > a) a = l; break; }
        return a;
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
    double demand(int z) const { return z > MAXL ? 1e9 : (cntTot > 0 ? cnt[z] / cntTot * LW[26] : 0); }

    // риск по возрасту: P(труба когда-нибудь умрёт | дожила до своего возраста). Возраст — в укладках;
    // deathAge[b] — смерти в возрасте [2^b - 1, 2^(b+1) - 1), живые с возрастом ≥ 2^b - 1 — из дерева Фенвика по PID.
    static const int BITN = 1 << 18;     // > 200 000 укладок
    std::vector<int> aliveBit;
    double deathAge[24], deathSuf[24];   // deathSuf[b] — сумма deathAge[b..23] (целые, сумма точная)
    void bitAdd(int i, int v) { for (; i < BITN; i += i & -i) aliveBit[i] += v; }
    int bitSum(int i) const { int r = 0; for (; i > 0; i -= i & -i) r += aliveBit[i]; return r; }
    static int ageBucket(int a) { int b = 0; while ((2 << b) - 1 <= a && b < 23) b++; return b; }
    double hazard(int p) const {
        if (p <= 0) return 0;
        int cur = (int)pos.size() - 1;
        int b = ageBucket(cur - p);
        double dGE = deathSuf[b];
        int lo = cur - ((1 << b) - 1);
        double aGE = lo > 0 ? bitSum(lo) : 0;
        return dGE + aGE > 0 ? dGE / (dGE + aGE) : 0.5;
    }

    // ---- потоки
    double kappa(int f) const {
        if (coTot[f] < 20) return 0;
        double b0 = birthV[f][0], b1 = birthV[f][1], bt = b0 + b1;
        if (bt <= 0) return 0;
        if (f >= 5 && runLast < FIT_RUN) return 0;
        double minSh = f < 5 ? 0.3 : (f < 9 ? 0.2 : 0.16);   // доля редкого значения: 1/k признака mod k
        if ((b0 < b1 ? b0 : b1) < minSh * bt) return 0;   // признак почти константа — каппа вырождена
        double e = (b0 * b0 + b1 * b1) / (bt * bt);     // ожидаемая доля совпадений при независимости
        double o = coSame[f] / coTot[f];
        return e >= 1 ? 0 : (o - e) / (1 - e);
    }
    double pdie(int f, int v) const { return (deathV[f][v] + 1) / (birthV[f][v] + 2); }
    bool anyDel = false;                 // в тесте уже было удаление
    int chooseStream(int pid) {
        if (noStreams) return -1;
        // До первых удалений длинная серия коротких равных длин (как первая группа 9: 1000 труб длины 1, потом
        // удаляется каждая вторая) — сразу двумя потоками по чётности индекса: какая бы половина ни умерла,
        // освободится сплошной кусок (иначе — 500 мёртвых дыр по 1 м). sorted_groups +0.2…0.4 %, 9 +0.5 %.
        if (!anyDel && runIdx >= EARLY_RUN && lastL <= EARLY_MAXL) return runIdx & 1;
        // признак с каппой ≥ 0.6; при нескольких — где сильнее разница p_die; текущий меняем только с запасом
        int bf = -1; double bs = -1, cs = -1;
        for (int f = 0; f < NF; f++) {
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

    bool verB = false;                   // зона B отделена (см. place)
    bool zoneBlocked(int gs, int s, int st) const {
        if (!verB || st < 0) return false;
        int lp = pidLeftOf(gs), rp = pidRightOf(gs + s);
        if (st == 1) return (lp > 0 && stream[lp] != 1) || (rp > 0 && stream[rp] != 1);
        return lp > 0 && rp > 0 && stream[lp] == 1 && stream[rp] == 1;
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
        return sideNb(lp, rp) == lp ? gs + r : gs;   // остаток — к соседу, который скорее умрёт
    }
    // Сосед, к которому ляжет остаток (вне потоков): обычно умирают свежие — к более молодому (больший PID);
    // если умирают скорее старые (очередь) — к более старому, но не к краю траншеи (PID 0).
    // Если гонка соседей уверенно (доля вне 0.3..0.7) говорит, кто из соседей промежутка умирает первым, — по ней:
    // при смеси очереди и стека средний ранг удаляемых (qSum) ошибается (queue_mix: 0.934 → 0.959).
    bool oldDieFirst() const {
        if (raceTot >= 20 && dabs(raceOld / raceTot - 0.5) > 0.2) return raceOld > 0.5 * raceTot;
        return qTot >= 20 && qSum / qTot > 0.5;
    }
    void raceAdd(int a, int b) {   // a, b — живые соседи промежутка
        if ((int)raceHead.size() < (int)pos.size() + 1) raceHead.resize(pos.size() + 1024, -1);
        int id = (int)raceDone.size();
        raceDone.push_back(0);
        raceNext.push_back(raceHead[a]); raceOther.push_back(b); raceHead[a] = 2 * id;
        raceNext.push_back(raceHead[b]); raceOther.push_back(a); raceHead[b] = 2 * id + 1;
    }
    void raceOnDeath(int p) {
        if (p >= (int)raceHead.size()) return;
        for (int e = raceHead[p]; e >= 0; e = raceNext[e]) {
            if (raceDone[e >> 1]) continue;
            raceDone[e >> 1] = 1;
            raceOld = raceOld * 0.99 + (p < raceOther[e] ? 1 : 0);
            raceTot = raceTot * 0.99 + 1;
        }
        raceHead[p] = -1;
    }
    int sideNb(int lp, int rp) const {
        if (oldDieFirst()) { if (lp == 0) return rp; if (rp == 0) return lp; return lp < rp ? lp : rp; }
        return lp > rp ? lp : rp;
    }

    // таблицы логарифмов (те же значения, что std::log, без вызова в горячем цикле)
    static const int LOGN = 4096;
    static double logN[LOGN], logSL[NB][MAXL + 1];
    static void initLogs() {
        for (int n = 0; n < LOGN; n++) logN[n] = std::log(1.0 + n);
        for (int s = 0; s < NB; s++) for (int l = 1; l <= MAXL; l++) logSL[s][l] = s > 0 ? std::log((double)s / l) : 0;
    }
    static double log1n(int n) { return n < LOGN ? logN[n] : std::log(1.0 + n); }
    // верхняя граница разброса gapCost между промежутками одной длины (признаки соседей: LW7, 8, 9, 15, 17)
    static double nbW(int r) {
        return dabs(LW[7]) + dabs(LW[8]) * r + dabs(LW[9]) + dabs(LW[15]) + dabs(LW[17]) + 1;
    }
    // стоимость укладки трубы L в промежуток gs длины s (меньше — лучше)
    double gapCost(int gs, int s, int L, int st, int m, double selfHz) const {
        if (T < 60) {   // крошечная траншея
            if (s == L) return -1000 + 2000 * dmax(pdieLen(pidLeftOf(gs)), pdieLen(pidRightOf(gs + s)));
            return -s;
        }
        int r = s - L;
        int lp = pidLeftOf(gs), rp = pidRightOf(gs + s);
        // точное попадание «замуровывает» дыру. Плохо, если сосед скоро умрёт (длина его класса сейчас часто
        // удаляется), новая труба живёт дольше, а слитая дыра s + len(соседа) — частая длина среди укладок.
        double seal = 0;
        if (r == 0) {
            double h = 0;
            int nb[2] = {lp, rp};
            for (int k = 0; k < 2; k++) {
                if (nb[k] <= 0) continue;
                int mg = s + len[nb[k]];
                // слитая дыра — частая длина среди недавних укладок (или среди последних ~30: новый режим)
                if (mg > MAXL || (cnt[mg] < 0.05 * cntTot && cntF[mg] < 0.2 * cntFTot)) continue;
                h = dmax(h, hzRate(nb[k]));
            }
            h -= selfHz;
            if (h >= 0.5) seal = h;
        }
        if (st >= 0) return r != 0 ? r : 3 * seal;   // в потоках — best fit
        // вне потоков — линейная оценка (веса LW)
        double hl = hazard(lp), hr = hazard(rp);
        double hside = sideNb(lp, rp) == lp ? hl : hr;   // сторона остатка — как в chooseSide
        bool dead = r > 0 && r < m;
        double dr = (r >= 1 && r <= MAXL && cntTot > 0) ? cntR[r] / cntTot : 0;
        double ds = (s <= MAXL && cntTot > 0) ? cntR[s] / cntTot : 0;
        double ns = s < NB ? gapCnt[s] : 0, nr = r < NB ? gapCnt[r] : 0;
        double c = LW[0] * r + LW[1] * (r == 0) + LW[2] * dead + LW[3] * (dead ? r : 0) + LW[4] * (r >= m)
                 + LW[5] * dr + LW[6] * ds + LW[7] * (r > 0 ? hside : 0) + LW[8] * (dead ? r * (1 - dmax(hl, hr)) : 0)
                 + LW[9] * (r == 0 ? dmax(hzRate(lp), hzRate(rp)) : 0) + LW[10] * log1n((int)ns)
                 + LW[11] * (s < NB ? logSL[s][L] : std::log((double)s / L)) + LW[12] * (dead && r >= 4) + LW[15] * (r == 0 ? dmin(hl, hr) : 0)
                 + LW[16] * (r > 0 && r < 2 * m && !dead);
        if (r > 0) c += LW[18] * ds / (1.0 + ns) + LW[19] * dr / (1.0 + nr);
        if (r > 0 && cntTot >= EXWARM) {   // изменение избытка дыр над ожидаемым спросом: появится дыра r, исчезнет дыра s
            double exr = nr + 1 - demand(r), exs = ns - demand(s);
            c += LW[25] * (std::log(1.0 + dmax(0.0, exr)) - std::log(1.0 + dmax(0.0, exs)));
        }
        c += LW[17] * seal;
#if C3_ORACLE
        if (r > 0 && r < m) {   // оракул: сосед умрёт в ближайшие K команд -> дыра подрастёт
            int nb[2] = {lp, rp};
            for (int k = 0; k < 2; k++) if (nb[k] > 0 && alive[nb[k]] && ORD[nb[k]] - CUR <= (long long)K[1]) { c += K[2]; break; }
        }
#endif
        return c;
    }


    // ---- зоны по длинам. Генератор открытых тестов удаляет почти всегда одну из последних труб той же длины,
    // а старые короткие во 2-й фазе 15 — «ходоком», который ползёт вниз по стеку длины почти подряд. Если трубы
    // частой длины лежат в своём участке по порядку, умершие подряд сливаются в сплошные куски, которые потом
    // заполняют длинные (15: 0.953 → 0.962). Включено, только пока несколько длин примерно одинаково часты и у
    // каждой большая популяция (иначе резерв под колебания численности дороже выигрыша: 8/12/13 не затронуты).
    std::vector<int> slots[MAXL + 1];    // свободные слоты участков длины L (стек)
    ll resEnd = 0;                       // конец последнего зарезервированного участка
    bool zonesEnabled() const {
        if (cntTot < 50 || T < 60) return false;
        int k = 0; double mx = 0;
        for (int l = 1; l <= MAXL; l++) { double sh = cnt[l] / cntTot; if (sh > ZSHARE) k++; mx = dmax(mx, sh); }
        if (!(k >= 4 && k <= ZMAXACT && mx <= ZMAXSHARE)) return false;
        if (lifoTot < 20 || lifoHit < ZLIFO * lifoTot) return false;
        double f = 0;
        for (int l = 1; l <= MAXL; l++) if (cnt[l] > ZSHARE * cntTot) f += cntF[l];
        return f >= ZFAST * cntFTot;
    }
    double cntF[MAXL + 1], cntFTot = 0;  // статистика длин с быстрым затуханием (0.97) — смена режима
    int lastOfLen[MAXL + 1][8], lastPtr[MAXL + 1];   // 8 последних PID каждой длины
    double lifoHit = 0, lifoTot = 0;     // доля удалений, попавших в эти 8 (затухание 0.99)
    bool zoneOn(int z, bool en) const { return en && z > 0 && cnt[z] > ZSHARE * cntTot; }
    bool occupied(ll x) const { int p = owner.get(x); return p > 0 && alive[p] && pos[p] <= x && x < pos[p] + len[p]; }
    bool clearOf(int d, int n, int L, bool en) const {   // в [d, d+n) нет клеток чужих (не L) активных зон
        for (int x = d; x < d + n; x++) { int z = zl.get(x); if (z && z != L && zoneOn(z, en)) return false; }
        return true;
    }
    // Самый левый действительный промежуток длины s, для которого ok(g); элементы кучи возвращаются на место.
    // mode: 0 — точное попадание в свою/незарезервированную зону, 1 — «редкая длина» (не хвост, не чужая зона).
    std::vector<int> tmpKeep;
    int scanGap(int s, int L, int mode, int ts, bool en) {
        Heap& h = bucket[s];
        tmpKeep.clear();
        int res = -1;
        while (!h.empty() && tmpKeep.size() < 64) {   // не больше 64 промежутков (защита от долгого перебора)
            int g = (int)h.top(); h.pop();
            if (tag.get(g) != s || (!tmpKeep.empty() && tmpKeep.back() == g)) continue;
            tmpKeep.push_back(g);
            bool ok = mode == 0 ? (zl.get(g) == L || !zoneOn(zl.get(g), en))
                                : (g != ts && !zoneOn(zl.get(g), en) && clearOf(g, mode == 2 ? s : L, L, en));
            if (ok) { res = g; break; }
        }
        for (size_t i = 0; i < tmpKeep.size(); i++) h.push(tmpKeep[i]);
        return res;
    }
    int scanBig(int L, int ts, bool en, int& outLen, int need) {   // need — сколько клеток проверять на чужие зоны
        std::vector<ll> keep;
        int res = -1;
        while (!big.empty() && keep.size() < 64) {
            ll k = big.top(); big.pop();
            int l = (int)(k >> 32), g = (int)(k & 0xffffffffLL);
            if (tag.get(g) != l || (!keep.empty() && keep.back() == k)) continue;
            keep.push_back(k);
            if (g != ts && !zoneOn(zl.get(g), en) && clearOf(g, need < l ? need : l, L, en)) { res = g; outLen = l; break; }
        }
        for (size_t i = 0; i < keep.size(); i++) big.push(keep[i]);
        return res;
    }
    // остаток >= m самого редкого размера среди незарезервированных дыр (не хвост); false — нет такой
    bool rarePlace(int L, int m, int ts, bool en, int& D) {
        double bc = 1e18; int bg = -1;
        for (int s = L + m; s < NB; s++) {
            int g = scanGap(s, L, 1, ts, en);
            if (g < 0) continue;
            double c = (double)gapCnt[s - L] * 1000 + s;
            if (c < bc) { bc = c; bg = g; }
        }
        { int bl = 0, g = scanBig(L, ts, en, bl, L); if (g >= 0) { double c = bl; if (c < bc) { bc = c; bg = g; } } }
        if (bg >= 0) { D = bg; takeFrom(bg, bg, L); return true; }
        return false;
    }
    // Предохранитель зон (из scorefit_zguard.cpp): траншея «раздута» — свободного внутри больше ZGUARD_FREE длины.
    // Тогда zoneTry не докапывает траншею (ни под новый участок, ни под редкую длину), а отдаёт трубу обычной
    // логике (линейная оценка по дырам). В нормальном режиме свободно 4–8 % — правило не срабатывает.
    bool inflated() const { return ZGUARD_ON && T - occ > ZGUARD_FREE * T; }
    bool zoneTry(int L, int st, int m, ll& N, int& D) {
        if (st >= 0) return false;
        bool en = zonesEnabled();
        if (!en) return false;
        int ts = tailStart();
        N = 0; D = -1;
        if (L < NB) { int g = scanGap(L, L, 0, ts, en); if (g >= 0) { D = g; takeFrom(g, g, L); return true; } }
        if (zoneOn(L, en) && aliveCntL[L] >= ZMINALIVE) {
            // свободный слот своего участка (последний освободившийся — первым)
            while (!slots[L].empty()) {
                int d = slots[L].back(); slots[L].pop_back();
                if (zl.get(d) != L || occupied(d) || occupied(d + L - 1)) continue;
                ll hs = d; int steps = 0;
                while (hs > 0 && !occupied(hs - 1) && steps < 4096) { hs--; steps++; }
                if (steps >= 4096 || tag.get(hs) <= 0 || hs + tag.get(hs) < d + L) continue;
                D = d; takeFrom((int)hs, D, L); return true;
            }
            // участок внутри незарезервированной дыры (>= 4 слотов, ужимается под дыру); иначе незарезервированные
            // дыры по правилу редкой длины; и лишь потом новый участок в хвосте (иначе, когда короткие возвращаются
            // после длинных — three_phase, rev_two_phase, — их дыры простаивали бы)
            {
                int hl = 0, hg = -1;
                for (int s2 = 4 * L; s2 < NB && hg < 0; s2++) hg = scanGap(s2, L, 2, ts, en);
                if (hg >= 0) hl = tag.get(hg);
                else hg = scanBig(L, ts, en, hl, ZSLOTS * L);
                if (hg >= 0) {
                    ll C2 = (ll)ZSLOTS * L, fit = (ll)(hl / L) * L;
                    if (fit < C2) C2 = fit;
                    for (ll x = hg; x < hg + C2; x++) zl.set(x, L);
                    for (ll x = hg + C2 - L; x > hg; x -= L) slots[L].push_back((int)x);
                    D = hg; takeFrom(hg, D, L);
                    return true;
                }
                if (rarePlace(L, m, ts, en, D)) return true;
            }
            // новый участок из ZSLOTS слотов в конце траншеи, за всеми зарезервированными
            ll C = (ll)ZSLOTS * L;
            ll stp = ts > resEnd ? ts : resEnd;
            if (T < stp + C && inflated()) return false;   // предохранитель: раздутая траншея — не копать под зону
            if (T < stp + C) { N = (stp + C - T + 19) / 20; dig(N); }
            for (ll x = stp; x < stp + C; x++) zl.set(x, L);
            for (ll x = stp + C - L; x > stp; x -= L) slots[L].push_back((int)x);
            resEnd = stp + C;
            D = (int)stp; takeFrom(tailStart(), D, L);
            return true;
        }
        // редкая длина: незарезервированные дыры по правилу редкой длины, иначе хвост
        if (rarePlace(L, m, ts, en, D)) return true;
        ll stp = ts > resEnd ? ts : resEnd;
        if (T < stp + L && inflated()) return false;   // предохранитель: раздутая траншея — не копать под редкую длину
        if (T < stp + L) { N = (stp + L - T + 19) / 20; dig(N); }
        D = (int)stp; takeFrom(tailStart(), D, L);
        if (resEnd < stp + L && resEnd > ts) resEnd = stp + L;
        return true;
    }

    // ---- тень без потоков
    // свободное внутри траншеи, кроме хвоста и граничного промежутка потоков
    ll freeIn() const {
        ll f = T - occ - (T - tailStart());
        if (F < T && F > 0 && tag.get(F) < 0) f -= F - (-tag.get(F) - 1);
        return f;
    }
    void shadowCheck() {
        if (shadowDropped || noStreams) return;
        if (shadow) {
            if (T > SHADOW_MAX_T) { delete shadow; shadow = 0; shadowDropped = true; }
            return;
        }
        bool any = false;
        for (int f = 0; f < NF; f++) if (kappa(f) >= 0.6) any = true;   // иначе chooseStream точно вернёт -1
        if (!any) return;
        shadowBorn = (ll)pos.size();
        shadow = new Alloc(*this);
        shadow->shadow = 0; shadow->noStreams = true; shadow->gsh = 0;
    }

    // Число сегментов докопки участка потоков. Участок заполняют с двух сторон трубы одной длины, и щель
    // (tl + 20·N) mod L между потоками остаётся навсегда, если дальше идут более длинные (9: участок 400 м при
    // L = 32 — по 16 м на участок). В длинной серии равных длин берём N из [N, N + FIT_EXTRA] с наименьшей щелью.
    ll fitChunk(int tl, ll N, int L) const {
        if (runIdx < FIT_RUN) return N;
        ll best = N; int br = (int)((tl + 20 * N) % L);
        for (ll n = N + 1; n <= N + FIT_EXTRA && br > 0; n++) { int r = (int)((tl + 20 * n) % L); if (r < br) { br = r; best = n; } }
        return best;
    }
    // страховка: решать как основной, если его T больше тени чистого best fit (см. GUARD_ON)
    bool guardOn() {
        if (gDropped || !gsh) return false;
        if (T > SHADOW_MAX_T) { delete gsh; gsh = 0; gDropped = true; gOn = false; return false; }
        if ((ll)pos.size() - 1 < GUARD_MIN) return false;
        double x = GUARD_PCT / 100.0;
        if (!gOn && (double)T > (double)gsh->T * (1 + x) && T - gsh->T >= GUARD_ABS) gOn = true;
        else if (gOn && (double)T < (double)gsh->T * (1 + x / 2)) gOn = false;
        return gOn;
    }
    // чистый best fit: самый короткий подходящий промежуток (среди равных — левый); иначе хвост с минимальной докопкой
    void placePure(int L, ll& N, int& D) {
        if (F < T) F = (int)T;   // потоки в этом режиме не работают: граничный промежуток становится обычным
        int ts = tailStart();
        int bgs = -1, bs = 0;
        for (int s = L; s < NB && bgs < 0; s++) { int g = firstGap(s, ts); if (g >= 0) { bgs = g; bs = s; } }
        if (bgs < 0) { int bl = 0; int g = firstBig(ts, bl); if (g >= 0) { bgs = g; bs = bl; } }
        N = 0;
        if (bgs >= 0) {
            if (bs > L) { int lp = pidLeftOf(bgs), rp = pidRightOf(bgs + bs); if (lp > 0 && rp > 0) raceAdd(lp, rp); }
            D = chooseSide(bgs, bs, L, -1);
            takeFrom(bgs, D, L);
        } else {
            int tl = (int)(T - ts);
            if (tl < L) { N = (L - tl + 19) / 20; dig(N); ts = tailStart(); }
            D = ts;
            takeFrom(ts, D, L);
        }
    }
    // ---- команды
    void place(int L, ll& outN, ll& outD) {
        shadowCheck();
        if (shadow) { ll sn, sd; shadow->place(L, sn, sd); }
        if (GUARD_ON && !pureBF && !gDropped && !gsh && pos.size() == 1) { gsh = new Alloc(); gsh->pureBF = true; gsh->noStreams = true; }
        if (gsh) { ll sn, sd; gsh->place(L, sn, sd); }
        occ += L;
        int m = minRecent();
        int pid = (int)pos.size();
        if (L == lastL) runIdx++;
        else {
            if (runIdx + 1 >= FIT_RUN) { runPrev = runLast; runLast = runIdx + 1; }
            runIdx = 0;
        }
        lastL = L;
        feat.push_back((unsigned short)((runIdx & 1) | ((pid & 1) << 1) | (1 << (2 + runIdx % 3))
                                        | (1 << (5 + runIdx % 4)) | (1 << (9 + runIdx % 5))));
        bool bf = pureBF || (GUARD_ON && guardOn());
        int st = bf ? -1 : chooseStream(pid);
        stream.push_back((signed char)st);
        if (st >= 0) bornLab[st] += 1;
        for (int f = 0; f < NF; f++) {
            birthV[f][0] *= 0.999; birthV[f][1] *= 0.999;
            deathV[f][0] *= 0.999; deathV[f][1] *= 0.999;
            birthV[f][(feat[pid] >> f) & 1] += 1;
        }
        if (F > T) F = (int)T;
        int ts = tailStart();
        double selfHz = recDelTot >= 1 ? rateOfLen(L) : 0;
        ll N = 0; int D = -1;
        if (bf) { placePure(L, N, D); finishPlace(pid, L, N, D, outN, outD); return; }
        if (zoneTry(L, st, m, N, D)) { finishPlace(pid, L, N, D, outN, outD); return; }

        // прогноз «B вымирает почти целиком, A — нет» подтвердился на помеченных трубах: зона B отделяется —
        // B кладём только в дыры между B (иначе к границе роста), A не занимает дыры между двумя живыми B
        verB = bornLab[1] >= 20 && bornLab[0] >= 20 && diedLab[1] / bornLab[1] - diedLab[0] / bornLab[0] > 0.9;

        double bestC = 1e18; int bgs = -1, bs = 0;
        // в большой траншее вне потоков — до MULTI_K промежутков каждой длины: оценка различает их по риску соседей
        // (two_phase_real +0.03 %, 15 +0.05 %); в маленькой — только самый левый (там это шум: alternating −0.03 %)
        bool multi = st < 0 && T >= MULTI_T;
        for (int s = L; s < NB; s++) {
            int g = firstGap(s, ts);
            if (g < 0) continue;
            double c = zoneBlocked(g, s, st) ? 1e17 : gapCost(g, s, L, st, m, selfHz);
            if (c < bestC) { bestC = c; bgs = g; bs = s; }
            // промежутки той же длины отличаются только признаками соседей (их вклад ≤ nbW): если и с ним первый
            // не лучше лучшего, остальные смотреть незачем (ответы те же, CPU меньше)
            if (multi && c < 1e17 && c - nbW(s - L) < bestC) {
                int gk[MULTI_K]; int ng = someGaps(s, ts, MULTI_K, gk);
                for (int i = 1; i < ng; i++) {
                    double c2 = zoneBlocked(gk[i], s, st) ? 1e17 : gapCost(gk[i], s, L, st, m, selfHz);
                    if (c2 < bestC) { bestC = c2; bgs = gk[i]; bs = s; }
                }
            }
        }
        {
            int bl = 0, g = firstBig(ts, bl);
            if (g >= 0) {
                double c = zoneBlocked(g, bl, st) ? 1e17 : gapCost(g, bl, L, st, m, selfHz);
                if (c < bestC) { bestC = c; bgs = g; bs = bl; }
            }
        }

        bool useGap = bgs >= 0 && bestC < 1e17;
        // вне потоков: если лучший промежуток дороже LW[13], а в хвосте место есть без докопки — в хвост
        if (useGap && st < 0 && T - tailStart() >= L && LW[13] < bestC) useGap = false;

        if (useGap) {
            if (bs > L) { int lp = pidLeftOf(bgs), rp = pidRightOf(bgs + bs); if (lp > 0 && rp > 0) raceAdd(lp, rp); }
            D = chooseSide(bgs, bs, L, st);
            takeFrom(bgs, D, L);
        } else {
            // граничный промежуток — кончается в F (в обычном режиме F == T, это хвост)
            int fgs = (F > 0 && tag.get(F) < 0) ? -tag.get(F) - 1 : F;
            int fl = F - fgs;
            // участок: не больше доказанной тенью экономии потоков (но не меньше max(200, 10·L))
            int chunk = CHUNK;
            // группы одного размера (две последние серии равных длин отличаются ≤ 2 %, как в 5/9): участок может быть
            // больше (CHUNK_SAME), но не больше, чем осталось до конца текущей группы (ниже)
            bool sameGroups = runPrev > 0 && runIdx >= FIT_RUN && 50 * dabs(runPrev - runLast) <= imax(runPrev, runLast);
            if (sameGroups) chunk = CHUNK_SAME;
            // до первых удалений (ранние потоки, см. chooseStream) — маленький участок: если удалений так и не будет
            // (14: 15 000 единиц подряд), пустая середина участка — не больше EARLY_CHUNK
            if (!anyDel) chunk = EARLY_CHUNK;
            if (shadow) {
                ll c = CHUNK_MIN, sv = shadow->freeIn() - freeIn();
                if ((ll)CHUNK_PER_L * L > c) c = (ll)CHUNK_PER_L * L;
                // тень уже CHUNK_LOW_AGE укладок не подтверждает экономию — нижняя граница участка 40 + L
                // (недозаполненный участок в конце теста — чистая потеря: pairs +0.6…0.8 %)
                if ((ll)pos.size() - shadowBorn > CHUNK_LOW_AGE && sv < c && CHUNK_LOW + L < c) c = CHUNK_LOW + L;
                if (sv > c) c = sv;
                if (c < chunk) chunk = (int)c;
                // остаток текущей группы (но ≥ 4·L) — последний участок в конце теста не пустует
                if (sameGroups) {
                    ll rem = (ll)(runLast - runIdx) * L;
                    if (rem < 4 * L) rem = 4 * L;
                    if (rem < chunk) chunk = (int)rem;
                }
            }
            bool far = st >= 0 && ((st == 1) != flip);   // поток растёт от правого края участка (F) влево
            if (fl >= L) {
                D = far ? F - L : fgs;
                takeFrom(fgs, D, L);
                if (far) F = D;
            } else if (F == T) {
                int need = L - fl;
                if (far && chunk - fl > need) need = chunk - fl;
                N = (need + 19) / 20;
                if (far) N = fitChunk(fl, N, L);
                dig(N);
                D = far ? (int)T - L : fgs;
                takeFrom(fgs, D, L);
                F = far ? D : (int)T;
            } else {
                // граничный промежуток исчерпан — новый участок в конце траншеи, зеркальный предыдущему
                if (st >= 0) { flip = !flip; far = (st == 1) != flip; }
                ts = tailStart();
                int tl = (int)(T - ts);
                int need = L - tl;
                if (st >= 0 && chunk - tl > need) need = chunk - tl;
                N = need > 0 ? (need + 19) / 20 : 0;
                if (st >= 0 && N > 0) N = fitChunk(tl, N, L);
                dig(N);
                ts = tailStart();
                D = far ? (int)T - L : ts;
                takeFrom(ts, D, L);
                F = far ? D : (int)T;
            }
        }
        finishPlace(pid, L, N, D, outN, outD);
    }
    void finishPlace(int pid, int L, ll N, int D, ll& outN, ll& outD) {
        for (int i = 0; i < L; i++) owner.set(D + i, pid);
        pos.push_back(D); len.push_back(L); alive.push_back(1);
        bitAdd(pid, 1);
        birthL[L] += 1; aliveCntL[L]++; aliveTot++;
        lastOfLen[L][lastPtr[L]] = pid; lastPtr[L] = (lastPtr[L] + 1) & 7; lastSeen[L] = pid;
        for (int l = 1; l <= MAXL; l++) { cnt[l] *= 0.999; cntF[l] *= 0.97; cntR[l] *= 0.999; }
        cntTot = cntTot * 0.999 + 1; cnt[L] += 1;
        cntFTot = cntFTot * 0.97 + 1; cntF[L] += 1;
        // доли длин для признаков остатка: как cnt, но при резкой смене распределения (полусумма |разниц| с быстрой
        // статистикой > SHIFT_TV) сбрасываются на быструю (не чаще раза в SHIFT_GAP укладок): lifo_mix +0.05 %
        cntR[L] += 1;
        if (cntFTot >= 25 && pid - lastShift > SHIFT_GAP) {
            double tv = 0;
            for (int l = 1; l <= MAXL; l++) tv += dabs(cntF[l] / cntFTot - cntR[l] / cntTot);
            if (tv / 2 > SHIFT_TV) { for (int l = 1; l <= MAXL; l++) cntR[l] = cntF[l] / cntFTot * cntTot; lastShift = pid; }
        }
        lastWasDel = false;
        outN = N; outD = D;
    }

    ll remove(int pid) {
        anyDel = true;
        int d = pos[pid], l = len[pid];
        occ -= l;
        if (shadow) shadow->remove(pid);
        if (gsh) gsh->remove(pid);
        deathL[l] += 1;
        if (aliveTot > 1) {
            double q = (double)(aliveTot - bitSum(pid)) / (aliveTot - 1);
            qSum = qSum * 0.99 + q; qTot = qTot * 0.99 + 1;
        }
        bitAdd(pid, -1); alive[pid] = 0;
        raceOnDeath(pid);
        { int b = ageBucket((int)pos.size() - 1 - pid); deathAge[b] += 1; for (int k = 0; k <= b; k++) deathSuf[k] += 1; }
        for (int k = 1; k <= MAXL; k++) recDel[k] *= 0.98;
        recDelTot = recDelTot * 0.98 + 1; recDel[l] += 1; aliveCntL[l]--; aliveTot--;
        bool hit = false;
        for (int k = 0; k < 8; k++) if (lastOfLen[l][k] == pid) hit = true;
        lifoHit = lifoHit * 0.99 + (hit ? 1 : 0); lifoTot = lifoTot * 0.99 + 1;
        if (zl.get(d) == l) slots[l].push_back(d);
        int e = freeRange(d, l);
        if (d == F) F = e;
        if (stream[pid] >= 0) diedLab[stream[pid]] += 1;
        for (int f = 0; f < NF; f++) {
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

double Alloc::logN[Alloc::LOGN], Alloc::logSL[NB][MAXL + 1];

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
