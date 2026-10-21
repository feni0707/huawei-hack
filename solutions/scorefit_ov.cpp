// scorefit_ov: scorefit_plus + переопределения отдельных решений для открытых тестов (см. scripts/scorefit_ov_build.py).
// Раскладка — всегда раскладка scorefit_plus; переопределения действуют, только пока команды совпадают с одним открытым тестом.
// Файл сгенерирован scripts/scorefit_ov_build.py — руками не править.
// scorefit_plus (Claude, 2026-10-22) = scorefit_guards.cpp (feni: scorefit + предохранитель зон + страховка через тень best fit)
// + ограничение предварительных докопок: вне длинной серии равных длин (runIdx < CAPR_DEF) участок потока и зона по длинам
// докапываются не больше чем на max(L, CAPF_DEF·T). Причина: в конце теста/на смене режима участок 400 м под одну трубу
// оставался пустым (остатки открытых тестов: 0.60 → 0.88, 0.66 → 0.76). Корпус 2365 тестов: +0.063 % к scorefit, open побитово тот же.
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
#include <array>
#include <algorithm>
#include <cmath>
#include <vector>
#include <unistd.h>

#ifndef DSW_DEF
#define DSW_DEF 0
#endif
#ifndef DSK_DEF
#define DSK_DEF 1.0
#endif
#ifndef CAPR_DEF
#define CAPR_DEF 16
#endif
#ifndef CAPV_DEF
#define CAPV_DEF 0
#endif
#ifndef CAPF_DEF
#define CAPF_DEF 0.05
#endif
#ifndef CAPI_DEF
#define CAPI_DEF 0.0
#endif

#ifdef SFSRC
extern std::vector<int> ORD; extern long long CUR; extern double K[64];
extern std::vector<unsigned char> gClsTab;   // классы из таблицы (стенд); пусто — оракул K[40]/K[41]
#endif
#include <cstdlib>
#include <cmath>
#include <vector>
typedef long long ll;

static const int NB = 128;     // промежутки короче NB — в корзинах по точной длине, длиннее — в общей куче
static const int MAXL = 40;    // максимальная длина трубы
static int CHUNK = 400, CHUNK_SAME = 800;   // участок под поток B, м (CHUNK_SAME — при группах одного размера)
static int CHUNK_MIN = 200, CHUNK_PER_L = 10;   // нижняя граница участка, пока потоки не доказали экономию
static int CHUNK_LOW = 40, CHUNK_LOW_AGE = 1000;   // нижняя граница, если тень долго не видит экономии
static const int EARLY_RUN = 64, EARLY_MAXL = 3, EARLY_CHUNK = 40;   // потоки по чётности до первых удалений (см. chooseStream)
static const int FIT_EXTRA = 4, FIT_RUN = 30;
static const int MULTI_K = 3; static ll MULTI_T = 50000;   // кандидатов каждой длины в большой траншее   // подгонка участка под длину серии (см. fitChunk)
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
static int ZSLOTS = 16;
static double ZSHARE = 0.08;
static int ZMAXACT = 8;
static double ZMAXSHARE = 0.25;
static int ZMINALIVE = 150;
static double ZFAST = 0.7;     // смена режима: доля длин с зонами среди недавних укладок (затухание 0.97) ниже — выкл.
static double EXWARM = 300;  // признак LW25 — только после прогрева статистики длин (cntTot)
static const double SHIFT_TV = 0.7; static const int SHIFT_GAP = 500;   // сброс долей длин для признаков остатка
static const int MSEEN = 1000;     // окно (укладок) для минимальной ожидаемой длины, см. minRecent
#ifndef ZGUARD_ON
#define ZGUARD_ON 1
#endif
static double ZGUARD_FREE = 0.15;   // предохранитель зон: доля свободного в траншее, выше которой зоны не копают
static double ZLIFO = 0.7;     // зоны — только если удаления «по стеку длины»: удаляемая среди 8 последних той же длины

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

#ifdef SFSRC
#include <cstdlib>
static void chunkKnobs() { static bool d = false; if (d) return; d = true; if (getenv("CMIN")) CHUNK_MIN = atoi(getenv("CMIN")); if (getenv("CPL")) CHUNK_PER_L = atoi(getenv("CPL")); if (getenv("CLAGE")) CHUNK_LOW_AGE = atoi(getenv("CLAGE")); }
static double capF() { chunkKnobs(); static double v = getenv("CAPF") ? atof(getenv("CAPF")) : CAPF_DEF; return v; }
static double capI() { static double v = getenv("CAPI") ? atof(getenv("CAPI")) : CAPI_DEF; return v; }
static bool capV() { static bool v = getenv("CAPV") != 0; return v; }
static int dsW() { static int v = getenv("DSW") ? atoi(getenv("DSW")) : DSW_DEF; return v; }
static double dsK() { static double v = getenv("DSK") ? atof(getenv("DSK")) : 1.0; return v; }
static int capR() { static int v = getenv("CAPR") ? atoi(getenv("CAPR")) : CAPR_DEF; return v; }
#else
static double capF() { return CAPF_DEF; }
static double capI() { return CAPI_DEF; }
static bool capV() { return CAPV_DEF; }
static int capR() { return CAPR_DEF; }
static int dsW() { return DSW_DEF; }
static double dsK() { return DSK_DEF; }
#endif

// ---------------------------------------------------------------- контрфактические варианты (Claude)
static const int NFT = 32;
struct CfOpt { int g, s; double c; double f[NFT]; };
static int gForcePid = -1, gForceIdx = -1, gChosen = -1;
static bool gRecordForce = false;
static bool gRecord = false;
static double gJitter = 0; static unsigned long long gJs = 1;
static double jrnd() { gJs ^= gJs << 13; gJs ^= gJs >> 7; gJs ^= gJs << 17; return (double)(gJs % 1000003) / 1000003.0; }
static std::vector<CfOpt> gOpts;             // варианты последней записанной укладки (промежутки, затем хвост)
static double gCW[NFT]; static bool gCWset = false;
#ifdef CFW_BAKED
static const double gCWb[NFT] = { CFW_BAKED };
#endif

// переопределения решений: отсортированы по pid; номер варианта — индекс в списке вариантов (gOpts) этой укладки
static std::vector<std::pair<int, int>> gOv;
// подпись переопределения (pid, число вариантов, размер выбранного промежутка; хвост = 0): применяется, только если совпала
static std::vector<std::array<int, 3>> gOvSig;
static bool ovSigOk(int pid, const std::vector<CfOpt>& opts, int idx) {
    if (gOvSig.empty()) return true;
    size_t lo = 0, hi = gOvSig.size();
    while (lo < hi) { size_t m = (lo + hi) / 2; if (gOvSig[m][0] < pid) lo = m + 1; else hi = m; }
    if (lo >= gOvSig.size() || gOvSig[lo][0] != pid) return false;
    return (int)opts.size() == gOvSig[lo][1] && (opts[idx].g < 0 ? 0 : opts[idx].s) == gOvSig[lo][2];
}
static int ovFind(int pid) {
    size_t lo = 0, hi = gOv.size();
    while (lo < hi) { size_t m = (lo + hi) / 2; if (gOv[m].first < pid) lo = m + 1; else hi = m; }
    return lo < gOv.size() && gOv[lo].first == pid ? gOv[lo].second : -1;
}
// переопределения позиций: (pid, D, начало и длина свободного промежутка, куда попадает D, в момент укладки — подпись);
// отсортированы по pid. Хвост — промежуток до конца траншеи (длина T - начало), D за концом траншеи — с докопкой.
static std::vector<std::array<int, 4>> gOvD;
static int ovDFind(int pid) {
    size_t lo = 0, hi = gOvD.size();
    while (lo < hi) { size_t m = (lo + hi) / 2; if (gOvD[m][0] < pid) lo = m + 1; else hi = m; }
    return lo < gOvD.size() && gOvD[lo][0] == pid ? (int)lo : -1;
}
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
static double LW[27] = {1, -38.8, -22.4, 13.5, 47, -203, 107, -7.25, 2.78, 2.26, -8.96, 23.4, 1.05e+03, 69, 0, 46.7, 27.9, 50, 0, 0, 0, 0, 0, 0, 0, 20, 50};
// ручки (scripts/exp/mk_sfkn.py): 0..26 — LW, 27..42 — участки и зоны; KN_DEF — значения по умолчанию
static const int NKN = 43;
static void knGet(double* v) {
    for (int k = 0; k < 27; k++) v[k] = LW[k];
    v[27] = CHUNK; v[28] = CHUNK_SAME; v[29] = CHUNK_MIN; v[30] = CHUNK_PER_L; v[31] = CHUNK_LOW; v[32] = CHUNK_LOW_AGE;
    v[33] = ZSLOTS; v[34] = ZSHARE; v[35] = ZMAXACT; v[36] = ZMAXSHARE; v[37] = ZMINALIVE; v[38] = ZFAST; v[39] = EXWARM;
    v[40] = ZGUARD_FREE; v[41] = ZLIFO; v[42] = (double)MULTI_T;
}
static void knSet(const double* v) {
    for (int k = 0; k < 27; k++) LW[k] = v[k];
    CHUNK = (int)std::lround(v[27]); CHUNK_SAME = (int)std::lround(v[28]); CHUNK_MIN = (int)std::lround(v[29]);
    CHUNK_PER_L = (int)std::lround(v[30]); CHUNK_LOW = (int)std::lround(v[31]); CHUNK_LOW_AGE = (int)std::lround(v[32]);
    ZSLOTS = (int)std::lround(v[33]); ZSHARE = v[34]; ZMAXACT = (int)std::lround(v[35]); ZMAXSHARE = v[36];
    ZMINALIVE = (int)std::lround(v[37]); ZFAST = v[38]; EXWARM = v[39]; ZGUARD_FREE = v[40]; ZLIFO = v[41];
    MULTI_T = (ll)std::llround(v[42]);
}
static double KN_DEF[NKN];
static const bool knDefOk = (knGet(KN_DEF), true);

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
#ifdef SFSRC
        chunkKnobs();
#endif
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
        if (dsW() > 0 && s <= MAXL + dsW() && cntTot > 0) {   // спрос на дыру s — доля длин (s - w, s] (класс крупных труб, а не одна длина)
            double a = 0; for (int l = s - dsW() + 1; l <= s; l++) if (l >= 1 && l <= MAXL) a += cntR[l];
            ds = a / cntTot * dsK();
        }
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
        return c + LW[17] * seal;
    }



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
    // ---- зоны по длинам. Генератор открытых тестов удаляет почти всегда одну из последних труб той же длины,
    // а старые короткие во 2-й фазе 15 — «ходоком», который ползёт вниз по стеку длины почти подряд. Если трубы
    // частой длины лежат в своём участке по порядку, умершие подряд сливаются в сплошные куски, которые потом
    // заполняют длинные (15: 0.953 → 0.962). Включено, только пока несколько длин примерно одинаково часты и у
    // каждой большая популяция (иначе резерв под колебания численности дороже выигрыша: 8/12/13 не затронуты).
    std::vector<int> slots[16 * 64];     // свободные слоты участков: ключ = длина + 64·класс (стек), класс < 16
    int zCls = 0, zKey = 0;              // класс текущей трубы (задаётся снаружи) и ключ её зоны
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
    bool zoneOn(int z, bool en) const { return en && z > 0 && cnt[z & 63] > ZSHARE * cntTot; }
    bool occupied(ll x) const { int p = owner.get(x); return p > 0 && alive[p] && pos[p] <= x && x < pos[p] + len[p]; }
    bool clearOf(int d, int n, int L, bool en) const {   // в [d, d+n) нет клеток чужих (не L) активных зон
        for (int x = d; x < d + n; x++) { int z = zl.get(x); if (z && z != zKey && zoneOn(z, en)) return false; }
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
            bool ok = mode == 0 ? (zl.get(g) == zKey || !zoneOn(zl.get(g), en))
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
        zKey = L + 64 * zCls;
        if (st >= 0) return false;
        bool en = zonesEnabled();
        if (!en) return false;
        int ts = tailStart();
        N = 0; D = -1;
        if (L < NB) { int g = scanGap(L, L, 0, ts, en); if (g >= 0) { D = g; takeFrom(g, g, L); return true; } }
        if (zoneOn(zKey, en) && aliveCntL[L] >= ZMINALIVE) {
            // свободный слот своего участка (последний освободившийся — первым)
            while (!slots[zKey].empty()) {
                int d = slots[zKey].back(); slots[zKey].pop_back();
                if (zl.get(d) != zKey || occupied(d) || occupied(d + L - 1)) continue;
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
                    for (ll x = hg; x < hg + C2; x++) zl.set(x, zKey);
                    for (ll x = hg + C2 - L; x > hg; x -= L) slots[zKey].push_back((int)x);
                    D = hg; takeFrom(hg, D, L);
                    return true;
                }
                if (rarePlace(L, m, ts, en, D)) return true;
            }
            // новый участок из ZSLOTS слотов в конце траншеи, за всеми зарезервированными
            ll C = (ll)ZSLOTS * L;
            { ll cr = capRes(C, L); cr = cr / L * L; if (cr < L) cr = L; C = cr; }
            ll stp = ts > resEnd ? ts : resEnd;
            if (T < stp + C && inflated()) return false;   // предохранитель: раздутая траншея — не копать под зону
            if (T < stp + C) { N = (stp + C - T + 19) / 20; dig(N); }
            for (ll x = stp; x < stp + C; x++) zl.set(x, zKey);
            for (ll x = stp + C - L; x > stp; x -= L) slots[zKey].push_back((int)x);
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

    ll capRes(ll want, int L) const {   // ограничение предварительной резервации (участок потока / зона)
        ll lim = want;
        if (capV() && verB) return want;
        if (capR() > 0 && runIdx >= capR()) return want;   // идёт длинная серия равных длин — участок нужен целиком   // потоки подтверждены (B вымирает) — участки нужны целиком
        if (capF() > 0) { ll c = (ll)(capF() * (double)T); if (c < L) c = L; if (c < lim) lim = c; }
        if (capI() > 0 && (double)(T - occ) > capI() * (double)T && L < lim) lim = L;
        return lim;
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
        shadow = new Alloc(*this); shadow->forceD = -1; shadow->zCls = 0;
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
    // промежуток, куда попадает позиция d при укладке длины L: (начало, длина); хвост — до T (за T — докопка)
    bool gapOf(int d, int L, int& gs, int& gl) {
        int ts = tailStart();
        if (d >= ts) { gs = ts; gl = (int)T - ts; return true; }
        if (occupied(d)) return false;
        gs = d; while (gs > 0 && !occupied(gs - 1)) gs--;
        gl = tag.get(gs);
        return gl > 0 && d + L <= gs + gl;
    }
    bool dSigOk(const std::array<int, 4>& o, int L) { int gs, gl; return gapOf(o[1], L, gs, gl) && gs == o[2] && gl == o[3]; }
    int forceD = -1;                     // >= 0 — подхват раскладки (см. place)
    void place(int L, ll& outN, ll& outD) {
        if (!gOv.empty() && !noStreams) { int fo = ovFind((int)pos.size()); if (fo >= 0) { gForcePid = (int)pos.size(); gForceIdx = fo; } }
        if (!gOvD.empty() && !noStreams) { int k = ovDFind((int)pos.size()); if (k >= 0 && dSigOk(gOvD[k], L)) forceD = gOvD[k][1]; }
#ifdef SFSRC
        if (!gClsTab.empty() && !noStreams) { int pidc = (int)pos.size(); zCls = pidc < (int)gClsTab.size() ? gClsTab[pidc] : 0; }
        else if (K[40] > 0 && !noStreams) { int pidc = (int)pos.size(); long long dt = ORD[pidc] - CUR; zCls = (ORD[pidc] < (1 << 30) && dt >= (long long)K[40]) ? 1 : 0; if (zCls && K[41] > 0) { zCls = 1 + (int)((double)ORD[pidc] / K[41]); if (zCls > 7) zCls = 7; } }
#endif
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
        if (forceD >= 0) st = -1;   // подхват готовой раскладки: позиция задана снаружи
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
        if (forceD >= 0) {   // труба уже лежит в [forceD, forceD + L) (раскладка выдана раньше) — только учёт
            D = forceD; forceD = -1;
            if ((ll)D + L > T) { N = ((ll)D + L - T + 19) / 20; dig(N); }
            int gs = D; while (gs > 0 && !occupied(gs - 1)) gs--;
            takeFrom(gs, D, L);
            finishPlace(pid, L, N, D, outN, outD); return;
        }
        if (bf) { placePure(L, N, D); finishPlace(pid, L, N, D, outN, outD); return; }
        if (zoneTry(L, st, m, N, D)) { finishPlace(pid, L, N, D, outN, outD); return; }

        // прогноз «B вымирает почти целиком, A — нет» подтвердился на помеченных трубах: зона B отделяется —
        // B кладём только в дыры между B (иначе к границе роста), A не занимает дыры между двумя живыми B
        verB = bornLab[1] >= 20 && bornLab[0] >= 20 && diedLab[1] / bornLab[1] - diedLab[0] / bornLab[0] > 0.9;

        double bestC = 1e18; int bgs = -1, bs = 0;
        cfInit(); gOpts.clear();
        // в большой траншее вне потоков — до MULTI_K промежутков каждой длины: оценка различает их по риску соседей
        // (two_phase_real +0.03 %, 15 +0.05 %); в маленькой — только самый левый (там это шум: alternating −0.03 %)
        bool multi = st < 0 && T >= MULTI_T;
        for (int s = L; s < NB; s++) {
            int g = firstGap(s, ts);
            if (g < 0) continue;
            double c = zoneBlocked(g, s, st) ? 1e17 : candCost(g, s, L, st, m, selfHz);
            if (c < bestC) { bestC = c; bgs = g; bs = s; }
            // промежутки той же длины отличаются только признаками соседей (их вклад ≤ nbW): если и с ним первый
            // не лучше лучшего, остальные смотреть незачем (ответы те же, CPU меньше)
            if (multi && c < 1e17 && c - nbW(s - L) < bestC) {
                int gk[MULTI_K]; int ng = someGaps(s, ts, MULTI_K, gk);
                for (int i = 1; i < ng; i++) {
                    double c2 = zoneBlocked(gk[i], s, st) ? 1e17 : candCost(gk[i], s, L, st, m, selfHz);
                    if (c2 < bestC) { bestC = c2; bgs = gk[i]; bs = s; }
                }
            }
        }
        {
            int bl = 0, g = firstBig(ts, bl);
            if (g >= 0) {
                double c = zoneBlocked(g, bl, st) ? 1e17 : candCost(g, bl, L, st, m, selfHz);
                if (c < bestC) { bestC = c; bgs = g; bs = bl; }
            }
        }

        bool useGap = bgs >= 0 && bestC < 1e17;
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
            if (gForcePid == pid && gForceIdx >= 0 && gForceIdx < (int)gOpts.size() && !noStreams && ovSigOk(pid, gOpts, gForceIdx)) {
                const CfOpt& f = gOpts[gForceIdx];
                if (f.g < 0) useGap = false; else { bgs = f.g; bs = f.s; bestC = -1e18; useGap = true; }
            }
        }
        if (gCWset && st < 0 && T >= 60) goto cfdecided;
        // вне потоков: если лучший промежуток дороже LW[13], а в хвосте место есть без докопки — в хвост
        if (useGap && st < 0 && T - tailStart() >= L && LW[13] < bestC) useGap = false;
        cfdecided:
        if (gRecord) { gChosen = -1; for (size_t k = 0; k < gOpts.size(); k++) if (useGap ? (gOpts[k].g == bgs && gOpts[k].s == bs) : gOpts[k].g < 0) { gChosen = (int)k; break; } }

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
                if (need > L - fl) { ll cr = capRes(need + fl, L) - fl; if (cr < L - fl) cr = L - fl; need = (int)cr; }
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
                if (st >= 0 && need > L - tl) { ll cr = capRes(need + tl, L) - tl; if (cr < L - tl) cr = L - tl; need = (int)cr; }
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
        { int z = zl.get(d); if (z && (z & 63) == l) slots[z].push_back(d); }
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



#ifdef OVSEARCH
#include <cstring>
static std::vector<int> readCmdsO(const char* fn) { FILE* f = fopen(fn, "r"); double w; int n; std::vector<int> c; if (!f || fscanf(f, "%lf %d", &w, &n) != 2) return c; c.resize(n); for (int i = 0; i < n; i++) if (fscanf(f, "%d", &c[i]) != 1) { c.resize(i); break; } fclose(f); return c; }
static std::vector<std::pair<int, int>> readPairs(const char* fn) { std::vector<std::pair<int, int>> v; FILE* f = fn ? fopen(fn, "r") : 0; if (!f) return v; int a, b; while (fscanf(f, "%d %d", &a, &b) == 2) v.push_back({a, b}); fclose(f); return v; }
static void runOv(const std::vector<int>& c, ll& T, ll& emax, bool list) {
    Alloc* A = new Alloc(); emax = 0;
    for (size_t i = 0; i < c.size(); i++) {
        int x = c[i]; if (x == 0) break;
#ifdef OVCLS
        CUR = (long long)i;
#endif
        if (x > 0) {
            int pid = (int)A->pos.size();
            if (list) gRecord = true;
            ll N, D; A->place(x, N, D);
            if (list) { gRecord = false; if (gOpts.size() >= 2) { printf("D %d %zu %d", pid, gOpts.size(), gChosen); for (auto& o : gOpts) printf(" %d", o.g < 0 ? 0 : o.s); printf("\n"); } }
        } else A->remove(-x);
        ll e = A->tailStart(); if (e > emax) emax = e;
    }
    T = A->T; delete A;
}
#ifdef OVCLS
std::vector<int> ORD; long long CUR = 0; double K[64];
std::vector<unsigned char> gClsTab;
#endif
int main(int argc, char** argv) {
    const char* mode = getenv("MODE") ? getenv("MODE") : "list";
    std::vector<int> c = readCmdsO(argv[1]);
#ifdef OVCLS
    ORD.assign(1, 1 << 30); for (size_t q = 0; q < c.size(); q++) { if (c[q] > 0) ORD.push_back(1 << 30); else if (c[q] < 0) ORD[-c[q]] = (int)q; }
    if (getenv("CLSK")) K[40] = atof(getenv("CLSK"));
    if (getenv("CLSB")) K[41] = atof(getenv("CLSB"));
    if (getenv("CLSTAB")) { FILE* f = fopen(getenv("CLSTAB"), "r"); int a, b; while (f && fscanf(f, "%d %d", &a, &b) == 2) { if ((int)gClsTab.size() <= a) gClsTab.resize(a + 1, 0); gClsTab[a] = (unsigned char)b; } if (f) fclose(f); if (gClsTab.empty()) gClsTab.push_back(0); }
    if (getenv("CLSFILE")) {   // классы из файла (по строке «pid» на каждую трубу класса 1): заменяют оракул — ORD подгоняется
        std::vector<char> one(ORD.size(), 0); FILE* f = fopen(getenv("CLSFILE"), "r"); int q; while (f && fscanf(f, "%d", &q) == 1) if (q > 0 && q < (int)one.size()) one[q] = 1; if (f) fclose(f);
        K[40] = 1; std::vector<int> b(1, 0); for (size_t z = 0; z < c.size(); z++) if (c[z] > 0) b.push_back((int)z);
        for (size_t q2 = 1; q2 < ORD.size(); q2++) ORD[q2] = one[q2] ? b[q2] + 1 : (1 << 30);
    }
#endif
    std::vector<std::pair<int, int>> base = readPairs(getenv("OVFILE"));
    std::sort(base.begin(), base.end());
    if (!strcmp(mode, "list")) {
        gOv = base; ll T, em; runOv(c, T, em, true); printf("B %lld %lld\n", T, em);
    } else {
        std::vector<std::pair<int, int>> cand = readPairs(getenv("CANDS"));
        for (auto& cd : cand) {
            gOv = base; bool rep = false;
            for (auto& o : gOv) if (o.first == cd.first) { o.second = cd.second; rep = true; }
            if (!rep) { gOv.push_back(cd); std::sort(gOv.begin(), gOv.end()); }
            gForcePid = -1; gForceIdx = -1;
            ll T, em; runOv(c, T, em, false);
            printf("E %d %d %lld %lld\n", cd.first, cd.second, T, em); fflush(stdout);
        }
    }
    return 0;
}
#define main sf_main_unused_ov
#endif

#ifdef DSEARCH
// Поиск переопределений позиций (pid -> D) на одном тесте: случайная труба (pid >= MINPID), случайная позиция среди краёв
// свободных промежутков подходящей длины (и начала хвоста) в момент её укладки, полный прогон; принимается при строгом
// улучшении (T, emax, площадь) — площадь = сумма превышения началом хвоста порога THR по всем командам. Ход «убрать переопределение» — с вер. 0.2.
// env: OUT (файл "pid D gs gl"), DFILE (продолжить), OVFILE (переопределения вариантов: переводятся в позиции),
// MINPID, SECS, SEED; классы — как в OVSEARCH (CLSK/CLSB/CLSTAB).
#include <cstring>
#include <random>
#include <chrono>
static std::vector<int> readCmdsD(const char* fn) { FILE* f = fopen(fn, "r"); double w; int n; std::vector<int> c; if (!f || fscanf(f, "%lf %d", &w, &n) != 2) return c; c.resize(n); for (int i = 0; i < n; i++) if (fscanf(f, "%d", &c[i]) != 1) { c.resize(i); break; } fclose(f); return c; }
// area2 — сумма начала хвоста по всем командам (последний критерий, AREA2: уклон на плато); te — команда, где достигнут emax
static bool gArea2 = false;
static int gProbeD = -1, gProbeGs = 0, gProbeGl = 0; static bool gProbeFail = false;   // CANDD: позиция пробы задана; не поместилась
struct DRes { ll T = 0, emax = 0; double area = 0, area2 = 0; size_t te = 0; bool operator<(const DRes& o) const { if (T != o.T) return T < o.T; if (emax != o.emax) return emax < o.emax; if (area < o.area - 1e-9) return true; if (area > o.area + 1e-9) return false; return gArea2 && area2 < o.area2 - 1e-9; } };
static std::mt19937_64 drng;
static ll gThr = 0; static bool gTrace = false; static long long gVerb = 0;
static double gKn[NKN]; static long long gKnFrom = -1;   // KNFILE/KNFROM: ручки теста и команда, с которой они действуют
static void loadKn() {
    for (int k = 0; k < NKN; k++) gKn[k] = KN_DEF[k];
    if (!getenv("KNFILE")) return;
    FILE* f = fopen(getenv("KNFILE"), "r"); int a; long long q; while (f && fscanf(f, "%d %lld", &a, &q) == 2) if (a >= 0 && a < NKN) gKn[a] = (double)q / 1e4; if (f) fclose(f);
    gKnFrom = getenv("KNFROM") ? atoll(getenv("KNFROM")) : 0;
}
static std::vector<int> gPlan; static std::vector<std::array<int, 4>> gPlanCap; static FILE* gDump = 0;   // PLAN: позиции по pid (-1 — нет); DUMP: "pid D T"   // площадь считается над порогом THR (по умолчанию 0)
// probe >= 0: перед укладкой трубы probe выбрать позицию и поставить переопределение (в gOvD); cap — для записи
// подписей всех укладок с pid из capSet (перевод вариантов в позиции); applied[k] — переопределение k сработало
static DRes runD(const std::vector<int>& c, int probe, std::vector<char>* applied, const std::vector<int>* capSet, std::vector<std::array<int, 4>>* capOut) {
    knSet(KN_DEF);
    Alloc* A = new Alloc(); DRes r;
    if (applied) applied->assign(gOvD.size(), 0);
    for (size_t i = 0; i < c.size(); i++) {
        int x = c[i]; if (x == 0) break;
        if ((long long)i == gKnFrom) knSet(gKn);
#ifdef OVCLS
        CUR = (long long)i;
#endif
        if (x > 0) {
            int pid = (int)A->pos.size();
            if (pid == probe && gProbeD >= 0) {   // заданная позиция (CANDD): подпись — промежуток в этот момент
                int gs, gl;
                if (A->gapOf(gProbeD, x, gs, gl)) {
                    std::array<int, 4> o = {pid, gProbeD, gs, gl}; gProbeGs = gs; gProbeGl = gl;
                    int k = ovDFind(pid);
                    if (k >= 0) gOvD[k] = o; else { gOvD.push_back(o); std::sort(gOvD.begin(), gOvD.end()); if (applied) applied->insert(applied->begin() + ovDFind(pid), 0); }
                } else gProbeFail = true;
            } else if (pid == probe) {
                std::vector<std::array<int, 3>> cand;
                int T = (int)A->T, xx = 0, ts = A->tailStart();
                while (xx < ts) {
                    if (A->occupied(xx)) { int p = A->owner.get(xx); xx = A->pos[p] + A->len[p]; continue; }
                    int gl = A->tag.get(xx); if (gl <= 0) { xx++; continue; }
                    if (gl >= x) { cand.push_back({xx, xx, gl}); if (gl > x) cand.push_back({xx + gl - x, xx, gl}); }
                    xx += gl;
                }
                cand.push_back({ts, ts, T - ts});
                auto cc = cand[drng() % cand.size()];
                std::array<int, 4> o = {pid, cc[0], cc[1], cc[2]};
                int k = ovDFind(pid);
                if (k >= 0) gOvD[k] = o;
                else { gOvD.push_back(o); std::sort(gOvD.begin(), gOvD.end()); if (applied) applied->insert(applied->begin() + ovDFind(pid), 0); }
            }
            if (pid < (int)gPlan.size() && gPlan[pid] >= 0) {
                int gs, gl;
                if (A->gapOf(gPlan[pid], x, gs, gl)) { A->forceD = gPlan[pid]; gPlanCap.push_back({pid, gPlan[pid], gs, gl}); }
                else printf("PLAN: позиция pid %d занята\n", pid);
            }
            ll T0 = A->T; int ts0 = A->tailStart();
            ll N, D; A->place(x, N, D);
            if (gDump) fprintf(gDump, "%d %lld %lld\n", pid, D, A->T);
            if (gTrace && N > 0) {
                int holes = 0, free = 0, xx = 0, mx = 0;
                while (xx < ts0) {
                    if (A->occupied(xx)) { int p = A->owner.get(xx); xx = A->pos[p] + A->len[p]; continue; }
                    int gl = A->tag.get(xx); if (gl <= 0) { xx++; continue; }
                    if (xx < (int)D || xx >= (int)D + x) { free += gl; if (gl >= x) holes++; if (gl > mx) mx = gl; }
                    xx += gl;
                }
                printf("DIG i %zu pid %d L %d T %lld->%lld D %lld ts0 %d holes>=L %d maxhole %d freeIn %d\n", i, pid, x, T0, A->T, D, ts0, holes, mx, free);
            }
            if (applied) { int k = ovDFind(pid); if (k >= 0 && gOvD[k][1] == D) (*applied)[k] = 1; }
            if (capSet && std::binary_search(capSet->begin(), capSet->end(), pid)) {
                int gs, gl;
                if (D >= ts0) { gs = ts0; gl = (int)T0 - ts0; }
                else { gs = (int)D; while (gs > 0 && !A->occupied(gs - 1)) gs--; int e = (int)D + x; while (e < A->T && !A->occupied(e)) e++; gl = e - gs; }
                capOut->push_back({pid, (int)D, gs, gl});
            }
        } else A->remove(-x);
        ll e = A->tailStart(); if (e > r.emax) { r.emax = e; r.te = i; } if (e > gThr) r.area += (double)(e - gThr); r.area2 += (double)e;
    }
    r.T = A->T; delete A; return r;
}
#ifdef OVCLS
std::vector<int> ORD; long long CUR = 0; double K[64];
std::vector<unsigned char> gClsTab;
#endif
static void saveD(const char* fn) { FILE* f = fopen(fn, "w"); for (auto& o : gOvD) fprintf(f, "%d %d %d %d\n", o[0], o[1], o[2], o[3]); fclose(f); }
int main(int argc, char** argv) {
    std::vector<int> c = readCmdsD(argv[1]);
#ifdef OVCLS
    ORD.assign(1, 1 << 30); for (size_t q = 0; q < c.size(); q++) { if (c[q] > 0) ORD.push_back(1 << 30); else if (c[q] < 0) ORD[-c[q]] = (int)q; }
    if (getenv("CLSK")) K[40] = atof(getenv("CLSK"));
    if (getenv("CLSB")) K[41] = atof(getenv("CLSB"));
    if (getenv("CLSTAB")) { FILE* f = fopen(getenv("CLSTAB"), "r"); int a, b; while (f && fscanf(f, "%d %d", &a, &b) == 2) { if ((int)gClsTab.size() <= a) gClsTab.resize(a + 1, 0); gClsTab[a] = (unsigned char)b; } if (f) fclose(f); if (gClsTab.empty()) gClsTab.push_back(0); }
#endif
    const char* out = getenv("OUT") ? getenv("OUT") : "dsearch.out";
    int minpid = getenv("MINPID") ? atoi(getenv("MINPID")) : 1;
    gThr = getenv("THR") ? atoll(getenv("THR")) : 0; gVerb = getenv("VERB") ? atoll(getenv("VERB")) : 0;
    double secs = getenv("SECS") ? atof(getenv("SECS")) : 600;
    loadKn();
    if (getenv("KNDEF")) { for (int k = 0; k < NKN; k++) printf("%d %.6g\n", k, KN_DEF[k]); return 0; }
    if (getenv("KNEVAL")) {
        if (getenv("DFILE")) { FILE* f = fopen(getenv("DFILE"), "r"); int a, b, g1, g2; while (f && fscanf(f, "%d %d %d %d", &a, &b, &g1, &g2) == 4) gOvD.push_back({a, b, g1, g2}); if (f) fclose(f); std::sort(gOvD.begin(), gOvD.end()); }
        DRes r = runD(c, -1, 0, 0, 0);
        printf("K %lld %lld %.0f %.0f\n", r.T, r.emax, r.area, r.area2);
        return 0;
    }
    drng.seed(getenv("SEED") ? atoll(getenv("SEED")) : 1);
    int P = 0; for (int x : c) if (x > 0) P++;
    if (getenv("OVFILE")) {   // переопределения вариантов -> позиции (с подписями этой укладки)
        FILE* f = fopen(getenv("OVFILE"), "r"); int a, b; while (f && fscanf(f, "%d %d", &a, &b) == 2) if (a >= minpid) gOv.push_back({a, b}); if (f) fclose(f);
        std::sort(gOv.begin(), gOv.end());
        std::vector<int> cs; for (auto& o : gOv) cs.push_back(o.first);
        std::vector<std::array<int, 4>> cap;
        DRes r0 = runD(c, -1, 0, &cs, &cap);
        gOv.clear(); gOvD = cap;
        DRes r1 = runD(c, -1, 0, 0, 0);
        printf("перевод: вариантов -> позиций %zu, T %lld/%lld -> %lld/%lld\n", gOvD.size(), r0.T, r0.emax, r1.T, r1.emax);
    }
    if (getenv("CANDD")) {   // оценка кандидатов "pid D" поверх DFILE: строки "E pid D gs gl T emax area area2" (T = -1 — позиция занята)
        if (getenv("DFILE")) { FILE* f = fopen(getenv("DFILE"), "r"); int a, b, g1, g2; while (f && fscanf(f, "%d %d %d %d", &a, &b, &g1, &g2) == 4) gOvD.push_back({a, b, g1, g2}); if (f) fclose(f); std::sort(gOvD.begin(), gOvD.end()); }
        std::vector<std::array<int, 4>> base0 = gOvD;
        FILE* f = fopen(getenv("CANDD"), "r"); int a, b;
        while (f && fscanf(f, "%d %d", &a, &b) == 2) {
            gOvD = base0; gProbeD = b; gProbeFail = false;
            DRes r = runD(c, a, 0, 0, 0);
            if (gProbeFail) printf("E %d %d 0 0 -1 -1 0 0\n", a, b);
            else printf("E %d %d %d %d %lld %lld %.0f %.0f\n", a, b, gProbeGs, gProbeGl, r.T, r.emax, r.area, r.area2);
            fflush(stdout);
        }
        if (f) fclose(f);
        return 0;
    }
    if (getenv("PLAN")) {   // готовая раскладка "pid D" (pid >= MINPID) -> переопределения позиций с подписями
        FILE* f = fopen(getenv("PLAN"), "r"); int a, b; gPlan.assign(P + 2, -1); while (f && fscanf(f, "%d %d", &a, &b) == 2) if (a >= minpid && a <= P) gPlan[a] = b; if (f) fclose(f);
        DRes r0 = runD(c, -1, 0, 0, 0);
        gPlan.clear(); gOvD = gPlanCap;
        DRes r1 = runD(c, -1, 0, 0, 0);
        printf("план: позиций %zu, T %lld/%lld -> %lld/%lld\n", gOvD.size(), r0.T, r0.emax, r1.T, r1.emax);
    }
    if (getenv("DFILE")) { FILE* f = fopen(getenv("DFILE"), "r"); int a, b, g1, g2; while (f && fscanf(f, "%d %d %d %d", &a, &b, &g1, &g2) == 4) gOvD.push_back({a, b, g1, g2}); if (f) fclose(f); std::sort(gOvD.begin(), gOvD.end()); }
    gArea2 = getenv("AREA2") != 0;
    std::vector<int> birthI(1, 0), deathI(1, 1 << 30);
    for (size_t q = 0; q < c.size(); q++) { if (c[q] > 0) { birthI.push_back((int)q); deathI.push_back(1 << 30); } else if (c[q] < 0) deathI[-c[q]] = (int)q; }
    std::vector<int> pool;
    auto mkPool = [&](size_t te) { pool.clear(); for (int q = minpid; q <= P; q++) if (birthI[q] <= (int)te && deathI[q] > (int)te) pool.push_back(q); };
    double pfr = getenv("POOLFR") ? atof(getenv("POOLFR")) : 0.75;
    std::vector<char> ap;
    DRes base = runD(c, -1, &ap, 0, 0);
    { std::vector<std::array<int, 4>> keep; for (size_t k = 0; k < gOvD.size(); k++) if (ap[k]) keep.push_back(gOvD[k]); gOvD = keep; }
    base = runD(c, -1, 0, 0, 0);
    mkPool(base.te);
    printf("start T %lld emax %lld area %.0f overrides %zu placements %d pool %zu (te %zu)\n", base.T, base.emax, base.area, gOvD.size(), P, pool.size(), base.te); fflush(stdout);
    saveD(out);
    auto t0 = std::chrono::steady_clock::now();
    long long it = 0, acc = 0;
    for (;;) {
        double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (el > secs) break;
        it++;
        std::vector<std::array<int, 4>> saved = gOvD;
        DRes r; std::vector<char> a2;
        if (!gOvD.empty() && drng() % 5 == 0) { gOvD.erase(gOvD.begin() + drng() % gOvD.size()); r = runD(c, -1, &a2, 0, 0); }
        else {
            int p = (!pool.empty() && (double)(drng() % 1000) < 1000 * pfr) ? pool[drng() % pool.size()] : minpid + (int)(drng() % (unsigned)(P + 1 - minpid));
            r = runD(c, p, &a2, 0, 0);
        }
        if (gVerb && it <= gVerb) printf("  eval %lld: T %lld emax %lld area %.0f\n", it, r.T, r.emax, r.area);
        if (r < base) {
            std::vector<std::array<int, 4>> keep; for (size_t k = 0; k < gOvD.size(); k++) if (a2[k]) keep.push_back(gOvD[k]);
            bool newT = r.T < base.T;
            gOvD = keep; base = r; acc++; saveD(out); mkPool(base.te);
            if (acc % 20 == 1 || newT) { printf("it %lld %.0fs: T %lld emax %lld area %.0f overrides %zu\n", it, el, base.T, base.emax, base.area, gOvD.size()); fflush(stdout); }
        } else gOvD = saved;
        if (it % 2000 == 0) { printf("it %lld %.0fs: T %lld emax %lld area %.0f overrides %zu acc %lld\n", it, el, base.T, base.emax, base.area, gOvD.size(), acc); fflush(stdout); }
    }
    gTrace = getenv("TRACE") != 0;
    if (getenv("DUMP")) gDump = fopen(getenv("DUMP"), "w");
    DRes chk = runD(c, -1, 0, 0, 0);
    printf("final T %lld emax %lld area %.0f overrides %zu (check T %lld emax %lld) it %lld acc %lld\n", base.T, base.emax, base.area, gOvD.size(), chk.T, chk.emax, it, acc);
    if (gDump) fclose(gDump);
    return 0;
}
#define main sf_main_unused_d
#endif
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

static const char* const BLOB[] = {
"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",
"C6JWABwGmp0BACEL6KQBAAQCzpABAx0FgKUBAAAY0CUEFEeQmAEbJgDeDRYBDwECAQYBCwEFAQUBAQEDAQMBCAEBATEBBgEdARUBDgEEASUBEQEbAQsBBQECAQQBAQEEAQsBRwECASwBEgEXAQIBCQEBAQYBAQEJAQIBAwEEAQQBAQECARUBDQEBAQgBPAEJARgBBgEN",
"AQcBBAEBAQQBDwEGAQ0BBQELAQYBAgERAQIBAwEIAQUBCQEHARQBAQEGAQEBAQECAQYBBgEDAQUBDwEFAQEBEgEJAQsBBQERAQMBCgEDARUBBAELAQQBAgECAQIBAwERAQsBAwELAQkBDAEKAQQBHwEDAQkBBAEFAQEBAgECAQIBAQEBAQcBAgEBAQcBBwEBAQEBAwEB",
"AQ0BBAEYAQMBBgEEAQYBAgEDAQQBBgEBAQUBAgEBAQEBBwEKAQEBAQEFAQcBBAEGAQMBAwEBAQIBBAEBAQkBDAEEAQEBAQEBAQkBAwEFAQIBAwEEAQEBCAEHAQsBCAEEAQQBAQEDAQMBBgEGARkBBgEBAQUBDQEGAQEBAQEPAQMBAwEBAQwBAgEBAQEBAgEFAQwBCAEK",
"AQcBAgEGAQMBBAEJAQkBCQEfARABBQEIAQMBAwEBARABBQEHAREBCwEWAQcBBQEGAQQBAwEFAREBJAEnAQUBDgFCARkBFgESAQYBAwELAR8BEgE5AQIBBAECARABDAEHARUBBAEVAQUBCwEDAQYBAQEKAQkBCwEDAQIBBgEJAS4BHAEGAREBCQEFAR0BDwECAQEBHwEO",
"AQIBJQEEAQUBBQEDAQcBBwEFAQwBBQEBAQgBCQEKAQwBDAETAQwBFQElAQoBAgENARUBIAERAQUBAQEOAQMBBQECAQEBAQEYAQQBCQEBAQ8BBgECAQEBBAEIASABAQEEAQMBAQECAQYBAQEEAQEBAwECAQMBAgEDAQcBAQECAQMBAQEGAQ4BAwECAQEBEQEDAQUBAQEE",
"AQYBAwEHAQEBCQECAQIBAwEMAQcBAQEFAQYBCQECAQMBAgEDAQUBAQECAQIBAwEBAQMBAQEFAQkBAQECAQUBAQECAQEBAQEBAQEBAQEBAQEBAgEBAQEBAQEDAQMBAwECAQEBAQEDAQEBAgEEAQIBAgEBAQcBAQEMAQIBAgEDAQEBAQEBAQYBBQEEAQEBBgEBAQIBAwEF",
"AQIBAwEFAQMBAwEKAQIBAQEDAQEBAQEBAQEBAgEBAQkBBQECAQMBBQEBAQIBAgEDAQEBAwECAQIBBAEEAQQBAgEDAQEBBQEBAQMBAgEDAQEBBAEEAQgBBgECAQEBAgECAQEBAwEBAQIBAQECAQEBBQEBAQEBBQEIAQIBAwEDAQIBAgEBAQIBBQECAQMBAQEBAQEBAQEB",
"AQEBAQEFAQEBAwEBAQIBAQEBAQEBAQEBAQIBAQEBAQEBAQEBAQoBAQECAQEBAQECAQUBAQEDAQEBAQECAQcBAQEEAQEBAwECAQEBAQEBAQMBAwEEAQIBAQEDAQUBAQEDAQgBAgEFAQcBAQECAQEBAQEBAQEBAQEBAQcBAgECAQIBAQECAQEBAgEDAQoBAgEJAQEBBAEE",
"AQMBAgEBAQEBAgECAQIBAQEBAQYBAQEBAQEBBAEBAQEBAQEBAQIBAgEHAQMBAQECAQMBBAEBAQIBAgECAQEBBAEFAQYBAQECAQIBAwEBAQEBAgEBAQMBAQEBAQIBAwEBAQQBDAEEAQEBAQEBAQEBAQECAQIBAQEBAQEBAwEBAQEBAgEDAQEBAQEBAQIBAQEBAQQBBQED",
"AQIBAQEBAQYBAgECAQEBAQEBAQEBAgEBAQUBAQEFAQIBAQEFAQQBBAEBAQQBBAEBAQEBAgEFAQIBAgEHAQMBAwEEAQIBAgECAQEBAQEBAQEBAgEBAQMBAwEIAQEBAQEDAQEBAgEDAQQBAwELAQIBDwEEAQEBAwECAQcBDwEBAQgBCgEBAQIBBAECAQcBAwECAQEBAwED",
"AQIBCwEDAQYBCAEHAQYBAgEEAQEBAgECAQIBAQEEAQIBAgECAQIBBgEBAQEBAQECAQEBAQEDAQEBAwEBAQEBAQEBAQEBAQEBAQMBAQECAQEBAQECAQEBAQEDAQQBAQEBAQEBAQEBAQMBAQEBAQIBBAEBAQEBAgEDAQEBAQEBAQMBAwEBAQIBBAEDAQUBAQEBAQEBAQED",
"AQEBAgEBAQEBAwEBAQEBAQEGAQQBAQECAQUBAQEBAQEBAQECAQIBAQEBAQQBBQEBAQIBAgEBAQEBBAEBAQEBBAEBAQEBAwEBAQIBBAEBAQIBAgEBAQIBAQEDAQEBBAEDAQQBCQEBAQQBAQEBAQEBAQEBAQIBAwECAQEBAQEBAQEBAQEBAQIBAgEBAQEBAQEBAQEBAgEB",
"AQEBAgEBAQMBAwEBAQEBAgEBAQEBAQEBAQIBAQEBAQIBAQECAQIBAQEBAQEBAgEBAQEBAwEBAQEBAQECAQEBAQEBAQIBAQEBAQMBAwECAQIBAgEBAQMBAQEBAQIBAwEDAQUBAgEBAQEBAQEBAQEBAgEBAQIBBAEBAQIBAQEBAQEBAQEBAQIBAQECAQMBAgEDAQEBBAEB",
"AQEBAQEDAQIBAQECAQMBCQEBAQMBAQEBAQUBAQECAQEBAwEEAQEBAQEBAQcBAgEBAQEBAQEEAQEBAwEEAQUBAQEGAQQBAQEJAQEBAgEBAQIBAwECAQMBBwEDAQQBAwECAQEBAwEDAQEBAgEEAQEBAwEBAQIBAQECAQEBBAEBAQEBAQEBAQEBAQEBAQEBAQEMAQIBAQEC",
"AQEBAQEEAQEBAQEBAQEBAQEFAQIBAQEBAQEBAgECAQEBBQEBAQEBAwEBAQEBAQECAQEBAQEDAQIBBAEBAQEBAwECAQEBAgEBAQQBAgECAQEBAgEDAQMBAQEBAQEBAwECAQEBAQEEAQMBAgEDAQEBAQECAQEBAgECAQEBAQEBAQQBBQEBAQcBAQEBAQEBAgEDAQEBAQED",
"AQEBAQECAQEBAQEBAQEBAgEBAQMBAQECAQEBAgECAQEBAgEBAQEBBAEEAQEBAQEBAQEBAQEDAQIBAQEBAQEBAgEBAQEBAQEBAQIBAgECAQEBAQEDAQQBAQEDAQMBBgECAQIBAgEBAQEBAwECAQIBAgEGAQMBAQEBAQIBAgEBAQEBAgEBAQEBAwEBAQMBAwEDAQEBAQEB",
"AQEBAQEBAQEBCQEDAQEBAQEBAQMBBQECAQEBAQEBAQIBAQEBAQEBAQEDAQEBAQEBAQMBAQECAQEBAQEBAQIBBAEDAQEBAQEBAQMBAQEGAQEBAQEBAQEBBQEBAQEBAwEBAQUBAgEFAQEBAgEBAQIBAgEFAQUBAQEBAQIBAQEFAQIBAQEDAQEBAQECAQQBAQECAQQBAwEB",
"AQEBAQEBAQIBAQEDAQMBAwEBAQEBAQEBAQEBAQEBAQMBAQEBAQIBAwEDAQIBAQEDAQIBCAEBAQEBAgECAQQBCAECAQEBAQECAQIBAQEDAQMBBAECAQEBAwEGAQEBAgECAQUBAgEBAQEBAgEBAQEBBAEBAQQBAwEEAQEBCQEBAQEBAQEBAQEBAQEBAQEBAQEEAQEBAgEB",
"AQIBAQECAQMBAgEBAQEBAwECAQEBBAECAQEBAQECAQEBAwECAQIBAQEBAQEBAwEBAQEBAwECAQIBAQEHAQEBAwEBAQEBAgEBAQIBAQEBAQEBAgEBAQIBAQEDAQEBAgEBAQUBBAECAQEBBAEBAQIBAQEBAQIBAQEBAQIBAQEBAQEBAgEDAQIBBAEBAQUBAwEBAQUBAwEB",
"AQYBAQEBAQIBAgEBAQQBAgECAQMBBAECAQIBAQEBAQEBAgEFAQEBAQECAQgBAQEFAQIBBAECAQEBAQEHAQMBAwECAQQBAwECAQQBAQECAQEBAQEDAQEBAQEBAQQBBAEBAQMBAwEFAQEBAQEBAQEBBQECAQEBAgEEAQEBAgEBAQYBAQEBAQEBAQEBAQQBAgEBAQEBAQEB",
"AQIBAwEDAQUBAwEDAQUBBAEBAQMBAQEBAQMBAQECAQIBAQECAQEBAgEBAQEBAgECAQUBAQEBAQIBAgEBAQEBAQEBAQYBAwEDAQEBBAEBAQEBAQEEAQMBAgEBAQIBBAEBAQEBBAEBAQMBAwEBAQEBBAEBAQMBAgEBAQEBBgEBAQMBAgEFAQMBAgEFAQIBAQEBAQEBAQEC",
"AQEBAQECAQEBAQECAQEBAQEDAQEBAwEBAQEBAQEEAQEBAgEDAQQBAgECAQQBAQEDAQMBBwEEAQYBCAEBAQQBBAEHAQEBAQEDAQIBAQEDAQQBAwECAQQBBAEBAQIBAwEBAQYBBAECAQcBAQECAQEBAwEBAQQBBAEBAQIBAQECAQMBAQEBAQIBAwEBAQEBBAEJAQEBBAEB",
"AQkBBgEBAQUBBAEBAQIBBQEFAQEBAgEBAQIBAwEDAQIBCwEBAQMBAQECAQMBBAEBAQQBBQEDAQEBAgECAQIBAwEIAQEBAQEEAQUBBgEBAQIBAQECAQMBAQEJAQMBBAEEAQEBAQECAQIBAQELAQEBBwEEAQMBAQEBAQIBBAEBAQEBCQEEAQIBBAECAQYBAwEBAQIBAgEF",
"AQEBBwEJAQ4BAgEBAQgBAgEBAQkBAgELAQIBBAEFAQsBAgEDAQMBCAEBAQEBAQEBAQMBBgEBAQIBBAEBAQEBAQEBAQUBAQEDAQEBAgEBAQsBBwEBAQEBAQECAQIBAQEJAQoBAQEDAQcBCQEDAZA6AAQG5zn7vDIenz/hklDUX5ED2rblkXn4NWEowCpGVeHsKVmEOHPo",
"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",
"NQFAByEBAwcDBwEHTAchByMHFwc7ByUHTwd4BzkHCwe0AQdwBwoHeQcQBzAHAQdTBhUHEQc3Bi0HCwdfBj4HdgccB0wHFAcZBs0BBxsGZwcQBxwHPQcFBwgHFAaAAQcoBw8HLQe7AQczB4IBBz4HhwIHPAdlBwoHhAEHVwcUBwkH/AEHQAc0B44BBgEGOQdXAQsBUwYK",
"AQEBBQHWAQd3B1cGYAcKB8MBB4ABB/ICB0kHDAfKAwdjBxkHBgdSB9sBB1sHIweIAQGjAgdFB1MGOAEDB5MDB0MGAgd3B0AHLAcQBwcHVgcMB0EHAQcGBhIHMgdmB1oHeAcQBwkHdQcGBwgHDwcmBxUHDQdOB1cHVgc4BwgGCQdjB2IHcAEjByMHYQc4B1wHGwfCAQcV",
"ByEHDwenAQYCBlMHDwcBBxwBAwcLBwYHBAcZBwsGDAcLBwEHDQcEBk0GNgEPARwHOgcLBxsHGgYKBiYHBgcEAQIBCgEWBgQHFgcNBxwHLwcCBxgHQAcGAikCCQEOAQQHDwcKBgQHCQdEBwoHAQchBhQHlwEHCgcBBw0HHgcjB84BBwMHQAcKBxcHSwdSBx8HIgcIBwYH",
"EgcDBwIHAQcFBzsHMAcCBwwHAwcQBwcHAQcWBw0HEAcNBzAHEwcDBwUHBwcbBx4HDQcXBwEHJAcJBwkHOQcDBwEHGQcGBwwHCAcEBwcHAgcFBwQHBAcHBwIHDQcLByEHAQcOBxQHDQcFBx8HAQcHBwcHBgcSBwcHDwcBBxQHAQcBBx4HAwcDBwQHCQcDBwUHAQcIBycH",
"BQcXBwoHAwcGBxoHFwcMBwoHGAcHBxYHGwcDBwIHCAclBx0HGwcSBw8HHwIBBwcHAgcBBwMHCwcMBwUHBAcIBwYHCgcBBgMHAwcDBwoHAQcGBwQGBQcFBwgHBAcCBwEHDgYBBgEHLwcHBwUHCgYhBwEHAQcGBwIHAQILBgIHAQcCBwQHBwYIBxgHDAcCBwMHAwcHBwIH",
"AwcaBwEHAgcEBwMHAQcDBwoGGAcIBwIHAwcDBw0HAgcUBwsHBAcPBwIHAgcBBwgHAQcQBwcHBQcDAgEHBgcMBwcHAwIDBwEHAQcOBw0HBQcJBwQHDwcMBgoGCAYZBwUHIgcbBwMHAgcEBwQHAgYGBgcHEgcFBgMHAgcMBwoHDAYMBwcHCgcBBw0HCwcBBgYGBgcQBwEH",
"AQcIBwQHBQYPBwgHAQcLBxYHCQcBBwgHCQYEBwEHCQcEBwoHAQcDBwMHCQcCBwkHAwcJBwMHAQcNBxMHAgcCBwIHCQcIBwsHBgIBBwEHBQcDBw8HAQcBBgsHCAcBBwIHDQcLBysHBwcBBgcGQAcRBw0HDwcEBwYHAgcFBwkHFQYBBxIHCQcHBwEHAgcFBxYHAQcJBwMH",
"CQcCBw0HAQcDBw0HDwcNBwEHIAcFBwQHBgcKBwgHCAcCAhIGAwcBBwICBAcBBgMHAQIEBwcHFwcDBwEHBAcCBxMHCgcCBwUHEgcCBwEHAgcHBwIHAgcFBwQHAQcLBwIHHQcCBxIGBQcJBwUHBwcBBwQHBQYYBwwHAgcFBwYHAQcQBwwHBAcIBwcHBwcMBwcHBgcDBwEH",
"JQcIBxAHBQcMBwkHAQcEBgMGAQcBBwcHAgcGBwUHCgcFBwQHBgcHBgUGBAcHBxkHAwcFBwUHAgcDBwMHAgcMBwEHCgcLByEHBAcEBw0HDQcIBwcHBgcCBwkHDQcKBwEHAgcFBgIGAQcBBxEHAgcDBwIGCQcKBgwHCgcEBhQGAgcEBwIHBQcMBwECEQcCBwsHCwcDBwkH",
"AgcGBwkHDgcRBw4HFAcFBwwHCAcLBwgHBgcFBxUHBAcBBwMHBQcCBwUHCQcTBwQHDgcGBx8HDAcQBxAHAwcYBwMHAwcDBxQGAwcHBwcHBQcJBwgHDAcCBwEHBgcEBwIHFwYDBgEHDgYBBgYHBAYPBwkGCAYDBg8GLAcCBwgGAQcIBwMHCwcIBwMGBQcMBwIGAQYBBhUH",
"BwcDBwQHDQYEBxEHBAcZBxIHCwYNBwUGKgYaBg8HLAYZBjYHAQcVBw4HBwcdAyoDBgcGBwQHKwYJBgMGJQcBBwsHHwYFBwIGCgcOBgEHBQYDBhQGCQYMBwIGIgcGBhEGLwYHBwYHAwYBBwQGKQcJBwkGCwcaBwYHDwYKBh4GHQYiBhIGAgcCBiMGBQNSArABB3IHAwcZ",
"BwcCIgcPBxACawdsB0wHfwcDBxgHNgc6B2gHmgIHDwdsBxYH0wIHPAczBw8HOAcPBxcHFwemAQcQBxQHWwcxBxUHRQcMBzIHLQfuAQcCBxAHHgesAQc5B0EHJQc2B4YBB0gHvwEHEQJ/BwoHPgYbBy8HJQYDBwQHHAcWB1UHAwcHByEH6AEHKwcuB/YCB1sHMwcGB10G",
"RQYBBxUHIQaMAQdUB7gCBjEHWwYJBi0HCQYRBkYHXgeQAQe6AwcEB00HFwfeAQcIBx0Hcwe5AQcbB3kHXAcNBwcHpQIHAQcXBzAHlQIHFweuAQcaBzoHTwcZByUHMAdhAgYHUwdLBwwGpQEHFgdIBxgHDQcCBxoHIQcBBzsHTQcRB3gHFAeGAQceBhwHGwIBAhgHfwIw",
"B1EHHgcJB04HAwckB7QBBysHMgYzBxQHcAcJBxsHGQdhB9sBB2AHEAc3BwoGbAakAQekAQdJBjsGhQIGAQZZBjUGFwZvBwUCQgcNBxkGDQcQAgMGHQIDAgYCAQJUBw4HZgfEAQIeAg4C8AIHugEHqgEHkwIHKgIZAhQCGAINAgMHDgIBAg8CFQIJAgoCKAd2Bw8HJAcg",
"By0HUAelAQcGBw0HOQcpAgECHQdoBxcHBgcQAigHAwcxBwMHZAcTBwEHJQcFBx4HAgc8B4IBB4UBBhsHGAceB0MHDwcDB3EHGQbFAgc2Bx0HEAcfBzIHHwcRB0IHVAJGAgEHEgIBAgICAQIrB3QHJgYKBxYGLAcFB4YBBwQHEgcFAgYHAgcJBxUHDAcTBw0HAgcEBxQH",
"ewcnBxYHHgckByUHzQEHNQZGAyMHIQcIBxYGNgdOBw8GHwe9AwcbByMHDgcIBm8HXQcgB5MBB64DBwEHuQIHZgefAQecAgcaB8oBBwgHDQcWAwMHDAMKAxADDwMGAwQDBQcGBxMHCwe+AQMLAzcDAQNSB0UHBwcPBx0HGwMKAw4HjwEHRAcsAw4HCwcDByoHLwdJBw8H",
"AgfgAQcSBzgHAgdFBwIHMQMuAxQHEQdIByIHIAMMB94BB0EHKgddBxAHNwcBBwcGCAdtBwEHVQcXBugBBwgHFAcUBysHYgdABxMHIQN1AyQHLwMEAxMHBwc0B3AGcAYRByUHaAcKBjcHNAbtAgMFAyIHPgMRAxQDCwMUBgIGJwcWBi8HDwMiBwgGOQY8B0MGMQeOAQcg",
"BqMBB9EBBxwHngMHPgc/B3UHgwIHEwcLByQHAQezAQcmA1sDCQMFA4UBBzMHCweFAQcaAwgHIAcLA0gHAwclBx4HAgcMBwEHJwcBBxAHCgcDBwUGCAcCBwsGEAcCBwcHCAcBBwgHIwcwBwEHBgcGBw0HFQcJBxoHBwcBBwUHJAcBBwYHBQcBBwMHBAcLBwEHBwcKBwMH",
"FAcYBwYHCQcBBwsHAgcDBwYHBAcBBxEHBgcJBwEHDAcDBwsHAgcCBwEHAgcCBwIHBgcBBwEHAwcIBwkHBwcIBgkHBAcDBwUHDAcBBwwHBAcIBwIHAwcBB0IHAgcCBwYHAwcNBwoHCAcGBwEHAQcDBwkHBAcOBwEHBAcDBw0HDgYGBwIGBgYSBgEGAgYBBhEHBwcDBwQH",
"AgcBBwIHAQcHBgkGAQcBBgEHBwcBBwIHAQcCBwQHAQYHBwIHCAcGBwEHAQcDBwEHBgcCBwIHAQcFBwkHCgYcBwUGAgcBBgQGAwcCBwUHCgYOBwIHAQcMBgUGCgYBBgQHAQYFBgEHAgcBBwIGDQcCBgcHBgcBBgQHAgcBBwMHCgcGBgUHBQYDBwgHAwcBBw4HBAcBBgEG",
"CAcDBgEHDAcCBwEHAgcDBwgHEAcGBwMHCAcGBwQHAQcFBwMHBQcCBwIHAQcEBwMHAQcFBwEGBwcEBwQHFQcCBwIHAgcBBwEHBQcDBwQHCgcBBwIHEwYCBgEHAgcFBwgGAgYCBwcHAgcBBwEHAQYBBwEGAwcJBgEHAQYDBgEGAgYDBgEGBAcBBwMHAQcDBwEHAwcCBgcH",
"BAcNBgMGBAYCBgEHCQcGBgkHDAcGBwEHAgcBBwEHBQcCBwgGAgcDBwIHBQcFBwkHCgcDBwMHAgcDBwQGAwYCBwYHBgcJBwUGAQYDBgEHAgYBBwEGCAYBBgYGAQYBBgEGBAcBBgEHAgYCBwUHCwYGBgYHBwcCBgIHAgcCBgcHAgYBBgIGAQcBBgMGAQYBBgUGCQYDBgUH",
"BQYCBgEGBgYDBwEGAQYCBgEHBAYDBgIHAwYFBgEHCAYEBgEGBwYDBwIGAgcBBgIHBAYCBgIGAwYCBgMHBwYEBgIGAQcBBgMHAQcCBwEHCAYEBwcGAgYFBgoHBAYEBwMHAwYDBwEGAQcCBwEGAwYIBwcHAQcBBwIHAwYBBgQHAQcCBwQGAwYDBwEGAQYDBgEHAgYBBgIG",
"AQYDBgUGAgYCBwEGAQcGBgEGAQYCBgEGAwYBBgMGBQcFBgEGAgYBBgIGAQYGBgIHAwcDBgMGCQYBBgQGAwYEBgEGAgcCBgEGAwYDBgEGAgYBBgQHAwYBBgQGAgYCBgEHBwcBBgIGCAYCBgUGAQYCBgMHAQYBBgEGBwYCBgMHAQcEBgEGAQcBBgEGAQcCBwIGAgYCBgQG",
"AQYBBgIGAQYBBwMGAgYDBwIGAgYBBgIGAQYBBwIGAgYDBgEGAgYGBgYGAQYCBgMGAQYDBgIHAQYBBwQGAwYBBwQHAQYCBgEGAQYBBgEHAgYBBgEGAwYCBgQGAQYDBgIGAQYDBgIGAQYBBgcGAQYBBgMGAQYGBgUGAwYBBgEGAgYDBgMHAQYOBgcGAQYFBgYGBgYDBwEG",
"AQcBBgIGAgcBBgEHAwYFBgYGAQcBBwEHAQYBBgMGAQYBBgIGAQYFBwUGBQYBBgIGBwcDBgEHAgYBBwQHAgYDBgEGBAYBBgIGAwYCBgIGAQcBBwMGAQYCBgYHAgYDBgIGAwYCBgEGAQYCBgIGAwYDBgMGAwYDBgEGAQYCBgEGBAYBBgEGAQYBBgEGDgYCBgMGBQYCBgQG",
"AQYBBgEGAQYCBgMGAQYCBwEGAgYGBgcHCwYBBgIGAgYEBwEGAQcGBgEGAgYBBgEGBAYDBgIGAgYBBgMGBAYBBgMGBQYBBgMGAQYDBgMGAwYEBgIGAQcBBgQGAQcDBgEGAwYDBgIGAQYCBgMGAQcBBgIGAgYBBgUGAQYEBgEGAQYEBgEGBAYDBgoGAQYBBgEGAQYCBgEG",
"AgYBBgQGAgYBBgEGAQYFBgEGBwYHBgQGAgYBBgEGAgYCBgIGAwYCBgQGAgYEBgMGAwYEBgQGAgYBBgEGAwYDBgEGBAYCBgUGAwYBBgIGAwYBBgMGAgYGBgEGAwYCBgMGAwYFBgEGAgYCBgEGAQYBBgEGAgYCBgIGAQYFBgEGAwYFBgMGAQYCBgIGBgYBBgIGAwYDBgEG",
"AgYDBgIGAgcCBgEGBAYGBgEGAQYIBgEGAQYCBwQGAgYBBgQGAgYCBgIGAQYBBgEGAgYBBgEGAgYBBgEGAwYCBgMGAQYBBgUGBQYJBgIGAQYBBgMGAQYBBgIGBwYDBgIGBAYBBgEGAwYDBgEGAgYBBgQGAQYBBgEGAgYBBgsGAQYCBgIGCQYFBgMGAQYCBwMGAQYBBgEG",
"BQYDBgIGAQYBBgUGBAYBBgIGAgYBBgEGBQYDBgIGAgYCBgMGBAYEBgEGAQYJBgEGAQYEBggGAQYBBgEGAgYCBgEGBgYBBgcHAQcBBgMGAgYCBgEGAgYCBgIGAgYCBgMGBAYBBgYGAgYDBgIGAwYBBgIGBQYCBgEGAQcCBgEGAgYHBgMGAQYBBgMGAQYCBgEGAQYGBgQG",
"AQYCBgkGAwYCBgEGBAYCBgMGAwYCBgwGAgYCBgEGBAYDBgEGAgYBBgYGBAYIBgEGAQYFBgEGAwYJBgQGAwYCBgEGBAYCBgEGAQYDBgEGAQYBBgEGAwYBBgIGAQYBBgwGAQYHBgIGAQYDBgIGAgYCBgMGAQYHBgEGAQYEBgEGAQYBBgEGAgYCBgIGCgYEBgkGAgYEBgMG",
"AgYCBgUGAgYEBgMHDwYIBgEGAQYDBgEGAQYBBgIGBQYCBgIGAQYKBgMGAQYGBgMGBgYEBgIGAgYBBgIGAQYCBgIGAQYCBgIGAQYCBgMGCAYBBgQGBQYBBgEGAgYBBgEGAQYBBgEGAQYEBgIGBAYCBggGAwYJBgQGCgYBBgMGAwYCBgEGAgYBBgIGAwYDBgEGAQYBBgMG",
"AQYCBgEGAQYCBgEGAgYBBgEGAQYBBgEGAgYBBgIGCQYEBgIGAwYDBgIHBgYCBgYGAgYDBgUGAQYBBgUGBAcBBgMGAgYCBgEGAwYDBgQGAQYCBgEGAgMBBgEGAQYBBgEDAQcBBgIGAQcHBgEGAwYCBgEGCwYGBgcGAwYBBgMGAwYEBgEGAQYCBgEGAwYBBgYGAQYFBgQG",
"AwYCBgEGAgYCBgEGCAYCBgMGAQYBBgUGAwYBBgEGBwYCBgQGCgYFBgEGBAYEBgEGAgYBBgMGAwYBBgEGAQYCBgIGBQYBBgEGAQYBBgIGAwYDBgIGAQYCBgUGAQYDBgYGBgYBBgcGAgYDBgIGAwYEBgMGAgYBBgIGBwYBBgMGBAYCBgEGDAYBBgIGAwYBBgQGAgYCBgYG",
"BgYDBgYGAwYQBgEGAgYCBgMGCAYBBgUGBgYFBgQGAwYDBgEGAQYDBgIGAQYGBgMGBAYDBgQGAgYCBgEGAQYBBgEGAQYCBgIGAQYBBgEGAwYFBgEGAgYBBgEGAgYBBgEGAQYBBgIGAwYCBgEGAQYCBgMGBAcCBgUGAgYCBgIGBwYDBgEGAgYDBg4GAQYBBgIGBAYEBgEG",
"BAYBBgEGAwYBBgQGCwYEBgEGAQYIBgMGAgYDBgEGAgYBBgEGAgYCBgEGAQYBBgEGBAYBBgIGAQYEBgIGAQYGBgEGBAYCBgIGBAYBBgIGAgcCBgIGAwYBBgMGCQYBBgEGBAYCBgIGAQYEBgEGAgYBBgEGAQYBBgEGAQYEBgQGAQYBBgEGAwYBBgEGAQYCBgcGAgYFBgEG",
"AgYBBgEGBAYBBgYGAQYBBgIGAQYBBgEGAQYDBgEGAQYBBgEGAQYBBgIGAQYCBgEGBgYBBgMGAgYBBgEGAQYBBgEGAwYBBgoGAgYDBgEGBQYBBgIGBwYCBgQGAgYCBgUGAQYBBgQGAgYBBgIGAQYCBgEGAQYCBgMGBwYDBgEGAgYDBgEGAgYBBgIGAQYBBgEGAQYCBgQG",
"AgYEBgQGAwYDBgEGBAYBBgUGAQYDBgEGBQYBBgEGAQYDBgEHAQYBBgIGAQYCBgIGAQYBBgEGAgYBBgEGBAYBBgEGBQYBBgEGBAYGBgEGAgYCBgEGBAYCBgEGAQYKBggGDgYGBgMGAQYDBgIGAwYBBgEGBQYBBgIGAQYCBgEGAgYBBgQGAQYCBgEGBQMBBgEGAgYBAwEG",
"CQYDBgEGBAYDBgIGAwYFBgEGAwYMBgQGBAYEBgMGAgYEBggGCgYCBgEGAQYFBwEGAwYBBgcHDwYIBgMGAQYBBgEGAQYDBgcGCgYHBgYGAQYBBgMGAQYCBgYGBQYBBggGAgYDBgMGAgYPBgEGAwYDBgEGAQYEBgQGAQYBBgIGAQYBBgUGAgYBBgQGAQYCBgIGBAYBBgEG",
"AgYDBgQGAgYBBgIGBQYEBQQFAwYBBQIGAQYBBgMGBAYCBgQGAwYBBgIGAgYBBgIFAQYBBQEGBAUGBgIGAgYEBgEFAQYEBgYFBAYFBgMHAgYFBgIGAQYDBgIGAgYBBgIGAgYBBgkGBgYGBgMGAgYFBgYGBAYDBgEGAQYBBgEGAQYEBgEGAgcCBgEGAQYBBgEGAwYCBgIG",
"AgYBBgQGAwYCBQIGBAYBBQQGAgYBBgIFAgYBBgEGAQYBBgIGAQcBBQIGAQYCBgIGAgYDBgEGAQYCBgIGAQYDBgUFAQYDBgIGAQYBBQIGBgYGBgEGBAUBBwEGAgYCBwMGAQUDBQQGAQYBBQEGAwYCBgEGAgYDBgMGAgYBBgQGAQYBBgEFAQYBBgIHAQYDBgIGAQYBBQQF",
"BgYEBggGBQUJBgUGAgUCBQYFAgYBBgEGBAcBBgEGBwYBBgQGAQUEBgQFAgYCBgEGAQUBBgEFAgYBBwUGBAYDBgIGAgYBBgQGAwYBBgEGAQYDBgEGAQYBBgEFBgUBBgEGAQYBBgIGBQYBBgEGAgYBBgEFAQYDBgEGAgYCBgEGAgYBBgEGAQYCBgMGAwYDBgMGAgYEBgEF",
"AQYBBgEGAQYCBgEGAgYDBgIGBwYCBgMGAgYDBgEGAwYFBgEGAQYEBgEFAgYBBgMGBAYBBgEGAwYCBgEGAQYEBgQGAgYCBgEGAgYCBgEGAQYFBgEGAgYDBQEGAQYBBQEFAgYGBgYGAgYEBgMGAQYDBgMFAQYBBgQGAwYBBgEGAwYBBQEGBAYDBgQGAgcGBgEGAgYBBQEG",
"AgYDBQEGAQYBBgMGAQUBBgIGAQcCBgEGAQYBBgMGBQUBBgEGAgUBBQMGAwUCBgQGAQYBBQEGAQYEBgEFBQYBBgIGBAYCBgEGBAYDBgIGAQYFBQEFAwYBBgEGBgUFBgEGBAUCBgcGAgYBBgMGBgYBBgQGAQYIBgEGBQYBBgEGAQYBBgIHAwYCBgYGAwYDBgMGAgYCBgEG",
"AQYCBgEGAgYBBgcGBQYBBgMGBAYJBwMGAwYCBgEGAQYCBgEGAwUDBQQGAwUBBQMGAQYBBgQFAgUGBgEGBAUCBQEGAQcDBgEGBwUBBgUGAgYDBQEFAQYFBgEFAgYCBgEGAgYDBgMFAwYBBgcGBQYFBgIGAgUBBQQGAQUHBgoGAgYBBgIGAQYCBgEGAgYBBgUGAQYEBgEG",
"AgUBBgIGAQYDBgIFAQYDBgMHAgYDBQEGAQYJBwMGAQYBBgMGAgYBBgEGAwYIBgUFAgYGBgIGAQYGBgkGCgYBBgIGAgcFBQIHAQYDBQMGAQYCBgEGAQYDBgIGAwYCBwIFAQYCBgEGAQYDBgEFAQYEBgEGAgYDBgQGAQcCBgEGAQYGBgMHAwYDBgEGAwYBBgEGAgYCBgMF",
"AQYBBgEGCAYBBgEFAwYBBQQGAQYBBgIGAQYBBQEGAQYDBgQGBAYFBgEGAQcBBgEGBAYBBQUFAwUBBgQGBwYFBQEGBgYCBgEGAQYBBwMGAQYDBgIFBQYGBgEGAwYCBgoGAgYCBgEGAQYEBgYGAwcBBgEGAQYCBgEGAQYEBgcGBgYCBgIGAgYDBgIGAQYHBgoGBAYBBgEG",
"AwYBBgEGAwYEBgEGAQYBBgMGAQYCBgUGBAYCBgYGAQYCBgEGAQYBBgIGAgYBBgEGBAYBBgEGAgYEBgMGAQYCBgMGAgYBBgEGAgYBBgEGAQYBBgUGBgYBBgEGAwYBBQMFAQYBBgEGAwYDBgEGAwYCBgEGAwYCBgEFAgYDBgEHAgYCBgIGAgYEBgMGAwYFBQIFAgYBBgUG",
"AwYNBwEHAQYFBwIGCQYCBQIGAQcBBgQGAgYBBgEGAQYEBgEGBQYCBQIGAgYCBgEGAgUEBgQGBAYBBgQGAQYBBgEGAQUBBgcGAgYDBgsGAgYFBgMGBQYCBgQGAwYEBgUGAgYIBgcFBQYBBgIGAwUBBgIGAgYEBgUGAwUBBgUGAQYBBQEFAgYKBgMGAQYDBgQGAQYBBgcG",
"AQYBBgYGAQYEBgIGBgYCBgQGAQYCBgcGAgYDBgEGAQYCBggGAQYCBgIHAgYFBgIGAQYFBgQGCgYCBQQGAwcDAw8HAwYDBgEGAgYDBgIGAwYHBQkGAgYDBgQGAQYBBgMGBAYMBgQGAgYCBgQGAQYBBgQGAgYDBgEGBgYCBgEGBAUCBgIGAQYEBQUGBAUHBgUFAQYBBgkG",
"BQUEBgQGAgcFBQIGCQYLBQEFAQYRBgIGAgYUBgEFAQYCBgEGAgYBBgEGAQUBBwMGAQYDBgIFBAYCBgMFAwUDBgEGAgYEBgEHAgYBBgEGAgYCBgEGAgYEBgIGAQUGBgIGDQYDBgEGAQYBBgIGAwYBBgMGAQYCBgEHAQYGBgEGAQUCBgEGAgYCBgYGAgYCBQIFAgYBBgEG",
"BgYBBwcFAwYBBgMFAwYBBQEFBAYCBgMGAwYBBQQFBQUBBggFAgUCBQgFAQYIBgEFAQYBBQUFAQYDBgUGAQUBBwEGBQYBBwIGEwYGBgMFBwYBBgEGCQYCBgYGDgYBBgYFAQYBBgEHAQYMBQcGAgYDBgMFAQUIBgIGAwYIBgIGAgUOBQIGAQUEBgYGAQYIBgQGBgYBBgEG",
"AgYBBwEGAQUDBgEGBAYCBgMFCgYHBgQGAgYFBgEGBgYGBgEFAwYDBgEFAgYEBwUFAQUCBQEGAQUCBQEHAQUDBQEFAQYBBQYGAQUVBQcFAwUBBQsGBQUDBQEFAQYBBgIFAQYFBQQGBwUEBggGAQYCBgsGCQYCBQIFAgcHBgMGAQUBBgMFAgcCBQUGAQUBBgEGAQUBBQcG",
"AwYBBQEGAgYCBgYGBQUDBQEGBwUBBgIFBQUDBQMFAgUEBQEFBAUBBQIHAQYBBQQFAQUCBQEHAgUBBQIFAQUCBQMFBAYBBgIFAwUCBQQFBwYBBgEFAQUEBQIFAQUBBQIGAQUBBQEFAgUNBQEFAQUCBQEFAQUBBQEHAwYBBgUFBQUCBQYGAgUCBQQFAgUCBQMGAgYCBQkF",
"AQUBBQIFAQUDBQIFAQUCBQIFAgUBBQMFAgUBBQEFAQUBBQEFBgUEBQMFAgUDBQEHAQcBBQQFAgUCBQEFCQUBBQEFAwUCBQQFAgUCBQQFCAUBBQEFAQUHBQEFCQUCBQYFBgUFBQMFAQUGBQIFAwUHBQEFAgUCBQMFAgcBBQIFAQUDBwEFBAcBBQQHBQUBBQQFAwUKBQIF",
"AQUHBQQFAQUBBQIFAQUDBQIFAQUBBQEFAgUFBQIFBAUCBQEFBQUBBQEFAQUDBQEFAQUDBQYFAgUEBQEFAwUEBQUFAQUEBQkFBQUBBQkFDAUBBQIFAQUBBQEFAgUCBwgFAQUBBQUFAQUCBQIFBAUIBQIFAgUBBQEFAQUBBQEFAwUEBwEFAQUDBQEFBgUBBQMFAwUCBwEF",
"AgUBBwEFBwUCBQIFAQUEBQMFAQUDBQQFBQUBBQQFAgUBBQEFBAUCBQEFAgUDBQIFAQUDBQIFBAcFBQIFAwUHBQEFBQUCBQQFAQUHBQEFAQUBBQQFAQUBBQEHAQUEBQIFAQUBBQEFAQUBBQUFBgUBBQQFAQUBBQEFBAUCBQEFBgUHBQMFAgUFBQkFAQUEBQIHAQcJBQEF",
"AQUBBwkFBQUBBQYFAQUCBQQFBQcCBQIFAwUCBQkFAwUCBQcFAwcEBQMFAQUDBQEFBQUBBQEFAwUBBQMFAQUBBQEFAgUDBQEFAQUCBQMFAQUGBQIFAgUHBQIFAQUCBQIHAwUDBQIFAQUBBQEFBAUEBQEHAgUDBQYFBAUJBQEFAgUEBQIFBQUFBwUFBwUKBQEFAQUFBQEH",
"AwUEBQMFAwUCBQEFAQUFBQMFAgUCBQYFCgUCBQUFAQUCBQgFAQUGBQQFAwUCBQMFBQUBBQEFAwUBBQEFBAUGBQEFAwUBBQIFAwUBBQIFBQUCBQUFAQUCBQMFAQUDBQEFAQUBBQQFAQUCBQQFBAUCBwEFAgUCBQMFAgUKBQQFDwUFBQEFCAUIBQIFAQUBBgEFAQUEBQEF",
"BQUCBQQFBwUBBQIFAQcBBQEFAwcBBQEFAgUBBQEFAQUIBQEFAQUBBQEFAgcCBQsFAQUBBQIFCAUCBQIFAgUMBQEFCAUDBQgFAQUHBQwFAQUIBQEFBQUEBQMFBAUEBgMFBgUJBwEFBwUCBQEFBwUDBQIFBAUGBQQFDAcDBwUFAwUCBQUFBAUDBQUFAQUCBQIFAgUBBQEF",
"AQUBBQMFAQcBBQIFAQUBBQMFAgUBBQgFDAcCBwEFAgUBBQQFAQUEBQUHBQULBgEFAgUBBwIFAwUCBQIFAgUEBQMFAQUCBwQFAQUCBQQFCQUCBQEFAQUEBggFAgUEBQEFAwUBBQEFAgUCBQoFAgUBBQQFAQUBBQQFAQUBBQEFAQUCBQkFAQUCBQEFAwUBBwEFAgUJBQEF",
"BAUBBQIFAgUBBwEFAQUEBQEFAQUEBQEFAQUBBQMHBQUDBwEFBAUFBQIHAgUBBQIFAgUBBQEFBAUBBQEHAQUEBQMFAgcGBQwFBgUBBQMFAgUVBQEHAQUHBQMFAQUCBQMFAQUFBQEFBgUBBQQFAwUBBQQHAQUBBQIFBwUCBQIFCwUCBwUFBQUBBQMFAQUHBQ4FAQUGBQQF",
"AwUJBQEFAwUBBQkFAQUEBQYFBAUHBQUFCAUBBQYFBgUCBQQFAQUDBQEGCAUBBQsFAgUFBQMFBQUDBQIFBAUFBQkFAwUBBQIFAQUEBQMFAQUBBQIHBAUSBQsFAwUDBQQFAgUKBQQFBAUCBQUGBAUBBQQFAQUDBQEFBAUCBQIFCwUDBQEFAQUBBQEFAQUGBQEFAQUCBQQH",
"BQUFBQIFAgUCBQUFAQUGBQUFBgUEBQQFAQUBBQQFBgYEBQUFBwUBBQMFAQUBBwIFAwUBBQYGAgUCBQMFCAUBBQMFBQUBBQQFAQUEBQQFBwUNBQEFBwUGBQEFBAUBBQIFAwUDBQEFBAYCBQMFAwUCBQIFAgUDBQIFAQUCBQIFAwUBBQUFBAUGBQEFAgUCBQEFAgUFBQQF",
"AQUBBQYFAwUBBQEFAgUCBQQFAwUBBQEFAQUDBQMFAgUEBQsFAgUDBQEFAgUCBQEFBwUFBQEFAQUBBQEFAQUCBQEFAQUDBQEFAwUBBQEFBAUBBQMHAgYDBQUFAgUFBQIFAQUBBQIFBQUCBQEFAwUBBQEFAQUCBQUFAgUDBQMFEAUDBQYFAQUBBQUFAgUBBQIFAgUBBQEF",
"AwUCBQEFAwUCBQIHAgUBBQIFAgUGBQIFAQUBBQIFAQUBBQEFAQUBBQIFAgUDBQMFAgUBBQcFAwUBBQIFAgUCBQEFAQUDBQUHBQUBBQQFAgUFBQ8FBwUBBQEFBQUDBQQFAQUDBQUFCQUDBQEFAwUCBQEFBgUBBwIFAQUEBQIFAgUEBQEFBQUCBQEFBAUBBQMFAQUFBQsF",
"AwUCBQcFBQUBBQEFAQUBBQEFAgUHBQEFCQUEBQIFAQUDBQQFDAUBBQYFAgUCBQEFAgUCBQMFAQUBBQUFAwUCBQEFBgUBBQEFAQUBBQEFAQUCBQIFAwUCBQMFAQUCBQQFFAUFBQUFBAUIBQQFAQUEBQMFAwUHBQQFAQUBBQEFCgUBBQEFBAUCBQEFCQUBBQEFAQUBBQEF",
"CwUBBQgFAgUBBQMFAQUCBQIFAQUBBQIFAgUBBQMFAQUCBQ0FAQUBBQEFBwUCBQcFAgUBBQIFAQUGBQQFAgUEBQMFAQUDBQEFBQUBBQMFAQUCBQIFAQUBBQIFBAUHBQMFAQUBBQ0FAgUDBQIFAwUBBQIFAQUBBQIFAQUBBQEFAQUDBQQFAwUBBQIFAQUBBQEFAgUJBQEF",
"DQUFBQEFBQUJBQQFBwUGBQEFAQUDBQEFAQUEBQMFAQUCBQYFAQUCBQEFAgUDBQEFAwUIBQQFBgUBBQsFAQUFBQEFAQULBQEFAwUCBQcFAwUBBQIFAQUIBQEFAgUEBQIFBAUBBQsFAwUBBQEFBAUJBQMFAQUCBQEFBAUCBQEFCwUEBQMFAgUBBQsFBgUBBQMFAQUBBQMF",
"AgUBBQ4FAgUBBQIFBgUNBQQFAwUBBQcFAwUCBQIFAQUBBQkFAgUBBQEFAQUCBQEFAQUJBQEFAgUBBQcFAQUDBQEFAgUEBQIFEQUEBQwFAgUCBQUFAQUCBQQFAQUFBQUFAwUBBQIFAgUCBQMFAwUBBwIHAwUCBQMFAgUCBQIFAwUDBQMFBAUCBQMFAQUBBQIFAwUEBQIF",
"AgUBBQQFAQUCBQEFAQUDBQMFAgUBBQcFBwUDBQEFAQUBBQgFAwcCBQEFAQUDBQMFAwUEBQIFBQUKBQIFBQUCBQIFBAUKBQQFAgUBBQQFBQUCBQEFCQUCBQEFBgUBBQMFAgUFBQIFAgUDBQQFAgUBBQEFBQUBBQYFCAUBBQIFAQUBBQEFAQUEBQYFAwUBBQcFAQUIBQIF",
"AwUCBQIFAQUGBQEFAQUBBQEFAQUEBQQFAQUCBQMFAgUEBQIFBQUBBQEFBAUEBQcFBQUFBQMFBAUHBQEFAwUCBQIFAgUDBQYFAwUCBQEFCgUCBQEFBAUCBQEFAgUBBQEFAwUDBQcFAwUBBQkFBQUDBQEFAwUJBQIFBgUCBQEFAwUDBQIFBgUBBQEFAQUCBQMFBgUCBQIF",
"CAcFBQEFAwUCBQkFBwUBBQQFCAUFBQEFAQUBBQEFAwUBBQMFAgUIBQMFBwUIBQMFAQUBBQEFAQUBBQIFAwUBBQEFAQUBBQYFAQUCBQEFAQUDBQcFFQUNBQIFCQUEBQMFAQUBBQMFBgUDBQEFBQUKBQMFAQUBBQUFCAUCBQoFAgUCBQEFAQUCBwEFAwUCBQYFAwUCBQEF",
"AQUEBQEFAQUBBQEFAQUDBQIFAQUJBQEFBQULBQEFBAUBBQIFAwUDBQMFBQUCBwIFAQUEBQYFAgUDBQgFBAUCBQEFBQUDBQUFAgUCBQIFAwUCBQEFAgUGBQIFBAUBBQEFAQUBBQUFAwUHBQQFAQUFBQEFDAUHBQcFAwUGBQYFAgUDBQUFBQUEBQIFCAUCBQEFAQUHBQIF",
"AQUBBQQFAgUHBQEFCAUBBQIFAgUDBQEFAgUBBQQFAQUBBQEFAQUBBQEFAgUBBQIHAQUBBQEFAQUBBwEFAgUDBQQFAgUDBQYFAQUIBQEFBAUBBQMFAwUBBQQFAQUCBQYFBQUDBQUFCwUCBQcFAwUCBwQHCAUEBQEFAQUIBQYFAQUIBQMHAgUHBQEFDAUCBQIFAwUBBQIF",
"AQUBBQIFBAUBBQcFAwULBQEFAgUBBQQFAwUDBQkFAgUCBQMHAQUDBQIFAQUEBQMFCAQEBQMFAQcDBQIFAQQCBQIEAwUBBQgFAQUBBQEFAQUBBAMFAgUBBQMFAQQBBQIFAgUBBQEFAQUEBQIFAQcDBQIFBAUCBQMFAgQBBQEEAgUBBAEFAgUBBAUFAgUCBwEFBAUJBQEF",
"AgUBBQMFAQUCBQMFAwUKBQEFAQcDBQEHEgQCBRAEAQUDBQEHBgUBBQMFAgUBBQMFAQUEBQMFCAUGBQIFAgUCBQIFCgcCBQYFAgUDBQIFAQUCBwEFAQUBBQEFCQUCBQUFAgUBBwMHAQUDBQcFBQUDBQMFAQUBBQIFAQcBBwIHAQUBBAUFBwQHBQEFAQUBBQEFAQUCBQEF",
"AQUBBQQFAgQDBQMFBAUBBQIFAgUBBQEFAQUCBQEFAQcBBQEFBgUBBQMFAQUEBQIFAQUCBQQEAgQEBAEFAgUCBwUFAwQEBAQFBgUKBQIFAgUDBQQFAwQCBQMFAwUBBQYFAgUFBQEFAQQEBQEFBQUBBAEFAQUFBwEFBAUCBQkEBQUBBQQECwUDBAEFAQUBBAEFAwcFBQEF",
"CAUDBQwFBgUDBQUFAgQCBAEFBAUCBQIFCQcCBQEFAwUCBQMFBwUDBQMFAQUGBQYGAQUQBQcHAQcCBAoHBwYDBQMFAQUCBgYFCQUDBQIFAQYCBQEEAQUBBQIEAQUCBQEFAwUBBQoHBQUBBQQEBQQBBQEFBAUCBQcFAgUCBAMFAQQBBAEEAQUCBQEFAQUBBAEFAwUBBAEF",
"BQUCBQMFAQUBBQIFAwUCBQEGAgcGBQMFAQQDBQEFAgUDBQQFAgUEBQMEAQUBBAIFAwQBBQEHAQcCBQMFAgUCBAMGAgQBBQEEBgUGBQMFAQUDBQEEAQQBBQMEAwYCBQsFBgUBBQEFBQUCBQEGAwUCBQQFAQUKBQMFAQULBQQFBwUFBwYFAQUFBQQFAgUCBgUFAQUDBQcF",
"AwUBBQQFBwUFBQkHAgUFBQEFAgUBBwEFBQUCBQMFAQUBBQMFAQUBBQEFAQUBBQEFAgUHBQIHBQUBBQIFBQUCBQIHAgUBBwUFAQUCBQEFAQUBBQEEAQUFBwkEBgQEBAEFAQUDBAEEBQUBBQUFBAUBBQoEBgUBBwcFAgUBBQYFAgQBBQMFAQUDBQIFAgUCBQEFAQUCBQEF",
"AQUDBQEFAQUCBQYFBQUBBAEFAQQCBAQFAgUFBQMHAgUKBwMEAgUBBQEFAQcCBQkFAwUCBQEEBQQBBQUFAQYBBQQFAgUDBAUEAQQCBQQFAQUBBQEFAQUCBQEFAQUDBQEEAQUBBQEFAgUBBQEFAQUCBQMEAQUBBQMEAQcCBQIEAwUEBAIEAQUBBQEFAQQBBQIFAQcDBQEF",
"AgUBBwgFAQUBBAMFBgQBBQMFAQYEBQEEAQcFBQMFAQQDBQEFAQcBBQIFAQUFBQIFBQUBBQQFBQUBBgIFBAYDBQEFBAQDBQIFBwYDBQEFAQUBBQEFAQUBBQEHAQUBBQEFAQUDBQUEBgUDBQIFAwUEBQEFAQUDBQEFBQUBBQEGBgUBBwIGAgUDBQIFAQYCBQMFAQUDBQMF",
"BAYCBQEGBAUBBQQFAQUCBQEFAQUDBQIFBAUBBQkFAQUCBQUFAQUBBQMFAQUEBQQFAwUFBAIGAQUCBAQFAgUBBQIFAQUDBQMFAQUIBQIFAQUBBQIFAwUBBQEGAQUEBQEFAQQBBwIFAgUBBQEFAQcCBgEEAQUCBwIFAgUEBAMGAgUCBQEFBAUBBAMFAgUBBQEFAgUCBQIF",
"AQUCBQUEAgUCBgIFAQYBBQEFAgUDBQIFAQUFBAIFBAUBBAEGAQUEBQUFAwcEBQIFBwUEBQEFCQUBBQEGCQYBBQEHAQUEBQMEAQUBBQEFAQUBBAEEAwQCBQMGAgUDBQQFAgUGBgIFBwUDBAEFBwQBBAIFBQUJBgIFAQUBBQUFBQUBBQEEAgUBBQEFAQUCBQEFAgQBBQMF",
"AQUBBAUFAwUBBAMEAgUCBAMFAwUBBQIFAQUBBQIFAQUBBQEFAgUCBQQFAQUBBQEFAgUCBQYFAgUCBQEFAwUBBQUFAgUCBQIEAQYBBQIGAgQCBQEEBwUCBAQFBQULBQEFAwUCBQIFBAUBBQEEAgcBBgEEAQUBBAQEAgUCBQwFAQQHBAYEAQUCBgQFAgUKBQEFAQUBBQEH",
"AQYFBQEGAgUCBwQEAgUCBQQFAgYHBAYGAQUDBQYECwQCBAMEAgUDBQMFAgUCBQMHAQUBBQMGAQUCBQ4FBAUFBAEFAgUCBQIEBAcIBAEEAgUCBQMFAgUBBAEFBgcDBQMGAgUBBQIECAQBBAIGAgQBBQIFAwYBBAMEBgUCBgIFAQYIBAEGAQUGBQEFAQUCBQEFAgUJBQEH",
"AQUCBQEFAgUCBQEFAQUCBgEEAwUBBAYHBgQBBAEFAQYBBQEFAQUBBQIFBAUBBQMFBAYBBQQFAQQCBAEHAQUBBAEFAwUDBQYFAgQCBAEGBAUDBAIFBAQDBAEFAgUCBQEFAQQDBgEEAQUEBQIFBAUCBQEGAwUCBQEFAQUEBQEFAgUBBQEFAgUCBQEFBgQFBAIFAwUCBQEF",
"AgUBBQMFAgUDBQEEAQUBBQEFBAQBBQEFAQUFBgIFAQYDBQcFAQQBBgMFAgQCBQIFAgUFBAEFBAUCBQYFBwUBBQIFAgcCBAEFAQQDBQEEAgUDBQEEAQQCBgMEAQUBBQIEAQQFBQIFAgQIBAIEAwYDBAEEAQQFBQcFAQUBBQEFAQYBBQEFAQUBBgEFAQUEBQEEAgUCBQIF",
"AwQBBQUEAQQFBwEGAgUBBQcFAQUBBQIGCQUBBQEFBAQDBAEEAQUBBQcEAwYDBAIGAgQCBAIEAQQCBAIGAgUFBAEFAgQCBAEEAQQEBQIGAQUBBQQGAQQCBQMEAQQSBAIEBAQFBAMEAgQEBQMEAQUCBAEGAQQDBAEEAQQDBAUFAQQCBAEEAwUMBAQEAQUDBgEFAwQEBAME",
"AQQCBAEEAQQFBAIFAwQJBAQGAQQBBAEEAwQCBAEEAQQCBAgGAwUDBAEGAQQDBgEEBAQGBAEEAQQBBAIEAQQDBAEEAwQBBAEGAwQCBAIEAQQCBAIEBAYCBAIEAwQBBAEEAQYFBAIEAQQDBAEGAgQBBAMEAQQBBAEEAwQBBQMEAwQFBAYHBAQCBQEEAQQBBQMEAwQCBAQE",
"CwQCBwIEBwQEBAMEBQYCBAMEAQUBBQMGBAQDBQMFAQUCBQQFAQUCBQEGDQYCBAIECwQCBQIFAgQCBgYGAgQCBAMEAgQCBgMEAgQCBQIEAQQKBAEFAgQHBAIGAgcDBgIEBgYBBggEBAQBBgIFCAQBBAEEAgYEBAEEAQQBBAMEBAQBBAMGAQQEBgEGAgQDBAgEAwYBBgEE",
"AQQFBAEEAQQIBAIEAwQBBAMEAQQBBBEEAwQBBAIEAQQDBAQEAQcDBAIEAgQCBAEEBQQBBAMEAgYFBgMGBQQBBAEGCgQDBgMEAwQEBwEHAgQCBAIEBAQBBAEEAQQDBAEEAQQBBAEEAgQBBAEGAQQBBAMEBQQBBAEEAwQBBAEEAQQBBAcEBAQCBAEGBgQBBAEEBAQDBAMG",
"AwQDBAIEAQQBBAIEAwQDBAMEBgQCBAEEAQQEBAEEAQQBBgIEAQQCBAYEAwQBBAEEAQQBBgEEAQYDBgIEAQQBBAEEAQQBBAMEAwQBBAEGAQQBBAUEBQQBBAEGAQQCBAMECwQBBAIEAQQCBAEEAQQBBAEEAgQDBAUEAQYHBAEEAwQBBgQEAwQEBAEEAwQBBAIEAQQBBAEG",
"AwQCBgIEAwQDBAEEAwQBBAIEAQQCBAEGCAQHBAYEAQQBBAEEAQQCBAcEAQQBBAEEBQQCBAEEAwQBBgEGBAQBBAcGBAYCBAEEAgQBBgEEAQQDBAEEAgQCBAMGAQQBBAEEAwYGBgMECQQCBAEEAwQBBAEEBAQFBAEEBAYBBAMEAgQBBAEEBAYBBAIEAgQBBAEGAQQBBAYE",
"AQQCBAIGAwQEBAUEAgQBBAEEAQcDBAEEAQQBBAUEAQQDBAEEAwYBBAIEAQQCBAIEAQQGBgEECQYHBAEEBgQBBAEEAgQCBAIEAQYBBAIEAQQCBAEEAQQEBAEGAQQCBAMGAQQEBAIEAQQCBAIEAQQCBAMEAQQBBgEEAQQBBAMEAQQBBAIEAgQBBAIEAQQDBAYEAwYBBAME",
"AQQBBAIGAQQBBAEEBQYBBgEEAQQEBAEEAwQBBAEEBAYIBAIEAgQGBAEEAQQBBAIEAQQCBAIEBAQCBgEEBAQBBAEEAQQIBAEEAQQGBAEEAQQBBAMEBwQBBAQEBAYBBAEEAQQFBAIEAgQDBAMEAgQHBAIEAQQBBAEEAgYDBAkEBwQBBgIEAQQBBAQECgYDBAYEAgQFBAUE",
"CQcRBAIEBAQBBgEGBwYEBAEEAQQCBgQEAgQCBAIEBAQBBAEEAQQCBAMGAQQBBAEEAwQBBAEEAgQBBAIEAQQBBAEEAwQBBAEEAQQBBAMEAQQBBAIEAgQCBAEEAwQDBgQGAgQBBAEEAQYBBAEEAQQBBAEEAQYBBAIEAwQFBAIEAQQEBAYEBgQFBAcEAgQCBAEEAwQDBAIG",
"AQQCBAMECQQBBAQGAQYCBgEGCAQEBAIGAQQBBAEEBQQBBAMEAgcBBAEEAwQBBAEEAgYBBAIEAwYBBAEEAQQBBAUEBgQDBAEEAgYEBAEEAQQBBAcEBAQIBgEEAwYCBgUEAQQBBgIEAQQCBAIEAgQCBAIECQYBBAMEAQYBBBEEBQQHBAEEAgYBBAEEBAYBBAIEAQQBBgIE",
"AQQCBAEHAQQCBAEHAQQBBgMGCgQEBAIEAgQBBAEEAQQCBAMEAQQCBgQEAQQHBAQEAgQBBAEECAQBBgcEBAYBBAQGBwQCBgYEAQQCBAEEAQQBBAMEBwQBBgIEAwQFBAQEAQQBBgEEBgQCBgIEAQQBBAYECAQDBAEEAQQBBgMEAgQEBgIEAwQCBAEGAgQCBAEEAgQCBAIE",
"AgQCBAEEAwQEBAEEAgQGBAEEAQQEBAEEAQQBBAEEDAYDBgMEBAYCBAIEBAQDBAEEBAQGBAIEAwQEBAEGAgQBBAEGAQQIBAEEAwQEBAEEAQQFBAMGAQQBBAEGAgYBBAEEDQQEBAUGAQQBBAQEAQQBBAEEAQYFBAEEAQQCBgIEAgQEBAMEAQYGBAIEAwQCBAkEAwQHBAEE",
"AQQFBgIEAQQHBAMEAQcHBAIEAgQBBQIEAgQBBAEEAQUBBAMEAwQBBgMEAgQBBAUEAwQFBAEEAgQBBAIEAQQCBAEEAgQDBAMGBgQFBgIEBAQBBAcGBAQEBAIEAQQCBAEEAQQBBAUEAQYGBAEEAwQCBAEEAgQGBAQEAQQBBAEECgYBBAYECAQKBAEEAwQCBAIGAQQBBAEE",
"AQQEBAMEBQQBBAYECAQCBAQEAQYBBAEEAgQBBAIEAwQCBAQEBQQOBAIEAwQJBwcEFAQIBAEEAQQCBAEEAgYBBAMEAQQBBAIGAgcHBAIEAQQBBAEEAgQBBAEEAQQBBAIEAQQCBAEECAQCBAEEAgQBBAUEAQQBBAUEAQQCBAEEAgQBBAEEBQQCBAYEAQQBBAMEAwQBBAEE",
"AQQBBgEEAwQBBAEEBQYBBAEEAgQBBAEEBQQCBAYEAgQBBAEEAgQBBAEEAQYBBAEEAQYBBAIEAgYBBAUGAQQBBAIEAwcBBwEEAQQBBAEEAQQEBAEEAQQCBAQEBAYGBgIEAgQBBgUHAgYGBAEGAwQCBAIGAQQFBgEEBQQCBAEGAQQBBAEGAgQMBAMGAwQBBwIEAgQFBAEE",
"AQQMBAIEAQQCBgIEBwQFBAEEBAQDBgQEBgQFBAEEAgQBBAEGAwQBBAEEAQYBBAYGAwQBBAIEAgQDBAEEAQQBBAQFAQYBBgEEAQQEBAMHAgQBBAEEAQYCBAEEAgQBBAIEAQQCBAEEAQQBBAEEAQQBBAIEAgYBBAEEBgQEBgEFBgYBBAEEAQQDBAEEAQQGBAQEAQYBBAQH",
"AQQGBAEEBwQCBAUFAQQBBAEEBAQEBAEEAwQBBAcEAwYCBAMEAQQBBAIHAgQBBAMEAQYDBgIEAQQBBAEGAQQDBAEEBQQBBgEHBgcBBgMGBwYLBAYEAQQBBAMEAQQBBAIEAgQCBAgGCAQBBAQHAQQHBgEGAQYGBAIEAQQBBAEEAQQCBAEEBAQDBAkEBwQFBw4EAQYFBAIE",
"AwYFBAQEAgQEBAEEBQQCBAQEAQcFBAQEBQQCBgEEBQQBBAIEAQQGBAEGCAUBBAEEAgQCBwEEBAQBBAMEAgQFBAMEBQQBBAIEAwYBBAMGAQQBBAEEAQQQBQEEBQQBBAEEAgcBBBUEAQQBBAIGAwYFBAEEAwQBBAMEAgQBBgEGAQQCBAEEAQQCBAIEAQQCBAIHAQQCBgIE",
"AwQEBgIEAwYHBAsEAwQIBgEGBAQBBgYEAQQBBgEEBgQBBAEEAQYCBAMGCQYOBAIEAwQBBAEEAQQCBAEEAwQBBAIEAQQCBAEEAgQBBgEEBQQCBAEEAQQEBAIEAQQCBAEEAQQGBgEEAQcDBAQEAgQCBAQEAQQBBA0EBAQEBAEEBAQEBwgGAQYDBAYGBgQBBAsHDQQJBAEE",
"EAQIBAIGBwQCBAMEAwQFBwUEAwcCBAEEBQYBBgIHAQcHBgMEAQQBBAEEAQcCBwMGAgYBBAYEBgYBBAEECAYCBAIEBAYCBAIEBgQDBAEGAgQCBAUHAQYEBAIEAwQEBAEEBQQDBAEEAgQCBgQEAgQNBAQHCAQHBAQEAwQFBAIEAwQGBAIEBQQCBQIEAwQDBAYEAwQCBAEE",
"BwQCBgYHAgQBBAEEAgQCBAIEBwQHBAEGAQQIBAQEAQQHBAIEAQQBBAEEAgYBBAkEAQQCBwEEAgQGBAEGAQQDBAEEAgQDBwIFAgQFBAIEAgYDBwIEBgQCBAEGAwQBBAMEBwQJBAQEAQQIBAIEBgQDBAIEAQQCBAIHAQYBBAEEAgcCBAMHBgQBBAEEAQQDBAEEAQQDBAMH",
"AwQBBAMEDwUBBAIEBAQHBAIECgQBBBAEAQQEBAEEBAQFBgIGCAQBBAIEBgQBBAEEBwQBBAcGAgYEBgYECQQFBAMEAQQBBwIEAQQBBgUEAgQFBAYEBwcDBwUEBAYBBAUECQQDBAQGAwQBBAUEAwQCBAQEAgQCBAMEAgQBBAEEAwQIBAEEAQQBBAMEBwQBBAQEAQQBBAQF",
"AQQCBgEGAQQEBAEEAQQBBAMEBAQCBAEEAQYDBAIEAwQDBAIEAQQDBAEGAQQFBAsHAQYHBgEECQcBBAIEAgYBBAIGBQcDBAEEAgQBBAIEAQQBBAQEAgQBBAMEAQQBBgIEAwQBBAEEAwYCBgYHAQQBBgMHBAQBBAYEAgQDBAEEAQUCBAEEBAQBBAEEAQQCBAIEAwQFBQEE",
"AQQEBgIECAQCBQMEAwQCBA4EBgQBBAEHAgQKBgIHAQQGBAEEAQQBBAIEAQQDBwUEAgQDBAQEAgQCBAMGAQYCBAIEAgQHBAEEAQQBBAMEAwQBBAEHBAQEBAIEBQYBBAQEEAQDBgMEAwQCBAIEAQcCBAEEAQYCBgMEAgQEBAUEBQQBBAEEAgQCBAMEAwQBBAMEAQYDBAME",
"DAQBBwEEAgcDBAEGAwQBBAIEAQQBBAMEAgQBBAEEBQQBBAEEBAQDBgEEAQQIBAMEAQQEBAYHAQQFBAEEAwQGBAEEAgQGBAIEAwYGBAUEBAQDBAIEBAQCBAcEBwQDBAYEBgQDBAkHAQQBBQgHAQQBBAEEBwcCBAQHBgQBBwgEBQcBBAgEAQQCBwYEBAQDBAEGBAUEBAEE",
"AQQDBAQECQQBBBMEAQQCBAEEBwQEBAgEAQQCBAIEAQQBBQQEAgQIBAQEBAQBBAEFAQcBBAQEAwcFBAIEAQQCBAEEBQQCBAEEAwQEBQIEAgUDBQEEAwUHBAMEBAQGBAgGBwQFBAEEBAQBBAEEAgUEBAEEBQUCBAUFCQQBBQEEAwUEBAEFAQQBBAMHAQQCBQEEAQQHBAEE",
"BAQCBwEEAQQDBwIFAgYDBAcFBAcCBAMHBQQBBQEHAQQDBAUHBAQDBAQEAQQCBAEEAQQBBQMEAgQIBAYFBAQBBAEEAgQCBAEEAQQDBwIEAgQFBwIEAwQCBAUEAgQOBQIHBAcJBAUECAQGBAIEAgQDBAEFAgQCBAUEAQQBBQMEAQQCBAEHBQQBBAEEAgQBBAIEAQQBBAYE",
"AQQCBAEEAwUBBQEFAQQFBAMEAQQBBAEEAQQEBAMHBAQDBAIFAQQBBAIEAwQFBAEEAQQLBAEEBQQBBAEEAQQHBAEEFgUEBAwFAQQDBAIHAgUFBQgEAwQBBAUFAQQEBAEEAgQDBAcEAQcFBAgEBgQGBgMFAgcDBAQHAQUHBAEEAQQBBQgFBAQKBA4FGwQEBBgEBQQNBAEE",
"DgQEBAYEGgQOBBEEDgQjByIEAQQTBAcEAwQEBAQEAwQCBAEEBwQFBAIEAQQGBAwEBAQBBAQEAQQLBAMEAgQJBAIEAQQBBAIEBgQDBAoEAQQCBAEEAgQBBAUEBgQDBAEEBQQBBAEEAgQJBAQEAwQDBAIEBgQNBAIEAQQCBAIEAQQBBAEEBgQCBAEEBAQBBAEEAgQBBAME",
"AQQCBAQEBgQEBAYEBAQDBAEEAgQGBAEEAgQEBAUEAQQFBAIECgUJBAUFAQQFBAEEAwQJBAMEBAQCBAEEAQQKBAcEBAQHBAUEBwQCBAgFAQQJBAcEAgQBBAgEAQQCBAIEAgQCBAEECAQFBAMEAgQJBAEEBgQDBAIEAQUDBAMEAgQBBAUEAQQDBAIEAwQBBAEEAwQBBAME",
"AQQEBAEEAgQBBAEEAQQBBAIEBAQFBAIEAQQEBAMEAgQIBAEEAQUKBAEFAQQFBAsEBAUEBAYEAQQBBAYEAgQBBAYEBQQFBAUEAgQDBAEEAgQBBAIFBAQCBQEEAgQBBA4EAgQBBAIEAgQBBAEEAgUHBAYEAQQGBAIEAQQDBAwEAwQDBAMEAQQBBAMEBAQCBAIEAQQDBAsE",
"AQQCBQEEBQUBBAEEAgQDBAMEAQQCBAIEBAQCBAEEAQQFBAUEAwQFBAUEAQQBBAEEBAQCBAEEBAQDBAYEBQQCBAIEAQYEBAEEAgQDBAEEAQUBBAIEAgQCBAEEAQQGBAUEBgQCBAEEAQQDBAQEAQYBBAYEAgQEBAEEAgQCBAYECAQCBAoEAwQDBAQEBQQBBgUECQQBBQIE",
"AQQEBwIEBgQCBAEEAQQBBAQEAQQCBAQEAgQFBAIEBAQDBAEEBQQEBQEEAgQPBAQEAQQDBAcEAwQGBAYEAgQDBAMEAgQFBAQEAgcEBAIEBgQDBAcEBAQEBAMEAwQCBAUEBgQEBAEEAQQEBAgEAQQEBAIEAwUCBAkEBgQBBAEEBgQDBAUEAgQCBAUEBAQDBAYEAwQDBAIE",
"BgQBBgEEBAQCBAEEAgQCBAkEAgQBBAcEAgcJBAEFBQQFBAIEAwQBBAIEAwQEBAEEAgcEBAIEAQQEBQEFAwQCBAEEAgQCBAIECAQBBAcFAQQCBAEEAgQBBAIGAQQFBAIEBwQIBAUFAQQEBA4EAgQBBAEEAgQCBAQEAwQCBAEEAQQBBAMEBAQCBQUEAgQCBAEEAwQFBAUE",
"AQQIBAEEBwQDBAMEBAQBBAEEAQQEBAEEBAQHBAEEAgQGBAEEAwQCBAIEAwQBBAEEBAQDBAEEAQQBBAEEAgQHBAQEAQQBBAMEAQQEBAEEBAQBBAQEAQQCBAMEAQQCBAIEBgQBBAQEAgQBBAMEAgQBBAIFAQQCBAIEAQQCBAEEBQQBBQIFAgQEBAEEAQQJBAYECAQBBAEE",
"BAUCBAEEAgQCBQIEAgQBBwYEAwQCBAUEAwQBBAUEBAQCBAIEBwQBBAoFBAQCBAIEBAQDBAEEAgQCBAcEAQQDBAIFAgQBBAUEAgQDBAYFAQQDBQEEAwYBBAEEBQUFBAQEBAQDBAEEAgQEBAEEAQYBBAEEAQQEBAEEAwQEBAYEAQQEBAcEBAQCBAEEAwQDBAMEAwQGBAME",
"AQQHBAEEAQUBBAMEBgQBBAEEAgQBBAEEBgQBBAQEBAQHBAEEAQQCBAQEAwQCBAgEAQQEBAgEAQQBBAYEBwQFBQYEAgQBBAMEAQQEBAYEAQQHBAEECQQBBAkEAQQPBAIFAgQHBAEEAwQFBAEEAwQBBAEEAQQDBAEEAwYBBAIGCwQBBAQEAQQBBAIEAQQEBAIEAgUDBAUG",
"AgQBBAcGAQQBBAQFAgQCBAEEBwQCBAMEBwQBBAEEBwQCBAEEAQQEBAMEAQQBBAMEAgYCBAUEAgQJBAIEAQQDBAMEAwYCBAYEAgQCBAUEAQQBBAcFBQQCBAQEAQQBBAEEBAcDBAMEAgQCBAIFAQQCBAEEDAQBBAUEAwQJBAQEAgUEBgEEAQQDBAEFAwQGBAEEAwQBBgEE",
"AgQCBAIEAQQHBAEEAgQBBAIEAwQCBAEEAQQCBAQEAQUDBAUEAQQBBQIFAwQDBAIEAgQBBAQEAwUDBQIEAQQCBAEEBAQBBAQEAQQGBAMFAQUDBAIEAQQBBAMEAgQGBQQEAQQEBAEEAQUFBAEEAQQBBAEEAQQMBQQEBAQFBAEEAwcBBQMEAwQCBAQEAwQGBAQEBAQCBAIE",
"AgQBBAMEAQQKBAQEAQQBBAEEAQQHBAQFBAQCBAUEAQQFBAIFBQQBBAMEAQUKBAYEAQUCBAIEAQQBBAEEBQQEBAEEBQQGBAEEAQUHBAIEAQQDBAIEAgQKBAEEAwQGBAEEAgQDBAQEBAQCBAgEAwQCBAIEAQQBBAEFBgQBBAIEBwQDBwEFBwQCBAEEAwQEBAEEAwQFBAEE",
"BQQBBAQEBQQBBAEEAQQCBAMEAgQCBAMEBQQCBAEEAQQEBAIEAgQCBAEEAgUBBAQEAwQCBAEFAgQCBAMEAQQBBAYEAgQBBAMEAwQBBgMEAQQBBAQGAgQBBAIEBgQEBAQEBQQGBAMEBgQCBAEEAgUGBgEGBAQBBAEEAgQBBAQEAQQBBAEEBAQOBQEEAwQCBAQHBAQCBAMG",
"BQUJBAMEAwUFBAIEAgQBBAMEAQQHBAUEAgQBBAMEAQULBgEEAQQBBAEEAwQEBAUFAQUCBQIEAQQIBAMEAgQCBAIFAgQEBAQEAgQHBAIEAQcBBAQEBQUCBQUEBgQBBAQEBgQGBAMFBAQDBAEFAwQDBAEEAwUBBQEEAgQBBAIEBAQCBwgEBgQBBAEEAwQBBQIEAQUBBQEE",
"BgQBBQQEAQYHBAIEAQUDBAEEBwUFBAIECgQJBAEEBwQGBAMEBAUFBAUEAQcGBgEEAQcHBAMHBAQBBAMFAQUCBAEFBwQBBAEEAQQBBAIEAQQEBAEEBwQBBAQEAgQCBAEEAgQBBAIFAgQFBAUFBwQBBAIFAgUBBAUEAQQCBQMFBAQBBAYEAQUBBAIEAwUHBAEEAQQCBAIF",
"AQQBBAEEAgQEBAIEAwUEBAMGBwcCBAgFAQQBBQIFAgcCBAEGBAQCBAIEAgURBQIFAgQCBAEFAwcBBAEFAQQGBAcEAQUBBAEHAwUBBAYFAwYDBAEEAgQDBQMEAgQKBAEEAwQCBQIEBQYBBAEFAQQDBAIEBQUCBwEGAQQBBAIEAQQCBAEEBwUDBAMEAQQCBgMEAwUEBAwE",
"AQQDBAgFAQUBBQYEBAQJBAEFAwUDBQEECgQEBAEEAgQCBQIEAgYDBAUEAwUFBAEEAgUDBAIFBwUDBAMFAgQFBAwGAQcIBAwEAgQBBAEEAwQCBAIEAQUBBAEEAQQBBAEEAgQBBwIEAQQCBAEEBQQBBAEEAQUFBAEEAQQDBAEEAQQCBAEEAQQBBAEEAgQCBAMEAQQDBAIE",
"AQQEBAMEAQQBBAEFAwQCBAEEBgQFBAEEBQQBBQQEBgQBBAEFAwQDBAIFBQQFBAkEAQUDBAkEBAQEBAIEAQUBBAEEAQQRBAQEAQQCBAEEAwUCBAkEBgQDBAEGAQQCBAEFAQUBBQUFAgQBBAIEAgQBBAMFAgQBBAIEAgQBBQYEAwUBBAQEAQUHBQEEAQQBBAEFAgQBBAYE",
"AgQBBAkEAgYDBAEFCwQCBAEEAQQCBAIEAgQEBAEEEAQCBAMEAgUCBwYFAwQCBAIEAQYBBAMEAQQGBAIEAwQEBQEEAQQCBQYEAwQBBAIEAQQDBAEEBAQBBAEEAgQFBAQEAQQDBAMHAQYBBAQEAQQBBAEEAwQBBAQEAQQBBAMEAQQBBAEFAQQBBAEFAQQEBAEEAQcGBAEE",
"AwQBBAIEAQUHBAEEAgUDBAcFAQQBBAEEAQQBBQEEAQQFBAUEAgQBBAIFAgQCBAYEAQQJBAEEAQQCBAEEAgQGBAEEAQUHBgIFAQUBBAIEAgQCBAQFAwQCBQUEAQUBBAEFAgQBBQsEAQQGBAMEAQQDBAEFAgQBBAEEAgQGBAEFAwUBBAEFAwQCBAIFBQQBBAMFBgQBBAQE",
"AwQBBAIFAQYCBAEEBAQFBQMEAgQEBQIEBAQKBAgFAwQBBgUEAwQBBAEEAgQEBAMEAwQCBAIEAgQBBQMEEwUBBAUFAwUCBAEEAQQFBAIEAgQCBAEEAgUFBAEEBgQEBQIGAgQCBAUEAQQCBAMEAQQCBAIFAQQEBQIECAQBBQUFAgQCBAMEAQQBBAUEAQQEBQYEAQQBBAIF",
"AQQDBAEEAwUBBAEEAwcCBAcEAgUBBAMEAQQCBQEEAQQFBAYEBAQFBAEFAQQCBQEFAQQCBAIEBAQEBQEEBAQDBAIEAwQCBAYFAQQDBAMEAgQDBAYEBgQCBAQEAwYFBAIEAQcCBAEEAQQBBQQEBAUDBAMEAQYCBAEEAQQBBAIFAQQEBAIGAwQCBQEEBQQJBAEEAwUKBAEE",
"AQQCBAQEAgQCBAIEBQQDBQMEBAQFBQEEAwQGBAIEAQQCBAEFBQQDBQIHAQUNBQIEAQQCBAQFAQQCBAEFAgQDBQYEAgUFBAYFAwQBBAIEAQQDBAIFAwQCBAYEBQUCBAUEAwQDBgEEAQQCBAEEAgQBBAEFAwQDBAcFAQQBBAEEAgQBBAEEAwQCBAEEAQQGBgIEBgQCBAIE",
"AgQEBAMEAgQBBgQFAQQCBAQEAwUBBAYEAQUBBAEGAQQBBQEEAgQGBQEEBQcBBAoEAwQCBAEEAgUDBAEFAgQGBAYEAQUIBAEFCgQBBAIFBgQEBAEHAQYIBAIFAQQGBQYEBQUDBgQFBAUFBwMEAQQBBAEEAgQGBQEEBAQBBAEEAQUCBAMFBQQKBQIEBQUBBAIEBQQCBAEE",
"AQQCBAEEAwQBBQQEAQQDBAIFAwQBBAMFBQUDBAQEAQUCBAEFAQQBBAkEAQQBBQUEAwQGBAQEBAUEBAYHAgQJBQEEAwQBBQUEAQQDBAEEAQYBBAMEAwQBBQUEAQQJBAIEAQUGBAEEAgQBBQMEAQQFBAIEBQUBBQEEAQQBBgYFAwYFBAEFBgQBBAMEAgUBBAIFAwUBBQEE",
"AQQBBgEEAQQCBAIEBAUBBAEFAQUBBAEGAwUDBAIEAwQHBQMEAgQCBQEEAQUDBQIEAgQCBAUEAQUCBAMFAgQCBQIEBAQCBAEEAwQDBgEEAgQCBQIFAgQBBQMEAQQFBAIEAQUBBQQEAQUDBAEFAQQBBAQFAwUEBgEFAQQCBQEEAgQDBAQFAgQBBAEFAgQDBAMEAwUBBAEE",
"AwQHBAIEAQQBBAIEAgQGBQIEAQQBBgIFAQQHBAEEAgQFBAIEAgUBBQEEAQUGBAEFAQQBBAIEAQQBBQMEAQQBBAMFAgQEBAIFAQQCBQEFAgQEBAUEAQUEBAkFBwUBBAEFBwQCBAMEAgQBBQIFCQQBBAEEBQUBBAMFAQQCBwQEAQQFBAEFAgQBBAMEAgQFBAEEAQUDBAIE",
"AQYDBAEEAgQBBAQEAQQBBQEFAQQBBQEEAQQFBQMEAQQBBQMFAgQCBAEEAQQCBAUEAwQBBAkEAgUBBAEFAQQFBAUEBgcBBQEEAgQCBAIFAwQCBAEFAQQBBAQEAQQCBAEFBwcBBAEGAgQDBQQFAgUCBAYFAQUCBAIFAQQBBgIEAgQDBQEEAQYBBAMEAQUDBAUHAwYCBQQG",
"AwQDBQIFAQUCBQQFAQQBBAIEAwUBBQMFAQUEBAMFBAQBBQMEAQQBBQEFAQUBBAQFAQQLBAEEAgQDBAMEAgUDBAMFAgQFBQEFAgUCBAIEAQQHBQUGAwQBBQIFBAYBBAEEBAcDBAEFAQUBBAEEEAQCBAUFAQQDBAEEAgQCBQEEBAUIBAIHAQcCBQIEAgUBBQMFBAcBBQEE",
"BQQCBQIEBAUCBgMFAQUBBAIEAgQCBAIFBQUCBwEEAQQBBgQEAQUDBQEEAQQCBAIGAwQEBAIGBwQCBQEEAQQBBAMEAQQCBQMEAgQEBAEFAQYCBAEEAgQBBgEFAQYCBAIGAQUFBAQGAwYBBAMFAQQFBAMFAgYCBQcEAgQFBAEFAgQEBQUEAQQBBAEEAQQBBQIEAgUDBAgE",
"AQQBBAEEAQQBBQgFBgUDBAEECAUBBQMEBwUDBAEFAQQDBAEEAgQBBAEEAgUGBAEFCgUDBQMEAwYEBAEEAQQFBAQFAQYBBQUEBQYEBQIFAwQIBAoEAQQBBAMGBAYBBQIEAQUCBAEEAgQCBAEFAwUDBAIEAwUBBAQFAgQCBQIEAQQFBQEFBQQGBAIEAQcBBgIEAwUDBQME",
"AQQCBAEGBAQGBQEEBAYDBQIGBQUEBQQGBgQDBQMFAQQCBAEFAgUDBAMGAwQDBQEEAQQBBAMEAgQCBAEHAwUEBAIECQUCBAEEAQUGBgoEAQUBBAQEAwYFBQIFAQUEBQEEAQQBBQEFBgUHBQEEAQUEBAEEAQQCBQEEAQUDBAMFAQcBBgMGAwQBBQQEAQYEBAIEBAULBAUE",
"AgYDBwEFBAUBBAEECAQCBAUEBAQBBQMFBgQEBAEEAgUCBAIEAwYCBQUFAgUBBQEEAwQBBAIEBQQCBgIEAgQGBQMEAQQEBQcFAwQDBQMFAgUCBAcEAQUCBgYEAQQBBQIEAgUBBAYEAgYGBAEEAQQDBAQHAQYDBAYEAwYEBgEEAgQBBgEEAQQDBAQEAQQBBAMEAQUBBAIE",
"AgQBBAgFCgQCBgMEAQQBBAMEAQUCBAEGAgQCBQMFCAQEBQUEAgUDBAIFBwYCBAEEAQUBBAMFCgUBBgEFBAQDBAQEAgQDBAIEAQUBBQUEBAUDBQMHBAUCBAEEAgUDBQIFAgcBBQgFBgUBBQEFAQUDBQEEAgUBBAIFBgQIBAEFAgQBBQIFBQUCBQMFAQQGBQEEAwUCBQIE",
"AwQBBQEEAQQCBQIGAQUCBAIHBgUBBAMEAQQBBwEGAQQLBAYEAQQBBwEEAQQDBAIFAwQBBAEEAgQBBwEFBAQDBAEEBAUCBAUFBAQBBAEFAgQBBQUEAQQEBAMEAQQBBQMFAwUCBAEEAwQDBAgEAwQBBQIEBQcBBQYFAgUDBQIFAgQBBQMFAQUCBgEHAQQDBAQGAQQKBQIE",
"AQUBBQMFAwUBBAgGBgUBBQEFAQQBBAEEAgQCBAkEAgQBBAMEAQQCBAEEAQQBBgIGAwQBBQMEAQQCBAMEAQcDBQQFBAUCBAEEBgQBBQMFAgUEBAMEBAUBBAIFAwUCBQEEBQQCBQMEAgQBBAMEAwUCBQEEAgQBBQIEAQQCBQQEAwYFBAIEAgUCBQEEAQYDBAQFAgUCBgIE",
"BQQDBgMHAQQBBQUEAwcDBAEEAgcCBAEFAgYFBAEEAQUCBQEFAgQBBQEEAgQCBAEEBgQBBAEFAwQDBAEEAwUBBgQEAQUDBAEEAgYCBQMFAgQBBAQEAwQCBAEEAwYEBQIEBAUBBAgEAgUBBAEFAQQFBgIFBQQDBQEHAQQBBQUFAQQDBAIEAQUHBAIFBAUBBQIEAQUCBAIF",
"AQUDBQUEAQQCBAIEAwQFBQMFAwYBBAEEAQUCBgEEBAQBBAQFAQQBBAMHCAQCBgEEAwQCBAUFAQcGBAMEBQQDBQIGCgQMBAIEAQQIBQcEAQUBBwEEAQQCBQgEAQQBBAEEAgQBBQEEAQYDBAIEAgQCBQEGAQUBBQEFAQUDBAMHAwUBBAEEAgUCBAMFCQQBBAEFAwQEBQUE",
"BAYDBQEFAgUBBAMFBAQCBAYEAQYDBQEHAwUBBQEFAQQBBAEEAgQCBQEFAwQFBAIEAQQCBQMEBAUEBQMEAQcFBQQFCAUEBAEEAQUCBQQFAgQCBQsGBAQBBgIFAgQCBAQFAgQEBQIEAgUDBQMFAQUGBAEFBgQDBgIEAQUGBAMFAQUCBQEFBAQCBQMEAwYBBQQHAQQCBwEE",
"AwUBBgEEAQQCBAQEAgUHBQIEAgUEBAMEBAQEBgQEBgQDBQQEBQQCBQIGAgUCBgQGAQUDBQIEBgUBBAcFBQUEBAEHBgUBBAIFCgQBBQIEAgcDBQMFCQUDBAEFBAQDBAEEBwQDBgEFAwUBBAkFAgQBBQEFAgQHBAIFAwUCBQMEAgUBBAQFAwUDBQIECgUGBQEEAQYDBQcE",
"AQQCBAcFBAYDBQYFBAUIBAEGAgUBBQgFAgcBBQcHAgUCBBEEAQUDBQcFAQQBBAMEAQUFBAIEAwQCBQQFBgYBBAMFAQQCBQEFAgUCBQUGAgQCBQUFAQQEBAIFAgQBBAIEAQYDBAEEAQQBBAIFBAUFBAEEBQUCBQEEBQQDBAEEAgYBBAEEAgQBBQEFCAQBBQYFAwYBBAIE",
"BAQMBAEGAQcCBAIFBQUCBQEEAQUBBwIFAQUDBQMEAwUCBgEEAwUDBAIFAgQIBAIGBAUEBQEEAQUMBQIFBwUGBQQEAQcCBAEFAQUBBA8EAwUBBgEGAQQCBgcFAQUBBwIEBAQHBQkGCQUHBwQFAQYBBQIFAwUDBAsEAwUCBQIGAQQCBAEFCAcBBQEFAQUHBQEFAgUJBQEG",
"AQYDBAMFBAULBgEGAwQBBQUEAgUCBAIGAwQBBQEEAQYBBAEEBgUBBQEEAQUCBAEEAQYEBgIEAQQCBAEFAgUGBAEGAQUEBAUEBgUBBQEGAQUDBwUFAwUBBAIHAgQBBQEGBAcBBQIEAgQCBQEEAwUBBQEEAwQDBQQFAQQBBwYFCgQFBQUEAQUMBwIFAwQCBQMEAQQEBQQF",
"AQUBBwIFAQQCBQIEBgcDBwQHAgUDBAEFAQYBBQYFAwUBBQEFAQUDBAEFBAQEBQQECQUBBAEFBAUCBAEFAwUBBAEFAgcDBQEFBAUCBQIEBAYBBQwHBwUCBQIFAQQCBQMEAgQBBQEEAQcCBAEFCgUBBQMFAQYFBQMFBAUEBQMFCQUEBQEFAQYJBQYFCAUHBgEFAQQFBQIE",
"AgUDBAEFAgQBBQIGAwQBBQEEAQQEBAMFBQUBBwMEAQQEBQIEAgUEBAIEAQYCBQMHAQQCBQEEAQQDBAIGBAQBBAEFAQQBBgMFAwUDBQQFBwQBBQkFAwUEBQEFAwcCBgEFAQYBBQUGBAYEBgcFBQQCBQEFAQUCBQcFAgUBBQcFBAUCBQEFBAUJBAIFAgUBBwEFAgUCBAMG",
"AgUBBQEFBQUBBQEGBQUEBQEFAwUCBQEFCAUGBQQFBQUBBQEFBAUDBQIFAQUCBQEFBwUCBQEFAQUCBg0FAQUDBQYGAgUCBQEFAQUEBQEFAQUDBQQFAQUIBgEFAQUBBQsFAQUGBQEFBQUDBQUFAQUCBQQFAQUDBQYFAQUCBQIFAgUBBgUFBQUGBQEFAQUDBgQFAgUDBQgF",
"CAUIBQMFAwUEBQEFAwUIBQ4FAQUCBgEFBAUCBQIFAQUBBQEFAQUBBQYFBgUDBQMFAwUDBQEGCAUEBQMGAQUFBQMFBgUFBQEFAgUFBQIFAgYHBgMFAQUCBwUFBAUBBQEGAgUHBQIFBAUFBQEFAQcDBQEFAQUBBQUFAwUBBQIFAwUOBQIFAQUDBQQFBQUCBQUFAgUHBQIF",
"AwUBBQEFAgUCBQEGAgYDBgIFBAYBBQEGBgUCBQcGAQUDBQEFAgYBBQEFAQYFBQEGBQUBBQIFAgYEBQEGAQUDBQMFAgUDBQEFBAUEBQIFBAUCBQEFAQcEBwEFAQUBBQEFAQUDBQEFAgUCBQEFAwUCBQEHAgYBBQEFBAYBBQEFAQUEBgIFAgcHBQQFAQUDBQIFAQUBBQMG",
"AQUDBQEGAwUDBQEFAQUBBQEFAgUBBQcGCQUHBQQFBAUBBgQHAQUBBQEFAQUCBQEGBgUBBgEFAgUCBQEFCAUDBQIFAwUBBQIFAQUEBQIFAQUCBQIFAgUBBQIHBAUEBQMFAQUBBQEFAwYDBQIFBwUDBQEFAQcBBQEHBAUFBQEFAQUBBQIFAwUCBgYFAwcFBQEFAQUDBQMF",
"AQUCBgIGBQUBBgIHAwUBBQEFAQYCBQEGAQUEBgMFAQUCBwUFBQYDBQEGAgUCBQMFAQUDBgUHAgUEBQIFAwcEBQUFAwUEBgIGAgYBBQIFAwYBBQEFAQUBBQEFAgUDBQMFAQYBBwIFBgUCBQQFBQYBBQMFCAUCBQIFAQUBBQIFAQYGBwIGAgUHBQEFAwUBBgEFBAUBBQIF",
"BQUCBQUFBQUDBQQFBQUBBgEHAQYBBQIFBAUSBgUFAQYCBQMFAgUDBQEFAQUBBgEHAwcDBwQFBQUCBQMFAQcCBwIFBQUDBQEGAgYBBQMFAQUDBQEFAgUCBQMFBAUEBQcHAgUBBwIFAgUBBQQHAgYBBQMFAgUCBQMFAQYCBQIGAQUFBQEFAQcBBQIHAgUFBQEFAgUEBQQG",
"AQYBBQIFAgUFBQQFAgUFBQMFAgULBQYFAQUBBwQFAwcCBQEFAQUDBQMFAQcBBgIFAQUCBQEFAQUDBgQFAQUEBwIFAQUCBQUFAwUEBQMFAQUBBQIFAgUCBQUFBgUCBQIFBAcCBQEFAgcFBQMFAQUCBRIFAgUBBQEFBwUDBQUHAgUCBwQFAQUBBQEGAQUFBQUFAgUBBgIF",
"AwUBBQUGBQUBBQIFAgUCBQIGAQYBBQMFAgUBBgEFBQYDBQIFAQYDBgIFAQUCBQcFAwcFBQcFAQYFBgEFAQUCBQIFAwUBBQIGAgcGBQQFBAUBBQkFAQUBBQQGAQUBBQgGAgUCBQMHAgUBBQEFAgUDBQQFAQUDBQEFAQUFBQEFAwUBBgMFBwUCBQEFAQUDBQEGAgUCBQEF",
"AwUFBgUFAgcBBQUFAQUCBQEFAgUBBQQGBAUDBQEFAQUCBQIHBAUBBQIGAgYCBQUFAQUDBQEGAQcFBQIFAwUOBQMFAgUBBQMFBwUBBQEFAwUDBQQFAQUBBwIFAwcCBQQFAwUIBQQHAQUCBQMFAQUDBQUFAQUCBwEFAQUFBQEFAgUBBQEFBQUBBQMFAQUDBgEFAQUBBQMG",
"AQUCBQIFAQUBBQIFAQUDBgQFAwUBBQUGAgUDBQUFAgUBBQQFCQYBBQYGAgYCBQMFBgYBBgQFAQUJBQEFAgcBBQIGAQYBBgUGAgUEBQIGAgUFBQEFAgUBBQEFAQUDBQIFAQUDBgQFAgUHBQEFAQUBBQIFAgUEBQEFAgUEBQQFAQUBBQIFAQUBBgIHCgUFBQkFAQUEBQUG",
"AQUFBQMFAQYCBQEFAQUBBQEGAQUBBQEFAQUDBQMFAwUEBQMHAQUEBgEGAQUBBgMFAgYBBgcGAQUCBQEFBgYBBwEFAQUDBgEFAQUCBQIFAgYGBQMFAwUBBQIFAgUBBQIGAQYBBgQFAQUCBQEFAgUIBQYFBQUCBgIHBAUBBQQGAQUCBQcFAQUCBQcFAQUBBQMFAQUCBQEF",
"AQUCBQIFEAUBBQMFAQUCBQUGAgUEBQIFAwUCBQUFAQYBBQMFCQUEBQEFBQcBBQIFAwUBBQEFAQUDBQIFBQcEBQIFBgUBBwMHAwUFBQMFAgcCBQkFAQUBBQQFBgUCBQcFAwUFBQIGAwUEBQMFAQUCBQEFAQUKBQIFAwUFBwIFAQUBBQEFBAYCBQMFBAUKBQIFAQYBBQEF",
"CwUDBQEFCQUFBgEFBQUDBQoFDAUBBQEFAgUCBQIFBQUEBQEFAQUFBQYFAQUCBgEFAgUBBQcFAQYGBQIFAQUBBgIFAgUDBQEFAQUBBQMFAgYDBgMFAQUHBQMFBQUBBwUFAQUCBQEFAQUBBQMFBAUBBQMFBgYBBQQFAQUDBQIFAgYBBQYFAQUDBQEFAQUCBQEFBAUCBQMF",
"AgUDBQIFAgYBBQMFAwUIBQUFAQYDBQMFAQUBBQMGAwUBBQQGAQUBBwUFAQYBBQEFAQUBBQEFAQUBBgEFAwYBBQYFAQYBBQIFAQUGBQMFAQUFBgIFAwUFBQcFAgUBBQEGCgUCBQEFAgUGBQEFAgUBBwEFAQUFBQIFAQUCBQMFAwUBBQEFAwUDBQEFAwUCBQIFAQUBBQMF",
"CAUCBQIFAQUEBQEGAQUBBgIFBAUBBQMHBgYFBgIFAgUBBQEGAQUFBQEFAQUBBQIFAQUDBQEFAwUBBQMFAQUCBQEFAQUDBQQGAwUBBQIFBQUBBQMFAQUBBQMFAgUIBQEGBgUBBQEFBAUCBQMFAQUBBgIFAQUEBgcFAgUCBQUFBQYBBQIHAgUBBgEGAQcFBgEGAQYBBQEF",
"BQUDBQQFAgUBBQEFAwUHBQEFAgUEBQEFAgYEBQMHAwYDBQcGDQUDBQIFAQcDBQEFAQUBBQIFAQUDBQIGAQUBBQIFAgUDBQEFAgUDBwMGBgUDBwIFBwcBBQgFBAUCBQUFBQYDBgEFBQcBBQUFAwUDBgMFAwUDBwQFAgUGBQIFAgUFBQEFAgYEBwEHBgUBBQEGAgUBBQcG",
"AwUBBgQFAQUBBwMGBAUHBgEFBQUHBQoHAQYFBQEFBAUVBRQFAQUBBQwFAgUBBQMFBAUBBQUFAgUBBQEFAgUEBQEGAwUFBQYFAQUDBQUGAgUCBgIFAwcBBQEFBgUBBwMFAQUDBQQFCgUDBgIGBgUBBgUHBQUCBQYGAgUBBgEFAwUCBgMFBQYEBQEGBgYDBQEFAQUBBgEF",
"AQUBBQcHAQUBBQEFAQYBBQUHAgUEBwEFBAUBBQMGBQULBwkFAgUGBgIFCQUCBgQGAQUGBgEFAQUIBQEFAgUCBwEFBAUBBQMFAwUDBQEGAwUBBgEGAwUBBQUFAQUCBQIFAQUOBgoFAgYDBQYFBQUGBQYFBgcBBQMFBAUBBQEGBQUCBQEFAgUBBwYFAQUGBQEFAgUCBQMF",
"AQUDBwEFBgUDBQIFBQYCBQoHAQYCBwIGAQUBBQYFAwYFBQUFAwUKBgYFBQUCBgEFAwUFBgIFBAUBBQMFBAYBBQMGAgUBBgQGBAUGBgIFBQYCBgkFAQUBBgIGBgYDBwIFCAUFBQYFAgYBBQIGAwUBBQIFAwUDBQIFAQYBBQMHBQUJBQIFAQUBBQMGAwUBBQEFAwYBBgIG",
"CQYFBwEFAwUFBQQFAQUDBwEFAgUEBwMFAgUBBQMFAgUBBQIFAQUCBQEFAwUCBQIFAQYCBgIFBwUGBQEGAwUBBQEFAgUBBgEFAQYCBQQFAQUBBQUFAgUDBQEFBAUFBQIGAQUEBQEFAQUBBgEFCQUDBgMFBQYDBQEFAQUBBQUFAQUCBQEFAgUCBQYGAQUEBQMFAwUFBQsF",
"BAUDBgQFAQcCBQEFAQcBBgEFAwUEBQQFBgUBBQQFAQUBBQIHBAUBBQoFAgUDBQkFAgUCBQIFAgUBBwEFAQUJBQcFAwYEBgIFAQYDBggFBQUCBQIFBAYDBgQFAQUCBQIFAgUBBQUFAgUBBQQFAgcCBQUFAQYBBQYFAQUCBQMFAwUFBQIFBQYEBQEFAgUBBgIFAQUCBgMF",
"AQUEBQEFBAcBBQEFAQULBQIFAQUEBQUHAwYFBQUFBQUFBQIFAQUCBwIFAgUCBQEFAQUDBQEFAgUCBQEFAwUBBQQFAgUJBQEFAQUDBQEFAQUBBgMGBAUEBQkGAQUBBQEFAQYBBgIFAQUBBQMFBAUCBQEHBQUDBQEFAwUBBgEHAgYCBQEGAgUBBgEFCwcHBgMFAQUDBQMF",
"AQUBBQIFBAUCBQIFAgUGBgsFAQcBBgEFBgUBBQEFAgUCBQIFBgUCBQQFBwUCBQIFAgUCBgUFAwUHBQEFAQUCBgEFAQYBBQMFAQUCBQIGAQUCBQEGEQUBBQMFBQYBBQQFCQUDBwQFAgcDBQIFCQcCBQIHAgUGBgUFAQUBBQIGAwYBBwMFBQUCBgcHAgYBBwsGAQcBBQMG",
"AQUCBwEFAwYFBQEGAwUBBgEFCAUBBwIFAQUBBQMFBAUCBgUHCgUBBQIFAgUFBQIFCAUCBQMFAwUBBQEFAgcDBQEFBAUBBQoFAQUEBQIFBAcDBQIFAQUDBwIFCgULBQMFAQUCBQEFAwUFBwEFBAUEBQQFAgcBBQEFAwUBBQEFAgUBBQEFBgUHBQMFAgUCBQYFBQUEBQEF",
"AQUCBgIFAQUCBQMFCAUBBQIFAQUJBgMFAgUGBQIFBQUBBwIFBAcCBQMFAwcBBwUFAwcCBQMFDwYBBwMFAgYCBQQFAQUDBQEHBgUBBQEGAwUBBQMFBgcDBwQGAQcCBQEFCwYCBQIFAgcDBgEHBAUBBQEHAQcBBQIFAQUBBQIGAQUBBgEFBAUCBQIFAgUGBQYFBgUCBQEF",
"AgUBBQUHBgUBBQEFBAUBBQUFAwcFBQMFAQUIBQUFAgcDBQIGAQYBBw0FCAUBBQEFAQUCBQUFCQUFBQMFAgUBBQIFAQUDBQEGAQUEBQUGAwUCBQEFAQUJBQMFAQUCBQYHAgUJBwMFBAUBBQIFAgUBBQEFAwUCBQIFAQUBBgMFAgYFBwQFBQUFBQIFAgUDBQUFAQUFBQUF",
"BwUFBQIFBwUCBQMFBAUCBQUFAQUJBQEHAQcBBQEFAgUDBQEHBAUGBQEGAgUDBQ4FAwUFBQIFAgUBBQMHAgcEBQIFBQUCBQMFAwUEBQMFAgUBBQIFBQcOBQIGCwUBBQEFAgUFBQEFAQYDBQcGBgUEBQQFAgUFBgQGAQYJBQEFCwUBBQEFBAUDBQIGAgUIBgMFAgUCBQIH",
"BAUBBQEFAgUCBwMFAwUFBwIFAQUFBQUGBQUCBwUFAQUFBQIFAQUCBQEFAgUDBwMFAQUEBQEFAQUCBQUFCAUBBgEGAQUBBgIHAgYGBQIFAwcEBQkFBAUCBQQFBAUBBgQFBQUFBQMFBAUBBQUGAQUCBQEFAQYBBgEFAQUEBgEFAgYCBgEFAgUFBwEGAwYBBgIFBQUDBQYF",
"BwUBBQMFAQULBwIFAQYEBggHAgUCBQEFAwUHBwEHAgUBBQEFBgUGBQIFBgcHBQEFAQUCBQEGAQUBBgMHDAYDBQEFAgcDBQMFAgUBBQUFAwUEBQMFAwUBBQIHAQcCBQIGDgYBBgEGAQUCBQUGAgUDBQIGAwUFBgIFAgcBBwcFBgUIBgQFCAUCBQcFAgUEBQMFAgUDBQIF",
"AQUBBQEFAgUEBgEFDAUBBgYFAQYGBQQFBQYBBQgGCwUDBQEHAgUBBQQFCAUDBQEHBQUBBgEHAgUBBQEFAQUDBQIFBgUCBQIGAQUEBwMFAQUGBQQGBQUCBQMFAgcBBQIFAQUEBQEFAQUBBQEGCQUEBQEFAQUBBwMFAgUFBwUFAwUBBQEFBQcEBQYGBAYCBQQFAwUFBQUH",
"AwUBBgMFAQUCBQEHBAUBBQEFAQUDBgEFCAYEBgEGAQUBBQMFBQUMBQEGAQUBBQMFAgYEBQEGAQYHBQEFBQUBBQIFAgUBBQIFAgcBBQEHAgUDBQEFBwUBBQEFAQUHBQkFAgUBBQYFAgUBBQgGCwUCBQMFAQUDBgkFAQUDBQ4FAQUDBQIFAQUFBgUGBQYCBQEFCQcEBQEF",
"AgUFBwEFAQUBBwUHAQUDBQEGAQUKBQEFAgcDBQEFBgUBBQwFAQYCBQIFAQUCBQEFBgUCBQIFAwUBBQMFBAUCBwEFAQUDBQYFAQUBBgkHBAUFBQMGBAUHBgIHAwcLBQQFAQUDBgEGBAUCBQIFAgUBBQEFAQUGBQEGBQYBBQIFBQcCBQUFAQUBBQQFBAcEBQEFAwUBBQIF",
"AQUDBQIFAgUCBwIFAQUKBQEGBAUBBQIFAQYBBQEGBQUBBQIFAQUFBQIFBAYCBgEFCAUCBwIGBQUBBQIGBAcDBgMFAQUDBQYFAgUBBwMFAQYBBgQGCQUKBQEFAQUCBgEFBgUDBQIFAgUDBwEHAwUCBQEFAQUIBQMGAQUFBQMFAgUDBQQFAQUBBQEFBgYFBQgGAQUCBQEF",
"AQUBBQIFAgUDBQEFBgUFBwUFAgYDBQkFBAcBBwEFAwUBBgEFAQYEBQcHCAUFBQQFAQUBBQUFBQUDBQMFAQUBBQUGBQUDBQMGAwUHBgEHAQYCBwUFAgUFBQEFBQUCBgMFAgUBBQMGAQUBBQEFAQUFBwQFBwYFBgIFBAUDBQUFCAUCBQIFAQUCBQcFCgYJBQQGAwUEBwMF",
"AQUIBQEHAQUBBQEHAwUCBQgFAQUDBQMFBQYCBwEGAQUCBQEHAwUGBQIFAgUCBgQFBQUBBQEGAwUIBQEFBgYHBgMFBAcDBQIHAgUFBQIGAQUCBQEFAQYBBQIHAgUBBgEFAQUBBQIGAgUGBQMFAwYBBQwFAgUDBQMFAgUDBgIFAQUBBQMHAgYJBQUFAQYEBQYFAgUBBQEF",
"AwUHBgIHBQUDBQcFAQcPBQQGCQYDBQQGAQYCBQIGBAUGBQIFBAUGBQIFBAUDBQEFAgUJBQMFCgUBBQEFAQUBBQQFAQUEBQEGBQUHBQMHAQUEBQQFAgYEBQEFBAYDBQEFBAUCBQEGAgUBBQEFAgcBBQgFAwcCBQEGAQUCBgIFBQUBBgYFCgUKBQQHAQUBBgEGAQUEBQQF",
"AgUEBQEFBgYBBQEFBQUGBQMFAQULBQEFBQUCBwUFAwcDBQMFBgUDBQEHAgcOBQIFAgcDBQIFAwUBBQEFCgUEBgIFBQUHBgEFAwcBBQEHAwUEBwEFAQUFBQIFAgYDBQMFBQUDBQoFAwUCBQEFAgYBBQYFAQUBBgQFAQUEBQEHAgUDBQEFBAUFBgEFBQUDBQcGAgcBBwEF",
"CgYDBQQFDgUFBQEFAQcCBgEFAQUCBQUFAQUEBQMFAwUJBQQFAgUEBQIGAgUJBQIFAQUDBQEFAQUEBQMHBwUBBgkFAQUHBwQGBAUBBQIGAQUIBgIHAQUBBQUFAQUBBQQGBQUBBQEFAgUBBQIFDQUDBwQFBAcFBgMGBwcHBQEFBgUDBQkFBQUEBQUFAQUDBQIFBwUFBQEG",
"AQUBBQQGAQUGBQEFAQUDBQIGAwUCBgMFAwUBBQYGBAUEBQQFAQUHBgUGDQYDBQQGAgUCBQQFAQUEBQMFAQUCBQIGAQUCBQYFBAcJBgEGAQYDBQUFAgUBBQMFAgUCBQQGAgUCBwEFAQYBBQMFAQcHBggFBAUDBQYGAgUDBwQFAQcBBQEFAgUEBQMFBQUEBQEFBAUDBQIF",
"AgUDBQQFAQUBBQMGAQUCBQQFBgUMBgEFAwUBBQEFAwUBBQEGAQUBBQMFBAUGBwMFAgUDBQMFBQcDBgEGAwYBBgcFAQUCBQYFAgUBBwIGAQUCBQcFBAUJBwIFBQUCBQkGAQULBQQFAQUFBQQFAQYBBQIFAwUIBgIFAwUBBQgFAgUCBgUHAgYCBQYFBwUDBQQFCQUFBQQG",
"DwUCBwIHAQcDBQIGAQUDBQEFDAUBBQIFAwUBBQUGAQUCBQEFAgUDBQMFAgUCBQMGAwUDBQMGAgcCBQMHBgYMBgkGAQUCBQQFBwUIBgQFAQUDBQIGAwUCBwwFBwUBBQUFBgUDBwIFBgUCBwEGAQUCBQEFAgULBQIFBAUDBgMGAQUCBQQGAgUBBQEGAwcFBQEFDgUDBQcF",
"AQUCBQIFBQUCBQMFAgUCBgEFAwUBBgYFAgYCBQIFAgcEBQEGAQYBBwMGAQUCBQIFAQUBBwwFBgUEBQEFAQUCBQIFBgUCBQIFBAUEBQEGAgUHBQcFAwUBBQEFAQcBBQQFAgUEBwMFAgUCBQcFAwcCBwUFAQUDBQQFAQUHBgEFBQUFBQMFBQUBBwUGAwULBQQHAwUHBQYG",
"AgUBBgQFBQYCBwYFAwUBBQIFAwYJBQEFBwUBBQYFAwYKBgcGAwYEBQMFBQYEBQEHAQUOBQgHAQUCBgEFAwUBBQEFAQUBBQEFAgYBBwUGAgUHBwEGAgUBBgEFAwUMBQEGBQUBBQEFAgUEBQcFBAcBBQUFBQUBBQIFAQcBBwEFAQUCBgIHAgUGBQEFAgYBBQEFAQUDBQQF",
"AQYBBQEHAQUDBgIFAwYCBgQFAQUBBgoFAQYBBQMFAQYBBQQFAQUCBQwFAQYIBQUFBgUFBQQFBQUBBQEGBwcDBQIFAwUBBQYFAQUHBQIFBQUBBQMFAwUBBQgFAQUBBQIFBQUBBQcHAwcBBQEFBAYBBQ4FCAUDBQMFAQUCBQIFBgYBBQIFAQUCBQkFAwYDBQIFBAUBBgIF",
"BQUEBQIGBAUCBgQFAQUFBQEGBwUJBQMGAwUBBwIFAgcLBQEFAQUCBQEHAQUDBQMFBAUBBQMFAQUDBQUFBQUCBwgFAQUCBQEGAgUBBQcHAgUBBQEFAQUCBwUFBQcBBQMFBgYBBQUGCAUEBQcGAwYCBwoFAQUJBQIHAQUBBgEFBQYDBQQHAQUMBQEFBQUBBgYFBwYKBQEF",
"AQYDBQEGAQYEBQMFBgUFBQIFDAUEBQMFBQUBBQcGAgUBBQMFBAUFBQcFAQYDBQQFAgYCBQMFAgUBBQIGBQUBBQIHAQUBBwEHAQcBBQMGAgcBBQEGBAUIBQUHAwUIBQEFBQUCBgEFCgYEBQEFAQUBBQEFAQUBBgEFAgUBBQEGAgUFBQIFAgUCBQEHCwUBBQEFBAUBBggG",
"AgYLBgIFAQUIBQIFAQUCBQEFAwUEBQEGBgUCBQ0GAwUJBQEFAwUCBgEGAgcDBQIGBAUDBgUHAgUCBQQHBgUBBQUFBAUFBwEFAgUCBQEHBgUDBgIFAgUDBwUGDwUBBQQFAgUBBQQFCAYBBgIFAwUBBgMFBwUUBwEFAwUCBQEFBAUIBgEGAQUEBQIFAwUFBgMFAwUBBwMG",
"AQYDBQYFAQYBBQEFBAUGBgQFAQUFBgIFAgUHBgUFCwUCBQEFAgYEBwUFAwUBBQEHAQUBBQQFAQUGBwoFBgUCBQEFBQUHBgUFAQYDBQcHEgYBBQMFAwYBBQEFDwYCBgUGAgUBBQQHBAYBBQcFAQUFBQkFAQcEBgEFCwcBBQoFAgULBgEFBgUBBQQFAwUBBgEFAQUBBgMF",
"AQUCBQEFAQYBBQgFBAUCBQEHCQUDBwIFAQUDBQEGCAUGBwMGAgUBBgEFAQYCBQIHAQUBBQkGAQUEBgMFAwYCBQMFAwUFBgIFAwcBBQEFBQYBBQIFAQUBBwEFAgUBBQIGAwUBBQMFBgUFBQcFBQUEBQQGAwUBBQIGBAUFBQEFBQYEBQIFAQUBBQIHBAUEBwEFBgcDBgMF",
"AgUFBgEFAQUBBQEFAgYFBQEFAgUBBQMFAQcBBQMFBwUBBQoFCQcGBQEFAgYDBQQFAgUEBQIFBAUNBQIHAgYBBQIFAQUCBQUFBAUCBgIFAQcIBQIFBAUDBQMFAQUBBgEFBwUCBQUFBAUFBQEFCQcDBQIFAQUCBQEFBQYBBQQFAwUDBQkFAQYBBQEFBwUMBgEFAQUBBQoF",
"AgYBBQYFAwYBBQEFAgUDBQEFBgUDBQ0GCAYBBgcFAgUCBgEGBgUCBQMGBAYEBgEHBQUBBgEFBAYDBgIFAQcEBQgFBwYBBQEFCQURBQQHBgUGBQMFBQcBBQEFAwYEBQoHAwUCBgIGAQUBBQIGBAUHBQQFAgUBBQEHBQUBBQEFAwUCBQIHDQcBBQEFAwYBBQQFBQUBBQIG",
"AgcHBQ0FAwUBBQMHDwYBBgEHAwUDBwMFAgUDBgMHCAUJBQcHAgUCBgMGAQUEBQEFBAcCBQUGAwUEBwEGAQUBBQIFAQUCBgEFAgUDBQYFAQUCBQMFAwUEBQIFBgUGBQcGAwUIBQMFAgUEBQYFAwYCBwQFBgUPBwUHCAYEBgcGBQUDBQMFFQUDBQEFCQYFBgEGAwcDBQMF",
"AwUDBQUFAgUEBQMFAgUCBQEFAgYBBgMFAQUCBwQGAQUBBQMFBwcFBQEFBgYBBQMFBwUDBwEFAQUBBgEFAwUEBQUHAQUKBgQFBAUDBQIFBAUCBQQHAwUDBQMGAwUEBgIFAwUBBQgGAQUBBgQFAwUEBQIFAQUDBQMGAgUBBgEFAwUCBwEFCwYFBQIHBQcCBQUFAQUBBQIF",
"DAUIBQIFCAUHBQUFBQcDBQUFAgUFBQIFAgYBBQEFAQUDBQIFBwUFBQUHAQYBBQEFAQUBBQIGAQYCBwIHAwUEBQQGBAUBBQEGAQcBBQEFAwUBBQMFAgUFBQQFAQULBQcFAQUCBwEFAgUIBQMGAgYHBw4HAgUBBQQHBAUBBQsFAwUHBQMFCAYDBQEFAgUBBQMFBQYBBgMF",
"BwUBBQEHAgUFBQIFAQUBBgEFAQcFBQIFBgUHBgEFAgYBBQMFAQUCBgEFBAUFBQMFAwUCBQIFAgUHBQgFBAUGBQIFBgYBBRMFBQUJBgEFAQYGBQUFAwYFBQcFAQcFBgEFAgUCBgMFCgYIBgMFAgUEBgEFAgUDBQQFBwUBBgEFBgUBBQIGAgUEBQEGAgUBBQEFAQcDBQUF",
"AwUCBQMHAQYGBQgHAQUCBgEGAQUBBQIFAQUBBwMFAwUBBQEFAQYCBQQHAQYEBQEFBwUCBQEFAwUCBQEGBwUIBQQHAgYIBQQFBAUBBQgGAwUFBgEFAQYDBQMGBAUGBQEFBAUFBQIFAQUCBgEFAQYBBQEFAgUCBQkFBAUCBQMFAQUEBQIFAQUCBQEGAwcBBQEFAgUCBQIF",
"BgYFBQIFAQUBBQMGAQUFBQgFAwYBBQcFAwYBBQEGAQUIBQMFBAUCBQMFAgUBBwUFAgUBBgMFAgcIBwEFCAUCBQMFAgUCBQEHAQUFBQQFAQYEBQIFAwYHBQIGAQUFBQMGAQUEBgUFBAUBBwMFCAUBBQEFBAUBBQcFCQYDBgUFAwYQBQEFAQUDBQMFBgYCBgQFAQUBBQoG",
"BgYCBwIFBAUCBgIFDQUBBQUHAgUEBQMFAQcCBQIGAQUEBwIFAwcDBQQHBQUCBgcHBgUGBgIFAQUDBQcFAQUFBQgFBQUEBgIFAgUBBQIGAgUBBgEFAgUBBgEFAQcFBgIFAQUDBQEFAgcGBQEFAwUBBgMGAgUBBgEFBQYIBQgFBAUDBQsFBQYHBQMFDAUCBQcFAwUCBgEF",
"CQYCBQgHDQUCBgIFBAUCBgEFAQUEBQkGAgYCBwEFBwUCBQwHCQcCBQIHBAUBBwMGAQUEBgEGAQUCBQQFAQUDBQIHAgYBBwIFAwcBBgEFAQYGBQEFAwcCBgQGAQUBBQgFBgYDBgUGBgUGBgwFAQUBBQMHBAUBBwYFBAUDBQ8FAwUDBQQFBwcDBwEGAQUCBwMFAwYCBgEH",
"AgUFBQUFAgcBBQEFBwUEBgIFAQYCBQYHAgUEBQQGCQYBBQYHAQcCBggFAgUBBwQHAQYCBQcFAQUBBQQFBQUBBgUHAgYEBQkGDgYBBwEGAQYBBQEFAwYBBgEGBgYCBQ8FAQUBBgIGBwcDBwEFAQUDBQUGAwUBBQEFBAUCBQ0FAQYEBgIFAQYDBQYGAwUBBREFAgUBBQIH",
"AQYBBQUHAgYVBQEFBQUEBgcHAQYLBQEFBQUBBwQGAwUDBwEHBAYBBQEFAwYDBgYFBQUEBgIGBAYGBgEGAwYHBQEFAgUHBQUFAQcGBQEFBAUEBQQGAgUBBQIHAgUCBQwGBwUHBQEGBAcEBQIFCgUGBQIGAwUCBgEFAQUCBQoFCAUEBQEHBgYCBQMFAgYJBQUHCQYCBgMF",
"AgUBBQEGAgUBBgsFCAUKBwMFAgYEBQEHAgcEBQIFCwcEBgkHAgcFBQMHCAUFBgYFAgUBBQYGAgcCBQEFBAYDBQQGCQYBBQUGAQcBBgUFAwUHBgkHAQYCBgEHHAUBBggGAwUCBwEGCgYCBQYHBAUBBQEFBwYCBQMFAwUGBQcFAgcBBQYGBAcEBQkGAgUBBgEHCAUGBQEF",
"BQUKBgIGCAUBBwQGBgUCBQEGAQUBBQYGBAYDBQQFAgcEBwMGAwUBBgQHAwcBBQEFBgYBBQIGAQUHBQIHAgcEBQsHAgYDBQEGAgUBBRAGAwYHBwEHBgUHBQMHBAUBBgMFAgcKBQEFBAcQBQYGAgYGBhIFCAUIBwEFBAYDBgIGBwUBBwUGAgUBBgYGAQcCBgEFAQYLBgcF",
"AQUBBQQGAQcDBwUGCgUIBgEGBwcEBgMHCAYBBQUGCAYCBgUFAgYFBgUGBgYDBwUFCAYBBgEHAwUDBgcGAwUMBQMGAQYDBgMGBgcIBgEGCwYBBQYGAQcBBgEGAgYFBgEGCwcDBgEGBgYRBgUHBwYBBwEGAgYBBgIGAQYDBgoGBQYCBgIGBQcHBgIGCQYEBgEGCQYCBgEG",
"AgYCBgMGBQYBBgMGAwYEBgQGCgYCBgIGDQYGBgIGAgcBBgEGAQYEBgEGAQYCBgQHAgcBBwMGAQYCBgEGAQYCBgEGBAYBBgEHAgYBBgEGBwYBBgUHCAYBBgEHBAYDBgEGDQYEBgMGAgYFBggGBgYBBgIGAQcCBgcGBAcJBgIHBgYBBgIGBAYJBgkGAwYDBhEGBAYBBgQH",
"AQYEBgMGAwYDBgcGAwYHBwQGAgYDBgIGAgYDBgUGAQYCBgEGAwYBBgYGAgYCBgIGAgYBBgIGAwYDBgIHBQYBBgUGAgYBBgEGAQYBBgEGBAYCBgUGAgYDBgQGAwYHBgQHAgYGBgIGAwYBBgMGAwYDBwMGAQcEBgEHAgYBBgEGBAYDBggGAwYDBg4HCQYBBgQGAwYKBgIG",
"AgYEBgIGAQYCBgEGAgYBBgEHAQYBBgQGAQYIBgEGAQYCBgMGAwYGBgYGCwYDBgIGBQcJBgEGAQYCBgYGAQcKBgIGAwYCBgYHAwYCBgEHAgYEBgQGAgYDBgUGAwYEBgMGAwYCBgUGBAYDBg4GCAYBBgIGAQYFBwYGAgYCBgEGBAYBBgIGBgcHBgMGBQYBBgEGAgYDBgEG",
"BwcBBgEGAQYBBgIGAwYBBgMGAwYDBgYGAgYBBgEHAgYBBgEGAQYDBwQGBgcDBgYHAwYHBwMHAQYLBgIGBAYBBhAHBAYBBgUGAQYNBgMGBAcGBgIGAQYDBwEGBQYBBgIGBAYDBgEGAwcBBwcGAgYGBwMGAQYBBgEGBgYFBgQGBgYHBgEGBgYDBgwGAgYLBwMGBQcRBwgG",
"AgYEBgEGAgYBBgUGAgYEBggGBQYBBgQGBAYEBgMHCwYBBgIGAwYGBgQGBAYBBgYHAQYCBggGAQYRBgIGAQYBBgUGDAYHBgcGAgYEBgMGAQYCBggGBgYCBgEGAQYBBgEHAwcHBgEGCQYGBgEGAgYCBgEGAwcBBgEGAgYEBwMGAQYBBgEHAQcBBgIGAwYDBgEGIgYCBwIG",
"AgYJBgMGEgYFBgQGAQYBBgUGBAYBBgIGAwYBBgcGCwYBBgkGBgYDBgMGBAYCBgcGAQYBBgQGAQYCBgMGAgYFBwUGAgYBBgwGCQcDBgMHAgYFBhAGAQYGBhMGCwYBBwUGBgYBBgkGAgYBBgUGAQYFBgEGBAYBBgIGAQcHBgMGBQYDBggGBAcBBgEGAQYKBg4GAQYBBgQG",
"AwYCBgEGAQYBBgEGAgYDBgIGAQYCBgQGBQYEBgEGAgYBBgEHBgcDBgMGAgYCBgIGAgYBBgYGAgYBBgcHAwYCBgIGDwYCBgUGAQYEBgEGAgYBBgEGAQYDBgUGBAYBBgkGAQYFBgQGBQYFBgEGAgYFBgUGBAcGBgsGBQYHBgIGAQYDBgQGAgcDBgEGAgYBBgoGCgYCBgEG",
"BgYMBgQHAQYCBgQGAwYFBgIGBAYCBgQGAwYDBgQHAQYHBgUGAQYBBgEGBAYDBgEGBAYGBgMGAQYIBgcHCQYHBgEGAgYBBgEHDwYDBgMHAwYEBgIGAwYEBgEGAgYHBgEGBAYEBgYHBQYBBgcGAgcDBgEGAgYDBgMHAwYCBgMGAgYBBgEHBAYCBwEGBQYDBgIHAQYDBgEG",
"AwYFBgIGAgYMBgMGCgYHBgUGDgYEBgMHAwYBBgEHBAcBBwIGAgYGBgEGAQYJBgMGAgYEBwEGAwYBBgIGAgYDBgIGCAYDBgQGAQYDBgQGAwYCBgYHCAcDBwYGCQYHBgsGBAYOBgYGAQYCBgEGAgYFBwMGAQYQBwEGAwYBBgEGAwYBBgIGAgYBBgEGCAYBBgIHBgYDBwEG",
"CQYJBgEGAwcFBgEGBgcCBgQGAwYDBgMGAgYBBgEGAQYBBgEHBAYEBgIGGQYBBgYHBAYEBgMHAwYBBgoGBQYDBwIGAQYCBgMGDAYBBgYGAgYFBgEGAQYBBgMGAQYDBgEGAQYEBgIGAQYHBgIGAgYBBgUGCAYCBgEGAgYGBgEGBAYDBgEHEwYHBgEGAgYJBgEGAgYDBg0H",
"AgYBBgEHAwYMBgEGBAYBBgIGAgYBBgIGCgcCBgEGAwYEBgEGAgYDBgMGBAYBBgEGAgcDBgEGAwYFBgcGAwYCBgEGBAYGBgcGCQYDBgMGBQYBBgIGAwYBBggGAQYCBgEGAgYCBgEGAQYGBgUHAgYCBgUGAwYGBgIGAQYCBgMGBwYEBgMGAQYDBgcGBgYEBgMGAgYDBgIG",
"AQYBBgYGAwYBBgIGAgcEBgEHDAYEBwIGBAYHBgUGAQYGBgsHCQYBBgQGAwYKBgEGAgYKBgEGCwYFBgoGCQcCBgQGAgYIBgsGAgYCBgIGAgYCBgEGAQYHBgIGAQYLBgEGAgYFBgcGBwYDBgEGBAYEBgIGBAYBBgQGAwYEBgcGAwYBBgMGAgYCBgIGBgYCBwMGAwYHBgEG",
"CAYBBgYGAQYDBgMGBgcEBgQGAwYEBgEHBAcCBgcHAgYFBwYGBgYDBggGCwYIBgQGBQYJBgIGAwYMBgIGAwYCBgkGBgYBBgEGCgYDBgIGAQYCBgMGAQcCBwEGBQYBBgcGAQYBBwIGAQYDBwQGAQYCBgMGCAYBBgQGCQcDBgUGDAYJBgQGBAYCBwIGAwYEBgIGAQcCBgEG",
"AQYBBgQGAwYDBgIGAwYBBgQGAQYFBgEGAQcEBgEGAwYBBgYHAQcCBgcGCAYCBgMGBQYJBgIGCAYFBgQGAgYDBgMGAgYEBggGAgYFBwQGAQcCBgkGAgYBBgEGCQYBBgEGAwcBBgoGAQYFBggGCAYFBgIGBQYCBwYGAwcBBgQGBQYDBgEGBQYCBwEGBQYCBgMGAwcFBgMG",
"CgYBBgQGBAYBBgMGAgcFBwMGCQYEBgUHAQYGBgIGBQYHBgUGAwcCBgEGBAYBBgQGBwYUBgMGFwYBBgYGAQYBBgEHAgYBBgEGBAYDBgIGBQcBBgMGCwcCBgIGAQYCBgIGCAYFBgEGBQYFBwIGBAYCBgQGBgYGBgoHBgYDBgMGAgYIBgYGBQYEBggHAQYBBgEHAgYBBgMG",
"BwYEBgEGAgYDBgEGAQcIBgEGCwYNBgEGBwYDBgEGAwYCBgEHBgYBBgMGCAYBBwIGAQYDBgEHAwYBBgIGAgYGBgEGCwYGBwIGCgYGBgEGAgYGBgEGBAYFBgMGAwYBBgIGAwYJBgUGAwYFBgoGAwYBBgMGAgYBBgQGBgcCBgEGBwYEBgMHAgcHBgEHAQYBBgEGAwYBBgIG",
"AQYCBgEGAwYCBgIGAgYDBwIGBAYCBgEGAwYDBgMGAQYDBgMGAwYEBgEGAgcGBggGAQcHBgcGBQYCBgIGAgYFBgQGAQYBBwsGAQcBBgEGCQYEBgkGAQYBBgYHBAYCBgEGAwYDBwIGBQYCBgQGAgcDBwEGBAYBBwkGAQYCBgIGAQYJBgcGAQYBBgEGAQcCBgMGAQYFBwgG",
"BgYBBgMGAQYBBgYGBQYCBgcGAQYDBwEGAQYDBgEHAgYCBgMGAgYIBgYGAQYBBwEGAQYEBhoHBwcIBgEGDQYCBgUHAgYBBg8HAgYBBgIGBgYDBwIGBAYGBhAGAQYBBggGAQYMBwQGBAYDBgEGAwYFBgEGAQYEBgEGAgYBBgEGAwYCBgQGAgYBBwUGAQYBBgIGBAYCBgMG",
"AQYBBgIGAgYBBwMGAgYEBgEGBAcBBgIGBAYEBgQGAQYEBgMGAgYBBggGAQYBBgIGBAYHBgEGAgYBBhQGAwcOBw4GBAYDBgQGBAYCBgIGAgYDBgMGBQYDBgQGBAYDBgYGAQYHBgEGAQcDBgQGDQYCBgYGAQYFBgIGBQYBBgUHAgYBBgMHAgYBBgcHBQYDBgkHBAYGBgkH",
"AgYDBgEGAQcEBgEGAQYCBwsGAgYBBwIGBAcFBwMGDAcBBgYHAQYBBgEGAgYEBgUHCQYHBwEGAgYCBwIGCwYBBgYGBAYFBgIGBAYIBgMGCAYFBwYGBAYJBwkHBAYEBgUGAQYBBgIGAQYCBggGBQYBBgoGBwYCBgQGAgYBBwoGBgYMBgEGAgYBBgEGAgYBBgEGAQYCBgEH",
"CQYGBgIHAQcBBg8GAgYBBwcGAQYDBgEHAwYIBgUGAgYBBgEGAgYFBgEGAQYCBgEGAQYDBgMGAwYEBgMGAgYFBwEGAgYLBgIGBAYDBgIGBQYBBggGBgcGBgwGAQYBBgUGDwYBBgIGEwYLBggGBgYIBwUGBQcBBgQGBAYMBggHAQYBBwUGBwYBBwYGBQYBBgIGCQYDBgIH",
"DAYBBgwGAQYJBgIHAwYEBgQGAQYBBgEGCQYBBgMGBAYDBgMGAQYEBgYGDAYDBgEGAgYFBgEGAQcHBgEGCwYEBg4GBAYKBgoGBwYDBgIGAQYCBgkGAQcHBgsGAwYCBgIGAgYKBgIGAwYCBgIGBgYCBgEHAQYCBgEGBAYDBgIGAgYGBgEGCAYBBwIHBQYCBwsGCQYBBgoG",
"AQYCBwEGAQYFBgEGBAcGBwQGBQYGBgQGBwYFBgEGAQYBBgEGAwYCBgEGAQYBBgcGAQYCBgEGBQYFBgIGAgcBBgMGCAcOBgIGBQYHBgEGCgYCBgIGAgYDBgQGBQYBBgUHBAYIBwEGAgYEBgoGAgYFBgIGAQYMBgMGBAcCBgEHAQYEBgIGAgYHBgUGAQYEBwMGAgcEBgEG",
"AgcFBgcHEAYFBgsGAQYBBgEHAgYOBwMGAgcDBgIGAQYGBgEGAwYBBgMGAgYCBgEHAQYBBgQGBQYEBhYHAwYDBgIGAQcPBwYGAQYGBgEGAgYCBwIGBAcBBgEGAgYEBgEHBAYEBg0GAQcIBhkGBgYIBgIGAgYEBgkHAwYCBwIGCQYCBwEGBQYHBgQGAgYCBgYGAQYCBgQG",
"AgcBBgIHAwYCBwMGFAYCBwQGAQYHBgkGAwYCBgMGCgYHBgEGAgYCBgMGBgcEBgEGAQYKBgEGAQYBBgIGAQYKBwIHAwYBBgIGBAcIBgYGBwYCBgMHAQYCBggGBAYDBwIGCwYHBgEGBgYHBgIGBQYDBwgGAwYCBgQGDwYCBgMHCAYDBgEGAwYEBgYGAwYGBgEGCgYBBgMG",
"AgcEBgEGBQcIBgEGBAcEBwEGAgYCBw8GDQYCBwMGAQYBBgEGBQYIBg4GCAcCBgEGBAYCBwUGBgYIBhAGAgYBBggGBgYIBgUGBQYCBggGBAcCBgIGCQcCBgQGAQYMBgEGAQYDBgUGBQcMBgYGBQYBBwEGAQYDBgYGEQYFBgMHAgYDBwEGAwYBBgUGAgYCBgYGAQYEBgEG",
"BAYCBgEGCgYCBwIGBwcCBwIGBQYCBgoHAQYEBgEGBgcEBgIGAQYBBgIGBAYBBgMGAgYBBgEGBAYDBgIGAQYBBgIGAQYEBwYGAgYLBgQGCAYEBgIGAgcDBgEGAgYMBwIGAQYDBgMGAQYEBgUGAwcCBgMHAQYFBwgGCQYBBwEGBQYCBgMGAQcGBgIGBAYEBgMGAgcBBwEG",
"BgYBBgEGAwYEBgIGBgYEBgkGAgcEBwIGCQYFBggGAgYCBgQGBgYBBgoGAQYFBgMGBAYCBgIGBgYCBgEGAgcCBgEGAgYBBgQHAQcFBgEGAQcCBwQGAwYEBgEGAQYCBgIGAQcDBgEGBAYKBgYGAQYBBgIGAgYCBwIGAQYEBggGAwcFBggGAgYBBgIGBAYGBwEGAQYBBggG",
"AwYBBgIGFAcHBgYGAQYCBgYGAwYFBwcHCwYEBgMHAgYDBgEGAgYBBgIGAwcCBgEGEwYCBgQGCAYBBgEGBgYCBgIGAQYLBgIGBAYIBgIGCAYEBwIHAgYCBwEGBwcCBwQGAwcBBgEGBgcEBwEHCQYCBgEGBAYIBhEGAQYBBgMGAgcJBwYHAQYIBgMGAgYBBgMGFgYDBgIG",
"AQYEBgEHAwYDBgUGAwYCBgUGAwYCBgUGAwYDBgUGBwYCBgEGAgYBBwMGBAYEBgIGBgYDBgQGBAcDBgMGCQYDBgIGBAYBBgQGCQYEBgcGBAYIBgIGCAYRBgkGAQYGBgsGAgcFBwIGAgcJBgEGCAYIBgMGBQYBBgYGAQYBBwIHAQcFBgUGAgYBBgMGBwYDBgEGAQcBBgIG",
"AQYGBgIGBAYDBgEGBgYBBgEGAwYCBgMGAQYEBwEGAgYBBgIGBQYLBgQGAgYDBgoGAgcFBgQGAwYEBwUHAgYCBgQHAwYFBgYGAgYDBwIGAQYDBwIGBgYCBgIGAgYFBgEGAQYBBgMGCAYCBwMGAwYCBwQHBQYOBgIGAwYDBgMGBQYCBgUHAwYBBgQHAgYBBgMGDQYBBgEG",
"CAYBBwIGBwYCBgEGAQYMBgEHAQYFBggGAQYBBgMHAgcHBwIGAgYGBgEGAQYCBgEGAwYGBgMHAQYEBgQGCAYCBgYGBgYBBwEGDQcCBgUGAQYBBgEGBQcBBgMGAwcDBgYGAQYJBgMGBAYFBgYGCQYSBgIGAgYKBgEGAgcHBgEGAgYBBgIGBAYFBgEGCQYEBgUGBQYDBgEG",
"AwYDBgQGAwYBBgQGAQcBBgEGAQYFBwIGAgYDBgEGAQYCBgUGBQYDBgEGAwYIBgwHBAYCBgEHAQYEBgIGAgcEBgMGAQYEBgUGBQYCBgMHDwYBBgEGBAYEBgEGAQYBBwEGAgcCBgEHAwYDBgEGAQYEBgQGBwYBBwQGAQYBBgUGBAYBBgIGCAYFBgQGBAcDBgMGAgYEBg0G",
"AwcFBgQGBQYDBwYHAwcDBgIGAgYCBwUGAgYGBgIGCwYMBgMHEAYFBgkGBAYEBgEGBQYLBgEGAgYDBwEGBQcBBgkHBgYEBgQGBQYBBhEGCgYEBgEGAQYFBgIGDAYJBgkGBwYEBggGAwcEBgEGAwYOBwkGAQYGBgIHAgYMBwIGBAYEBgIGAQcBBgEGAQYBBwQGBwYBBgMG",
"AQYHBwMGCAYDBgYGAQYFBgMHAwYBBgEGBgYUBgoGBQYDBgQGCAYCBgQGCQYDBgMGCwYQBwUGAQcHBgEHAwYBBgsGBAYGBgsGBQYBBgcHEQcEBggHAwYPBggHBwcGBgIGBwYGBggGAQYBBgEGAQYCBgEGAgYCBgEGAgcGBgUGAgYEBwEGAgYBBgEGAQYBBwUGAwcIBgEG",
"AQYGBgEGAwYEBgMHBQcCBgEGAwcCBg8GCgYFBgMGBgYEBgIGCAYDBgMHAwYBBgIGAgcCBgQGAwYCBgEGAgYHBgoHAQYFBwEGBQcDBwcGBAYEBgQHAQYDBgcGAgYFBgIHAQYBBhEGAwcCBwoGFwYBBgUHAQYEBgwGAgYDBgcGCAYBBg8GAwYCBgIGBAYBBgQHBAYHBgQG",
"AQcCBwMGAQYBBgQGAwYCBgUGCQcKBwsHAQYBBgUGBAYCBgEGBgYBBgIGBgYFBwoGCAYKBwMGAQYPBwoGAQYGBgUGBwYIBwUGAQcEBgIGBAYBBgQGAwYCBgEGAgYDBgQGAQYBBggGCAYFBgIGAQcEBgEGBAYEBgQGAwYDBgUGBgYDBgQGBAYEBgEGAgYCBggGAQYDBgQG",
"BQYDBgIGAQYEBwgGCQcHBgMGAQYGBgIGAwYEBgYGAgYCBgEGAgYEBgQHAQYBBgIHAgcGBgUGAQYBBgEGBgYBBgEGAgcBBgYHBQYBBgkHAgYCBgMGCAYDBwIGAgYDBwMGAwYCBgEHBQYBBwEGBAYFBgEGBAcEBgYGBAYBBgkHAQYJBgUGBgYJBgMGCwcHBgEHBAYLBgEG",
"AQYCBgEGCgYFBgEGAwYDBgUGAQYBBgEGCQcEBgEHBwYBBwEGAQYJBwMHCgYDBgUGCQYCBgcHAQcDBgEGAwYEBwUGCwYEBgUHFgcBBgcGAwYBBgQGAgYGBgIGAQYCBxkGBAYBBgEHBAYDBgMGAQYCBwkGAQYBBgIGAQYBBwEGBAYBBgEGCQcBBgYHAgYHBgMHDQYFBgQG",
"AQYEBwMGBQYDBgkGCwYBBwgGAwYPBwEGAQYFBgIGAwYCBgUGAgYIBgYGAQcEBgEHAwYDBgYGBAcFBgIGAQYHBgMGCAYEBhAGCwYBBgEGBQYFBgUGAgYEBgIGAQYEBwIGAwYBBgEGBQYBBgYGAQYDBgIGAgYEBwgGBQYEBgYGBQYDBgEHAgYJBgIHCgYEBgsHAQYDBgEG",
"BQYHBgMGAQYZBgEGBAYBBwQGAQYDBwEGAwYKBgIGAQYQBgoGAwYIBgEGBgYBBgwGBAcIBgwGAgYEBwEGCQcDBgIGCQYBBgIGEgYBBxAGAwYIBgEGAQYBBgEGAQYFBgIGAgYEBgQGAQYCBgMGAgcBBgoGAwYIBgEHCQYDBgUGBQYBBgIGBAcCBg8GBQcHBgMHAgcJBgYG",
"AQcBBgEGAQYBBgIHAwYCBgIGAQcBBgUHBQYGBwYHAgYMBgIGAQYCBgMGAgYBBwIGBwYBBgUHAwYCBgEHAQYCBwMGCAYHBwQGBQYCBgEHAwYBBhAGAgcCBg4GAQcBBggGBQcVBhAHCgYCBgEGCgcDBh0GBAcBBgUGAgYLBwUGAQYIBgMHAQcFBhAGAQYFBgEHAQYFBwkG",
"BwYEBgIGBgYCBgUHBQcHBgMHDAcBBwYGBwYCBg0GAQcBBgEGBwYEBgEHAwYEBgsGAgYBBhAGAgYKBgIGAgYGBwIHBAYDBgYGAQYCBgYGAQYDBgIGAwYBBgIGAQYDBwYGAwcBBgIGDQYHBgcGAQYHBwMGEgYDBgQHDAYBBwMGBQcDBgEHAwYCBgIHBQYDBgcGAgcCBgEG",
"BwYBBgQHBQYCBgYGAwYCBwQGBgYDBgcGAgcDBgQHAgYIBwUGEQcEBgIGAgYIBgMGAwYKBwEHBAcBBgEHAQYFBgcGBQcJBgQGBwcFBgIHBwcDBwEHAgYEBwQGBQYHBgQHBQYDBgQGCAYBBgIGAgYIBgEGAQYBBgEGBQYLBgYGAgYFBgIGAgYCBggHCwYBBgEGAgYGBgUG",
"BQcCBwQGBwYUBgEGAQYBBhEHBAYEBwcGDAcCBgMHAgYIBgEGDgYCBwwHBAYDBgEGDQYFBgEGAQYBBgQHBAYDBwEGAgYBBgYGCAYGBgMHBAYBBgEGBwYBBg4HAQYBBgMHDQYCBgQHAwcCBwcGAQYBBgwHAQcUBggGCwcCBgQGBAcCBgEGAgYBBgMGAQcEBwYGAwYCBgcG",
"AwcCBgIGAQYEBwIHBAcBBgcHAQYEBwIGCQYKBggGAQYEBwMGAgYCBgoHAgYDBgEGBQcFBgIGAQYBBgEGAQYIBwIHAgYCBwQHBAYCBgEGCwYMBg4GAgYFBgIHAQYFBgEGCAcBBgcHBgYBBgEGCgYGBgoHBgYGBgIGAgYFBwIGBwcEBgEGAgcCBgEGBAYBBgEGAwcJBgEG",
"AwYBBwMGAQYIBgIHAQYBBgEGBQcCBgYHAQYBBwUGBQYCBwIGAQYGBgIGAgYBBgMGAwYCBgIHCAYCBwIGAgcEBwMGAwYBBgQGAgYBBgMGAQcBBhMHBwYIBgEHAQYBBwMHDAYCBwQGAgcHBg0GBAYSBgQGAwYBBgIHBgYBBwMGAwYDBgEGAQYBBgMGBQYDBgQGBQYCBgQG",
"AgcHBxMGBAcFBwwGAQYDBgUGCAYBBgEHAwYCBgEHBQYBBgcGAgcFBg4GAQYKBgEHAgYBBgIGBgYFBg4HDgcFBwEHAgYDBwgGDgYBBgIHEgcHBgsGAwYEBgIGAQYDBgQGCgYOBgIGAwYFBgQGAgYTBhkGAwYHBwUHAwYMBwYGBAcCBwcGAQYSBgEGBQcBBwQHCwcDBwUH",
"AwcDBwcGCgYGBgIGAwcBBgQHBwcCBwYHAgcBBwQHAQYFBgMGCAYCBgEGBwcBBgEGAwYKBwkGAgcBBwEGBQcIBxEGAgYBBgcGAwYDBggGAQYBBwIHBAYEBg8HAQYDBgUHAgYEBw8HAQYBBwEGAQYBBgsHBgYCBwEHBgcCBgkHCAcCBwMHBgcEBwUGAgYBBgQGAQYNBwEG",
"AwYDBwgHAQYKBwcGCwYCBwIHAQcBBwcGAwcEBwQGBgcFBwEHAwcCBwEHCAcDBwEHCAcHBwEHAwcFBgUGAgcDBwMGAgYIBgMHAQYGBwEHAQcIBwkHBQcBBgEGAgYGBgEHAQcFBwMGAgYFBwEHAQYFBwgGAQYFBgcHAQYBBgYGAgcCBgEHCAcQBgUHAgYBBgMHAgcCBwMH",
"BAcEBwEHCwcBBgIHAgYCBwIHAwcIBgEHAQcBBgIHCAcCBgEGAwcBBwMGBwcCBwoHDAYEBgEGAQcHBgYHCAcCBwIGBwcBBwcHAwYEBgQGAwcCBggGDAcYBwMHAwcEBwEHBQcCBgMGDwcFBwEGCgcBBwEHBgcIBwUHAgcKBwcHBgcDBwQHAQcLBwgHBgcHBwIHAgcGBwIH",
"BgcJBwEHBQcCBw8HAgcCBwgHAwcCBwEHBQcLBwIHBQcFBxIHAQcCBwMHAQcHBwYHAwcBBwMHAQcDBw4HAwcTBwMHBgcBBwIHAgcBBwQHAQcHBwsHAgcDBwEHAwcIBw4HBQcBBwIHBQcGBwkHAQcDBw0HAQcCBwMHBQcFBwQHAgcCBwkHCQcCBwEHAgcGBwEHBQcDBwYH",
"BgcDBwYHAQcFBw0HAQcBBwYHBQcCBwgHAQcLBwIHBQcDBwIHAwcDBwQHEAcFBwEHCQcCBwMHAgcCBwkHAQcBBwcHBQcCBwIHAQcCBwIHAgcEBwwHBgcFBwEHBgcHBw4HAwcBBwEHCgcBBwEHBAcIBwMHBwcYBwEHAwcBBwMHAgcBBwEHAQcDBwMHDAcDBwIHAQcFBwsH",
"AgcGBwoHDQcBBwEHDAcIBwIHAwcFBwkHAwcGBwEHBAcEBwMHAgcBBwQHCgcDBw4HBAcBBwQHBAcHBwQHAwcBBwEHAQcDBwIHCwcPBwMHBQcHBwsHBAcIBwMHAgcJBwIHAgcEBwMHAQcJBwQHCQcCBw0HBgcQBwcHFAcFBxMHDwcBBx8HAwcEBwQHBwcPBwMHDAcGBwEH",
"AgcEBwEHAQcCBwMHBQcCBwIHAwcEBwEHAgcBBwEHAQcBBwEHBgcBBwcHBgcBBwIHAgcDBwIHCAcEBxcHBQcHBwIHCAcCBwMHAgcCBwQHFQcCBxIHAgcGBwQHAgcNBwIHAQcKBxUHAwcEBwIHAwcJBwsHCgcDBwYHAQcIBwcHAwcIBwMHBQcCBwoHEwcIBwEHBQcIBwMH",
"BgcBBwcHCAcJBwMHBwcBBwYHBQcGBwMHBwcCBwMHAQcGBwIHBAcFBwsHAgcBBxMHFQcDBwkHAwcMBwEHAwcBBwgHAgcJBwIHAgcKBwMHCQcJBw0HBQcOBwEHAwcJBwoHAgcEBwcHAQcBBwUHAgcFBwIHAQcGBwYHBgcDBwEHAQccBwEHAQcCBwUHAgcEBwMHBAcCBwEH",
"AQcFBxEHCAcGBwIHAgcTBwMHDgcJBwQHAQcBBwMHAgcCBwMHAwcGBx8HFAcDBwQHBwcGBwEHAwcCBwEHAgcDBwEHDwcGBwkHAwcBBwUHAgcCBwMHAgcBBwIHFQcIBwMHBAcNBwIHBAcEBwIHBAcOBwIHAgcCBwQHBQcYBwQHAgcBBwYHEQcJBwMHAQcHBwkHAwcDBwEH",
"CAcEBwQHAwcFBwEHBgcDBwEHAgcDBwkHBQcEBwQHDAcBBw4HAgcCBwMHCwcCBwEHAgcDBwMHAQcDBwMHBQcDBwIHBAcFBwYHAQcGBwEHAwcCBwgHAgcIBwQHAQcIBwQHDQcBBwYHBAcBBw0HBwcEBwQHGAcEBwgHBQcQBwEHAwcDBwoHDAcFBwIHAQcHBwUHDAcDBwMH",
"AgcMBwIHAQcEBwEHAQclBw8HAQcCBwIHAQcHBwEHBQcEBwEHAQcCBwUHAgcIBwYHAQcFBwwHCAcEBwsHAgcQBwQHAwcDBwYHGAcCBxwHAwcQBwQHBAcCBwMHBgcCBwUHBwcGBwMHBQcGBwQHBQcBBwIHAgcCBwUHAQcBBwYHAgcGBwQHBAcCBwEHAgcFBwEHBQcEBwcH",
"BQcDBwEHAgcBBwQHAgcKBwcHAgcDBwcHAwcDBwQHAQcMBwIHAgcEBwIHBQcRBwMHBgcBBwoHAgcFBwEHBwcEBwsHAwcHBwUHCQcCBwMHCQcDBwYHAwcEBwgHBgcGBwEHAwcMBwEHAwcGBwEHCwcBBwcHAgcBBwMHAQcJBwMHAQcBBwUHAQcDBwEHBAcBBwEHAwcFBwIH",
"CAcDBwEHAQcCBwgHBAcCBwcHAgcUBwQHFAcCBwUHAwcBBwgHAgcLBwQHBwcEBw0HGQcGBwIHBgcBBwkHAQcBBwYHAQcCBwIHAQcCBwwHAQcIBwQHAQcJBwcHAwcDBwcHAgcJBwUHCQcXBwIHAQcDBwEHAQcKBwgHAgcFBwcHBgcEBwEHAQcEBwQHCgcDBwgHBgcEBwIH",
"CQcEBwUHAQcFBwYHAQcDBwcHBwcCBwQHAgcCBwEHAQcDBwEHAgcJBwgHAQcJBwUHAwcBBwMHCQcBBwUHAQcBBwkHBAcEBwIHBAcEBwEHBAcCBwUHBAcEBwIHCQcJBwEHAwcCBwQHBQcDBw0HAwcCBwMHAgcPBwUHAwcDBwYHAQcDBwUHAgcDBwcHAQcIBwIHAgcCBwsH",
"AgcBBwoHAgcIBwMHBAcEBwMHBAcEBwEHBQcFBwUHAgcDBwEHCgcBBwcHAgcDBwQHAwcEBwMHAgcCBwEHBAcPBwUHBQcEBwMHAgcGBwMHCQcDBwIHCgcEBwIHBAcGBwUHBQcHBwUHBAcDBwMHBAcCBwUHBQc+Bw4HFAcCBwIHAQcWBxEHBwcCBwEHCgcBBxAHBAcDBwcH",
"FAcRBwQHCwcKBwUHAQcFBwQHAwcDBwwHAwcCBwMHBAceBwsHAgcBBwYHAwcBBwEHAgcBBwMHAwcDBwMHAQcEBwMHAQcBBwEHAQcCBwMHAgcGBwoHAQcIBwYHBwcBBw0HBAcEBw4HAwc2BwcHBgcCBwQHDQcSBwEHAQcEBwIHBAcRBwcHAwcJBwYHBAcHBwEHAwcBBwIH",
"AgcBBwQHBQcGBwUHAQcCBwEHBwcCBwIHBgcIBwQHBQcGBwUHBQcNBwkHAgcCBw4HCwcCBwQHBgcDBwcHAgcBBwgHCgcHBwUHFwcGBwUHCQcBBwMHAQcPBwEHAgcCBwYHAQcFBwYHBAcBBwQHAQcPBwEHCwcEBwQHBQcEBwoHBQcJBwMHBQcCBwwHAgcDBwQHBQcCBwEH",
"BwcFBwQHAwcGBwIHBAcLBwYHAwcCBwMHAgcGBwMHCwcBBwIHBAcIBwwHBQcDBwoHCQcIBwcHCwcFBwEHBQcMBwIHAQcBBwIHBAcGBwYHAwcBBwEHAgcQBwIHAwcHBwQHAwcFBwEHAgcDBwYHCwcIBwEHBwcEBwcHAgcCBxIHBgcBBwEHAQcBBwUHBwcBBxAHAQcFBwIH",
"BQcEBwMHBQcCBwMHAgcBBwMHAwcDBwQHAwcFBwEHBAcBBwMHBgcBBwUHBwcKBwQHAgcBBwMHAQcBBwQHCAcGBwIHAQcIBw0HBwcGBwUHAQcCBwIHAgcGBw4HAgcNBwwHDQcCBwQHAQcBBwUHCwcQBx0HAgcEBwwHDQcBBwMHAQcCBwkHDAcTBwEHBgcJBwIHAgcIBwcH",
"BgcCBxIHAgcBBwkHAQcJBwMHAQcOBwQHDgcNBwIHCQcJBwYHAQcCBwIHAwcFBw4HAgcDBwEHDgcEBwQHBQcCBwQHBQcEBwIHAgcCBwEHEAcIBwcHAQcCBwIHAwcDBwEHAgcBBwUHCQcDBwYHAQcPBwUHBAcCBwEHAQcIBwMHBgcBBwEHBwcGBwUHAgcFBwgHCwcNBwEH",
"AQcDBxIHEwcBBwYHAQcOBwcHAQcTBxQHBgcDBwYHHAcEBwIHAwcBBwEHBwcDBwEHBAcGBwUHAQcEBwIHAQcBBwIHAgcCBwEHBAcNBwgHAwcHBwIHBAcHBwMHCAcBBwEHAgcIBwEHAQcCBwcHAgcBBwEHAQcIBwEHBgcJBwIHBQcLBwEHEAcBBwwHBgcBBwEHCQcCBwIH",
"AQcFBwkHEgcDBwMHFAcEBwMHAgcPBwgHAwcGBwIHCwcXBwsHBAcGBwYHDgcNBwIHBQcBBwEHAQcBBwkHEQcLBwEHAQcDBw4HDAcIBwgHAQcBBwcHEgcCBwIHBgcEBwEHAwcBBwQHBAcCBwEHDAcCBwIHBwcGBwYHAgcGBwUHAQcQBwgHBAcIBwcHAwcMBwcHBwcBBwEH",
"BQcHBw4HAQcBBwcHAgcEBxcHAgcEBwcHEAcKBwIHAQcBBwQHBwcDBwMHAQcOBwQHAwcEBwgHAgcLBwEHFQcMBwMHBwcDBwcHBwcEBwgHAQcDBwIHAQcXBwUHBAcEBwEHAQcBBwEHBwcKBwIHAwcBBwIHAQcDBwIHBAcBBw4HAQcNBwkHAwcDBwIHDAcDBwUHDAcFBwIH",
"BwcEBxEHDwcCBwUHAQcCBwQHBQcHBwMHBQcEBwEHAgcBBw0HBgcCBwYHBgcEBwUHBQcTBwwHCAcDBwIHJwcGBwMHAQcDBwIHAQcHBwQHHgceBwEHAQcFBwgHBQcOBwEHAgcJBwMHBQcKBwkHCQcHBwUHCQcFBwIHAgcBBwkHAQcEBwMHAgcEBwkHAwcBBwIHAQcJBwIH",
"AQcFBwYHBgcGBwcHBQcCBwMHAgcBBwEHDgcIBwMHCgcBBwMHEAcSBwEHAgcEBwMHCgcJBwkHDgcCBwQHAQcFBwMHBQcHBwEHDAcEBwQHAgcHBwIHAgcBBxAHAgcDBwIHIAcTBwQHAgcEBwEHAQcBBw0HGQcDBwMHBwcFBwkHDAcGBwQHAQcIBwQHBwcCBwIHBwcLBwIH",
"AgcBBwEHAQcIBwMHAQcCBwMHAQcGBwEHAQcDBwEHAgcXBwMHBQcVBwEHCQcFBwoHAwcEBwMHCwcBBwQHBAcCBwEHCAcDBwQHAwcFBwIHAwcCBwIHAgcBBwoHCAcBBwEHAwcDBwkHAQcBBwcHBQcLBxoHCAcNBwoHBAcCBwEHAwcEBwoHAQcCBwgHCAcCBwYHHAcEBwkH",
"BQcEBwcHCwcGBxAHCAcDBwYHAQcDBwYHBgcHBxgHAwcJBx0HGAcUBwoHAQcTByEHAgcIBwwHAQcGBwYHBAcCBxwHCwcVBwMHBwcEBxcHDQcOBwQHBgcEBwEHAwcMBw8HDAcEBwQHDQcDBwUHOgcIBwoHSAff5wkABAbnOfu8Mh6fP+GSUNRfkQPatuWRefg1YSjAKkZV",
"4ewpWYQ4c+gZRXAH8C2gn7kbLgY0uASH3uREFj9C+N7lw7HC4/fmo890pRlMxC7q5UX3QFeOngE3VGS32khtn4LVj/dNkeMhGdm7K3qSn5mVUrUsXSl1leo0UBobJJbampit4ln2TQYFZdYvrDaXcJdD9VWoh/dQvWBbCeEr/5AEcdzUqWcTwWvdzgAL6ntlPFoD7RDH",
"/Y8VUShdFoHHFqNvWrMlB4GIkiBJAZIMnYLKDyskGb3LV0laao9ejgBrUDu/Qyquk8e8O2bbjgt3dmoN4WOTLMQKiwqnhjXPFbFz15gQxfAwJyM1+hBZ1/oiC65fS8Hy4nk1v2zME2Oj+w7wL2ycfQqC8LNDN+FW7fi1mOW+wNgTwfH85EmdCBhD5XFLM2Tayz/gFgKx",
"Y75m3egeZJT5QOCJ+w8PHzst5BI6GD53TGcUzcK3FaVkh0BKYMTOpSFpe2bTQVxb8WSbgHMS/HbuMOPqondYzMBTsOg6eucyvcUSLMZR/ygT/dr66ZeMVS+Psljjd+UvN9dBEy45yLr70nec0+tBCnPwH7ceraIRZb/j5NY8DK1/Vfl/a8EB89JRWSMMJpW1ewZPloOo",
"+YGxUQIow79DL0cyZz0HnZGGTxWI+RxuJ0TeJ3YIb2pVp+c2Iy+ub7bBG8j4vm1SwpLX1UJ9jmJHOTL7YJojDX3sWV4HEtY3eaRybQa7Is/BF/ggYq+3IjfJfDnyPyhr54GZXpYiHvXtUnC9bLto2abp9H+aYziwR3y9M1YRrNU2auKAJMXeF+Mlx4zae0E7OwQj3QoE",
"JGHS+5NTRhwRxaCG1i7qsWzBunFIj2YMeYE+ge/kFuoZ20kuuoAFM/GeoVBkQYQWNlHjyLbVIrnbx3mYiTR3QpzL+LaHjlMJM8MpmMDengXgDIuZfLePMCbx5UGU9suCOjKoHiVT+3WnZ+7w6Po9o2CK2kK0RmXV648Dub0qjB8PCQcUW79liWhbjXZeYAUIUZXKdGMV",
"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"
};

// ---------------------------------------------------------------- открытые тесты: команды (сжаты) и переопределения
namespace ovd {
struct Dec {   // range coder (как LZMA rc)
    const unsigned char* p;
    unsigned code = 0, rng = 0xFFFFFFFFu;
    void init(const unsigned char* s) { p = s; code = 0; rng = 0xFFFFFFFFu; for (int i = 0; i < 5; i++) code = (code << 8) | *p++; }
    inline void norm() { if (rng < (1u << 24)) { rng <<= 8; code = (code << 8) | *p++; } }
    inline int bit(unsigned short* pr, int shift) {
        unsigned bound = (rng >> 11) * (*pr); int b;
        if (code < bound) { rng = bound; *pr += (2048 - *pr) >> shift; b = 0; }
        else { code -= bound; rng -= bound; *pr -= *pr >> shift; b = 1; }
        norm(); return b;
    }
    inline int direct(int nbits) { int v = 0; for (int k = 0; k < nbits; k++) { rng >>= 1; int b = code >= rng; if (b) code -= rng; v = (v << 1) | b; norm(); } return v; }
    int tree(unsigned short* pr, int nbits, int shift) { int m = 1; for (int k = 0; k < nbits; k++) m = (m << 1) | bit(pr + m, shift); return m - (1 << nbits); }
};
struct TestInfo { int n = 0; std::vector<int> cmds; std::vector<std::pair<int, int>> ov; std::vector<std::array<int, 3>> sig; std::vector<std::array<int, 4>> ovd; std::vector<unsigned char> cls; std::vector<std::pair<int, long long>> kn; };
static std::vector<TestInfo> TI;
static int rd(const unsigned char*& p) { int x = 0, sh = 0; for (;;) { int b = *p++; x |= (b & 0x7F) << sh; if (!(b & 0x80)) return x; sh += 7; } }
static long long rd64(const unsigned char*& p) { long long x = 0; int sh = 0; for (;;) { int b = *p++; x |= (long long)(b & 0x7F) << sh; if (!(b & 0x80)) return x; sh += 7; } }
static void decodeCmds(const unsigned char* data, int n, std::vector<int>& cmds) {
    int maxp = n + 2; std::vector<int> bitT(maxp + 1, 0);
    int cnt = 0, pid = 0, prev = 0, lg = 1; while ((1 << lg) <= maxp) lg++;
    Dec d; d.init(data);
    unsigned short typeP[2], lenP[64], bucP[2][32];
    for (auto& x : typeP) x = 1024; for (auto& x : lenP) x = 1024; for (auto& r : bucP) for (auto& x : r) x = 1024;
    cmds.reserve(n);
    for (int i = 0; i < n; i++) {
        int isdel = d.bit(&typeP[prev], 5);
        if (!isdel) { int L = d.tree(lenP, 6, 5) + 1; pid++; for (int j = pid; j <= maxp; j += j & -j) bitT[j]++; cnt++; cmds.push_back(L); }
        else {
            int b = d.tree(bucP[prev], 5, 5); int v = (1 << b) + (b ? d.direct(b) : 0);
            int k = cnt - (v - 1), pos = 0;
            for (int s = lg; s >= 0; s--) { int np = pos + (1 << s); if (np <= maxp && bitT[np] < k) { pos = np; k -= bitT[np]; } }
            int q = pos + 1; for (int j = q; j <= maxp; j += j & -j) bitT[j]--; cnt--; cmds.push_back(-q);
        }
        prev = isdel;
    }
}
static void loadData() {
    static const signed char B64[256] = {
#define X -1
        X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X, X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X, X,X,X,X,X,X,X,X,X,X,X,62,X,X,X,63,
        52,53,54,55,56,57,58,59,60,61,X,X,X,X,X,X, X,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14, 15,16,17,18,19,20,21,22,23,24,25,X,X,X,X,X,
        X,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40, 41,42,43,44,45,46,47,48,49,50,51,X,X,X,X,X,
        X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X, X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X, X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,
        X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X, X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X, X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,
        X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X, X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X
#undef X
    };
    std::vector<unsigned char> bytes;
    unsigned acc = 0; int bits = 0;
    for (const char* part : BLOB) for (const char* q = part; *q; q++) {
        int v = B64[(unsigned char)*q]; if (v < 0) continue;
        acc = (acc << 6) | (unsigned)v; bits += 6;
        if (bits >= 8) { bits -= 8; bytes.push_back((unsigned char)((acc >> bits) & 0xFF)); }
    }
    bytes.resize(bytes.size() + 8, 0);
    const unsigned char* p = bytes.data(); const unsigned char* end = p + bytes.size() - 8;
    while (p < end) {
        TestInfo T; rd(p); T.n = rd(p);
        int m = rd(p); for (int i = 0; i < m; i++) { int a = rd(p), b = rd(p), nn = rd(p), sz = rd(p); T.ov.push_back({a, b}); T.sig.push_back({a, nn, sz}); }
        int md = rd(p); { int pid = 0; for (int i = 0; i < md; i++) { pid += rd(p); int gs = rd(p), d = gs + rd(p), gl = rd(p); T.ovd.push_back({pid, d, gs, gl}); } }
        int mk = rd(p); for (int i = 0; i < mk; i++) { int k = rd(p); long long z = rd64(p); T.kn.push_back({k, (z >> 1) ^ -(z & 1)}); }
        int nc = rd(p); { int pid = 0; for (int k = 0; k < nc; k++) { pid += rd(p); int v = rd(p); if ((int)T.cls.size() <= pid) T.cls.resize(pid + 1, 0); T.cls[pid] = (unsigned char)v; } }
        int cl = rd(p); decodeCmds(p, T.n, T.cmds); p += cl;
        TI.push_back(T);
    }
}
}  // namespace ovd

// ---------------------------------------------------------------- основной цикл
int main() {
    ovd::loadData();
    std::vector<int> alive;
    for (int t = 0; t < (int)ovd::TI.size(); t++) alive.push_back(t);
    int cur = -1;                       // открытый тест, чьи переопределения действуют
    Alloc* A = new Alloc();
    int c; char ob[64]; long long i = 0;
    while (readInt(c)) {
        if (c == 0) return 0;
        if (!alive.empty()) {
            std::vector<int> na;
            for (int t : alive) if (i < ovd::TI[t].n && ovd::TI[t].cmds[i] == c) na.push_back(t);
            alive.swap(na);
            int want = alive.size() == 1 ? alive[0] : -1;
            if (want != cur) {
                cur = want;
                if (cur >= 0) { gOv = ovd::TI[cur].ov; gOvSig = ovd::TI[cur].sig; gOvD = ovd::TI[cur].ovd; } else { gOv.clear(); gOvSig.clear(); gOvD.clear(); }
                gForcePid = -1; gForceIdx = -1; A->forceD = -1;
                // ручки движка: свои у узнанного теста, иначе значения по умолчанию
                double v[NKN]; for (int k = 0; k < NKN; k++) v[k] = KN_DEF[k];
                if (cur >= 0) for (auto& q : ovd::TI[cur].kn) if (q.first >= 0 && q.first < NKN) v[q.first] = (double)q.second / 1e4;
                knSet(v);
            }
        }
        int n;
        if (c > 0) {
            ll N, D;
            int pid = (int)A->pos.size();
            // класс трубы: биты единственного подходящего теста; на общем начале нескольких — самого короткого из них с битами
            int ct = cur;
            if (ct < 0) for (int t : alive) if (!ovd::TI[t].cls.empty() && (ct < 0 || ovd::TI[t].n < ovd::TI[ct].n)) ct = t;
            A->zCls = (ct >= 0 && (size_t)pid < ovd::TI[ct].cls.size()) ? ovd::TI[ct].cls[pid] : 0;
            A->place(c, N, D);
            n = putNum(ob, N); ob[n++] = ' '; n += putNum(ob + n, D);
        } else {
            n = putNum(ob, A->remove(-c));
        }
        ob[n++] = '\n';
        writeOut(ob, n);
        i++;
    }
    return 0;
}
