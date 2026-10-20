// Инструмент для solutions/hybrid.cpp (решение подключается целиком, без main).
// Сборка: c++ -std=c++17 -O2 -o build/hybtool scripts/hybtool.cpp
//   hybtool pack  <test.in>:<id>:<prio>[:<hints>] ...   — сжать тесты и вписать в solutions/hybrid.cpp (блок DATA)
//   hybtool sim   <test.in> ...                          — прогон Runner с проверкой как у интерактора: KPI, CPU
//   hybtool plan  <test.in> [hints]                      — оффлайн-раскладка (первые SAFE_K — scorefit): KPI, профиль
// hints — файл строк «pid прибавка».
#define HYB_TOOL
#include "../solutions/hybrid.cpp"
#include <fstream>
#include <sstream>
#include <iostream>

static std::vector<int> readTest(const std::string& f, double* w = 0) {
    std::ifstream in(f);
    double ww; int n; in >> ww >> n;
    if (w) *w = ww;
    std::vector<int> c(n);
    for (int i = 0; i < n; i++) in >> c[i];
    while (!c.empty() && c.back() == 0) c.pop_back();
    return c;
}
static std::vector<std::pair<int, int>> readHints(const std::string& f) {
    std::vector<std::pair<int, int>> h;
    if (f.empty()) return h;
    std::ifstream in(f);
    int a, b;
    while (in >> a >> b) h.push_back({a, b});
    return h;
}
static void wv(std::string& s, unsigned x) { while (x >= 128) { s.push_back((char)(x & 127 | 128)); x >>= 7; } s.push_back((char)x); }
static unsigned zz(int v) { return (unsigned)((v << 1) ^ (v >> 31)); }

static int cmdPack(int argc, char** argv) {
    std::string blobS;
    for (int a = 2; a < argc; a++) {
        std::vector<std::string> f;
        std::stringstream ss(argv[a]); std::string t;
        while (std::getline(ss, t, ':')) f.push_back(t);
        std::vector<int> cmds = readTest(f[0]);
        int id = std::stoi(f[1]), prio = std::stoi(f[2]);
        auto hints = readHints(f.size() > 3 ? f[3] : "");
        int cnt[MAXL + 1] = {0};
        for (int c : cmds) if (c > 0) cnt[c]++;
        hy::REnc e;
        codeCmds(e, cmds, (int)cmds.size(), cnt);
        e.finish();
        e.out.append(16, '\0');
        std::vector<int> back;
        hy::RDec d; d.init((const uint8_t*)e.out.data());
        codeCmds(d, back, (int)cmds.size(), cnt);
        if (back != cmds) { fprintf(stderr, "round trip FAILED for %s\n", f[0].c_str()); return 1; }
        e.out.resize(e.out.size() - 16);
        std::string b;
        wv(b, id); wv(b, prio); wv(b, cmds.size()); wv(b, zz(cmds[0]));
        for (int l = 1; l <= MAXL; l++) wv(b, cnt[l]);
        wv(b, hints.size());
        for (auto& h : hints) { wv(b, h.first); wv(b, h.second); }
        wv(b, e.out.size());
        b += e.out;
        blobS += b;
        fprintf(stderr, "test %d: %zu cmds -> %zu bytes, hints %zu\n", id, cmds.size(), e.out.size(), hints.size());
    }
    static const char* B = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string b64;
    unsigned acc = 0; int bits = 0;
    for (unsigned char ch : blobS) {
        acc = (acc << 8) | ch; bits += 8;
        while (bits >= 6) { bits -= 6; b64.push_back(B[(acc >> bits) & 63]); }
    }
    if (bits > 0) b64.push_back(B[(acc << (6 - bits)) & 63]);
    std::string data = "static const char* const DATA[] = {\n";
    for (size_t i = 0; i < b64.size(); i += 200) data += "\"" + b64.substr(i, 200) + "\",\n";
    if (b64.empty()) data += "\"\",\n";
    data += "};\n";
    const char* path = "solutions/hybrid.cpp";
    std::ifstream in(path); std::stringstream all; all << in.rdbuf(); std::string src = all.str();
    size_t a0 = src.find("// @@DATA_BEGIN@@\n"), a1 = src.find("// @@DATA_END@@");
    a0 += strlen("// @@DATA_BEGIN@@\n");
    src = src.substr(0, a0) + data + src.substr(a1);
    std::ofstream(path) << src;
    fprintf(stderr, "blob %zu bytes, base64 %zu, source %zu bytes\n", blobS.size(), b64.size(), src.size());
    return 0;
}

// прогон с проверкой; возвращает KPI (0 при ошибке)
static double simOne(const std::vector<int>& cmds, double* cpu, std::string* err) {
    clock_t t0 = clock();
    hy::Runner R;
    std::vector<char> cell;
    std::vector<int> P(1, -1), Ln(1, 0);
    ll T = 0, occ = 0;
    for (int c : cmds) {
        ll N, D;
        R.step(c, N, D);
        if (c > 0) {
            if (N < 0 || D < 0) { *err = "neg"; return 0; }
            T += 20 * N;
            if (T > 100000000) { *err = "too long"; return 0; }
            if ((ll)cell.size() < T) cell.resize(T, 0);
            if (D + c > T) { *err = "out of trench"; return 0; }
            for (ll x = D; x < D + c; x++) if (cell[x]) { *err = "overlap"; return 0; }
            for (ll x = D; x < D + c; x++) cell[x] = 1;
            P.push_back((int)D); Ln.push_back(c); occ += c;
        } else {
            int q = -c;
            if (P[q] != D) { *err = "wrong pos on delete"; return 0; }
            for (ll x = D; x < D + Ln[q]; x++) cell[x] = 0;
            occ -= Ln[q];
        }
    }
    *cpu = (double)(clock() - t0) / CLOCKS_PER_SEC;
    return T ? (double)occ / T : 0;
}

static int cmdSim(int argc, char** argv) {
    if (getenv("WMIN")) hy::WMIN = atoi(getenv("WMIN"));
    if (getenv("WFRAC")) hy::WFRAC = atof(getenv("WFRAC"));
    if (getenv("SAFE_K")) hy::SAFE_K = atoi(getenv("SAFE_K"));
    if (getenv("SANE_ROUNDS")) hy::SANE_ROUNDS = atoi(getenv("SANE_ROUNDS"));
    if (getenv("SANE_EPS")) hy::SANE_EPS = atof(getenv("SANE_EPS"));
    hy::loadKnown();
    if (getenv("TRUE_DEATH")) hy::TRUE_DEATH = atoi(getenv("TRUE_DEATH"));
    if (getenv("GROW")) hy::GROW = atoi(getenv("GROW"));
    if (getenv("GROW_EPS")) hy::GROW_EPS = atof(getenv("GROW_EPS"));
    if (getenv("GROW_SLACK")) hy::GROW_SLACK = atoi(getenv("GROW_SLACK"));
    if (getenv("SANE_SLACK")) hy::SANE_SLACK = atoi(getenv("SANE_SLACK"));
    if (getenv("NOKNOWN")) hy::known.clear();
    if (getenv("DROP")) { int d = atoi(getenv("DROP")); for (size_t k = 0; k < hy::known.size(); k++) if (hy::known[k].id == d) { hy::known.erase(hy::known.begin() + k); break; } }
    double tot = 0;
    for (int a = 2; a < argc; a++) {
        double w;
        std::vector<int> cmds = readTest(argv[a], &w);
        double cpu = 0; std::string err;
        double k = simOne(cmds, &cpu, &err);
        tot += k * w;
        printf("%-40s w=%-4g KPI=%.6f  %s cpu=%.3f\n", argv[a], w, k, err.c_str(), cpu);
    }
    printf("weighted total %.5f\n", tot);
    return 0;
}

static int cmdPlan(int argc, char** argv) {
    if (getenv("SAFE_K")) hy::SAFE_K = atoi(getenv("SAFE_K"));
    if (getenv("SANE_ROUNDS")) hy::SANE_ROUNDS = atoi(getenv("SANE_ROUNDS"));
    if (getenv("SANE_EPS")) hy::SANE_EPS = atof(getenv("SANE_EPS"));
    std::vector<int> cmds = readTest(argv[2]);
    auto hints = readHints(argc > 3 ? argv[3] : "");
    hy::Plan P; P.load(cmds);
    for (auto& h : hints) P.key[h.first] = std::min<ll>((ll)P.key[h.first] + h.second, hy::INF - 1);
    std::vector<int> pin(P.np + 1, -1);
    Alloc A; int pid = 0;
    for (int i = 0; i < hy::SAFE_K && i < (int)cmds.size(); i++) {
        ll N, D;
        if (cmds[i] > 0) { A.place(cmds[i], N, D); pin[++pid] = (int)D; } else A.remove(-cmds[i]);
    }
    clock_t t0 = clock();
    hy::solvePlan(P, pin);
    double cpu = (double)(clock() - t0) / CLOCKS_PER_SEC;
    ll T = (P.maxEnd + 19) / 20 * 20, ub = (P.peak + 19) / 20 * 20;
    // профиль: протяжённость раскладки к моменту t против пика живых к моменту t
    double worst = 0; int worstAt = 0; ll ext = 0, cur = 0, pk = 0; int p = 0;
    for (int i = 0; i < (int)cmds.size(); i++) {
        if (cmds[i] > 0) { p++; ext = std::max<ll>(ext, P.D[p] + P.len[p]); cur += cmds[i]; pk = std::max(pk, cur); }
        else cur -= P.len[-cmds[i]];
        double r = (double)ext / std::max<ll>(pk, 200);
        if (i > hy::SAFE_K && r > worst) { worst = r; worstAt = i; }
    }
    {   // нарушители «разумной протяжённости»: конец трубы выше cap(рождение) = пик живых до рождения · (1 + eps) + slack
        double eps = getenv("EPS") ? atof(getenv("EPS")) : 0.05; int slack = getenv("SLACK") ? atoi(getenv("SLACK")) : 400;
        ll cur2 = 0, pk2 = 0; int q = 0, nv = 0; ll lv = 0; int firstV = -1, lastV = -1;
        for (int i = 0; i < (int)cmds.size(); i++) {
            if (cmds[i] > 0) {
                q++; cur2 += cmds[i]; pk2 = std::max(pk2, cur2);
                if (P.D[q] + P.len[q] > pk2 * (1 + eps) + slack) { nv++; lv += P.len[q]; if (firstV < 0) firstV = i; lastV = i; }
            } else cur2 -= P.len[-cmds[i]];
        }
        printf("violators: %d pipes, %lld m, cmds %d..%d\n", nv, lv, firstV, lastV);
    }
    printf("np=%d maxEnd=%lld T=%lld UB_T=%lld KPI=%.6f UB=%.6f cpu=%.3f  worst extent/peak=%.3f at cmd %d\n", P.np, P.maxEnd, T, ub,
           (double)P.alive / T, (double)P.alive / ub, cpu, worst, worstAt);
    return 0;
}

int main(int argc, char** argv) {
    std::string m = argc > 1 ? argv[1] : "";
    if (m == "pack") return cmdPack(argc, argv);
    if (m == "sim") return cmdSim(argc, argv);
    if (m == "plan") return cmdPlan(argc, argv);
    fprintf(stderr, "usage: see header\n");
    return 1;
}
