// Точность предсказателя времени жизни: для каждого горизонта j — доля верных классов «смерть за H_j» (после разгона).
//   exp_pred <test.in>...
#define SIM
#include "../solutions/portfolio.cpp"

#include <fstream>
#include <iostream>

int main(int argc, char** argv) {
    for (int ai = 1; ai < argc; ai++) {
        std::ifstream f(argv[ai]);
        double w; int n;
        f >> w >> n;
        std::vector<int> cmds(n);
        for (auto& c : cmds) f >> c;
        int np = 0;
        for (int c : cmds) if (c > 0) np++;
        std::vector<int> death(np + 2, 1 << 28), birthT(np + 2, 0);
        for (int i = 0; i < n; i++) if (cmds[i] < 0) death[-cmds[i]] = i;
        Predictor P;
        long long tp[Predictor::NH] = {0}, fp[Predictor::NH] = {0}, fn[Predictor::NH] = {0}, tn[Predictor::NH] = {0};
        int pid = 0;
        for (int i = 0; i < n; i++) {
            int c = cmds[i];
            if (c > 0) {
                pid++;
                float p[Predictor::NH];
                P.onPlace(c, p);
                int life = death[pid] - i;
                for (int j = 0; j < Predictor::NH; j++) {
                    bool pred = p[j] >= 0.5f, truth = life <= Predictor::H(j);
                    if (i < n / 5) continue;  // разгон
                    if (pred && truth) tp[j]++; else if (pred) fp[j]++; else if (truth) fn[j]++; else tn[j]++;
                }
            } else {
                P.onDelete(-c);
            }
            P.tick();
        }
        printf("%s\n", argv[ai]);
        for (int j = 0; j < Predictor::NH; j++) {
            long long tot = tp[j] + fp[j] + fn[j] + tn[j];
            if (!tot) continue;
            printf("  H=%-6d acc=%.3f  T-truth=%.3f  prec=%.3f rec=%.3f\n", Predictor::H(j), (double)(tp[j] + tn[j]) / tot,
                   (double)(tp[j] + fn[j]) / tot, tp[j] + fp[j] ? (double)tp[j] / (tp[j] + fp[j]) : 0.0,
                   tp[j] + fn[j] ? (double)tp[j] / (tp[j] + fn[j]) : 0.0);
        }
    }
}
