// Baseline: каждая труба получает свои новые 20-метровые сегменты в конце траншеи.
// Нужен только как эталон протокола ввода-вывода. Score на открытых тестах ≈ 4.91 из 32.1.
#include <cstdio>
#include <vector>

int main() {
    long long trench = 0;                 // текущая длина траншеи, м
    std::vector<long long> pos(1);        // pos[pid] — где лежит труба, pid с 1
    int c;
    while (scanf("%d", &c) == 1) {
        if (c == 0) return 0;             // конец — ничего не печатаем
        if (c > 0) {
            long long n = (c + 19) / 20;  // сколько 20-метровых сегментов докопать
            long long d = trench;
            trench += n * 20;
            pos.push_back(d);
            printf("%lld %lld\n", n, d);
        } else {
            printf("%lld\n", pos[-c]);    // на удаление — позиция трубы
        }
        fflush(stdout);                   // обязателен после каждого ответа
    }
}
