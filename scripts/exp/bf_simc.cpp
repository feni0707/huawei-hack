// Адаптер bestfit.cpp под стенд simc (только для стенда).
#include "../../solutions/bestfit.cpp"
struct Alloc : BestFit {
    struct Z { long long get(long long) const { return 0; } } zl;
    int tailStart() const { return (int)tail; }
};
