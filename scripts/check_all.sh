#!/usr/bin/env bash
# Полная проверка решения portfolio (Windows + MSVC): сборка стенда и решения, корректность (--check) на всех наборах,
# прогон настоящего интерактора на открытых и граничных тестах, замер CPU/памяти.
#   scripts/check_all.sh            # быстро: стенд + интерактор на open и fuzz
set -e
cd "$(dirname "$0")/.."
cmd //c "$(cygpath -w "$PWD/scripts/build_sim.bat")" 2>&1 | grep -Ei "error|ошибка" && exit 1 || true
cmd //c "$(cygpath -w "$PWD/scripts/build_sol.bat") portfolio" 2>&1 | grep -Ei "error|ошибка" && exit 1 || true
[ -d build/fuzz ] || python scripts/gen_fuzz.py
[ -d build/fuzz2 ] || python scripts/gen_rand_suite.py 6
echo "== стенд: корректность укладок (--check) на всех наборах"
for d in build/fuzz build/fuzz2/*; do ./build/bin/sim.exe --port --check "$d" | tail -1; done
echo "== стенд: сводка по классам (suite5)"
scripts/suite5.sh --port --check --line
echo "== настоящий интерактор: open + fuzz (CPU/память на Windows включают ввод-вывод ≈2.1 с на 200k команд)"
python scripts/win_bench.py build/bin/portfolio.exe materials/tests build/fuzz -j 1
