#!/usr/bin/env bash
# Прогон решения на открытых тестах через официальный materials/runner.py.
#   scripts/run_tests.sh <solution_binary>             — все тесты из materials/tests
#   scripts/run_tests.sh <solution_binary> <test.in>   — один тест
# Решение должно быть ОДНИМ исполняемым файлом (runner не поддерживает "python3 sol.py" на POSIX —
# оберните в shell-скрипт с shebang). Вывод runner'а: weighted score по каждому тесту и Total.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
[[ $# -ge 1 ]] || { echo "usage: $0 <solution_binary> [test.in]" >&2; exit 2; }
SOL="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
INTERACTOR="$("$ROOT/scripts/build_interactor.sh")"
TEST=""; [[ $# -ge 2 ]] && TEST="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"
WORK="$ROOT/build/run"; mkdir -p "$WORK"; cd "$WORK"   # runner пишет score.txt/answer.txt/report.txt в cwd
if [[ -n "$TEST" ]]; then
  python3 "$ROOT/materials/runner.py" --interactor "$INTERACTOR" -i "$TEST" "$SOL"
else
  OUT="$(python3 "$ROOT/materials/runner.py" --interactor "$INTERACTOR" -d "$ROOT/materials/tests" "$SOL")"
  grep '^Test' <<<"$OUT" | sort -b -k2,2n   # runner перечисляет тесты в порядке os.listdir
  grep -v '^Test' <<<"$OUT" || true
fi
