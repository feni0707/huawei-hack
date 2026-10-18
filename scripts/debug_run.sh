#!/usr/bin/env bash
# Один тест БЕЗ runner.py, с видимым stderr интерактора (runner.py его глушит в DEVNULL).
# Показывает "Score = occupied/trench * weight = ..." или точную причину Wrong Answer
# (коллизия труб, выход за траншею, неверная позиция при удалении, ошибка формата).
#   scripts/debug_run.sh <solution_binary> <test.in>
# Транскрипт ответов решения (то, что прочитал интерактор): build/debug/tout.txt
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
[[ $# -eq 2 ]] || { echo "usage: $0 <solution_binary> <test.in>" >&2; exit 2; }
SOL="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
TEST="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"
INTERACTOR="$("$ROOT/scripts/build_interactor.sh")" || exit 1
WORK="$ROOT/build/debug"; mkdir -p "$WORK"; cd "$WORK"
rm -f fifo score.txt; mkfifo fifo
"$INTERACTOR" "$TEST" tout.txt <fifo 2>interactor.log | "$SOL" >fifo
cat interactor.log
rm -f fifo
