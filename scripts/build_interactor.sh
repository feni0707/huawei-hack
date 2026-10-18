#!/usr/bin/env bash
# Собирает интерактор из исходников под текущую ОС в build/interactor.
# Нужен на macOS (materials/interactor — Linux ELF), на Linux тоже работает.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/build/interactor"
mkdir -p "$ROOT/build"
if [[ ! -x "$OUT" || "$ROOT/materials/interactor.cpp" -nt "$OUT" ]]; then
  "${CXX:-c++}" -std=c++17 -O2 -I"$ROOT/materials" -o "$OUT" "$ROOT/materials/interactor.cpp"
fi
echo "$OUT"
