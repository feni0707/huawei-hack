#!/usr/bin/env bash
# Стандартный набор для оценки: open + gen2/gen по seed 1..3. Использование: scripts/suite.sh [--pol K | --port] [доп. опции sim]
cd "$(dirname "$0")/.."
SETS="materials/tests"
for s in 1 2 3; do
  SETS="$SETS build/gen2_s$s/twophase build/gen2_s$s/sorted2 build/gen2_s$s/invert build/gen2_s$s/tiny build/gen_s$s/alternating build/gen_s$s/periodic build/gen_s$s/phases"
done
./build/bin/sim.exe "$@" $SETS
