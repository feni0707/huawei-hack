#!/usr/bin/env bash
# Расширенный набор №3: suite2 + sortedbig (gen3 seed 1..2). Аргументы — как у sim (без путей). Бинарник: $SIM (по умолчанию build/bin/sim.exe)
cd "$(dirname "$0")/.."
SETS="materials/tests"
for s in 1 2 3; do
  SETS="$SETS build/gen2_s$s/twophase build/gen2_s$s/sorted2 build/gen2_s$s/invert build/gen2_s$s/tiny build/gen_s$s/alternating build/gen_s$s/periodic build/gen_s$s/phases"
done
for s in 1 2; do for c in sawtooth bimodal fifo lifo desc sortedbig; do SETS="$SETS build/gen3_s$s/$c"; done; done
${SIM:-./build/bin/sim.exe} "$@" $SETS
